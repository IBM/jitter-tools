#define _GNU_SOURCE
#include <mpi.h>
#include <pthread.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <errno.h>
#include <string.h>
#include <numa.h>

#define MAX_CPUS_PER_NODE 176

#define POWER9 0x4e  // witherspoon

/* Processor Version Register */
#define SPRN_PVR 0x11F  

/* Version field */
#define PVR_VER(pvr)    (((pvr) >>  16) & 0xFFFF)       

int get_pvr(void)
{
   unsigned long ret;
   int msr, pvr;
   asm volatile("mfspr %0,%1" : "=r" (ret): "i" (SPRN_PVR));
   msr = (int) ret;
   pvr = PVR_VER(msr);
   return pvr;
}

//===========================================================================
// This routine will bind threads only if the env variable BIND_THREADS=yes.
// This version uses a simple cpulist : 0,1,2,..., max_cpus - 1.
// default : spread out MPI ranks and threads evenly
// options : BIND_STRIDE=number   ... stride per process not per thread
//           BIND_SKIP=number of logical cpus to skip
//           BIND_CPU_LIST="cpu1,cpu2,cpu3,...,cpuN" a specific list
//===========================================================================

#pragma bindthreads_=bindthreads
void bindthreads(int *ppn, int *lrank)
{
  char * ptr;
  int envlist[MAX_CPUS_PER_NODE];
  int uselist[MAX_CPUS_PER_NODE];
  int socket0_cpus[MAX_CPUS_PER_NODE];
  int socket1_cpus[MAX_CPUS_PER_NODE];

  int cpu, max_cpus, inc, ndx, socket;
  int i, j, k, myrank, nranks, tid, skip;
  int max_ranks_per_node, my_base_ndx, cpus_per_rank, verbose;
  int nthreads, rc, use_envlist, badlist, processor_version, smt_width;
  int available_cpus_per_socket, cpus_per_socket;
  int system_cores, available_cores_per_socket, socket1_base_cpu;
  int half_smt_width, num_halfs, cpu_inc;
  char * list_ptr;
  char delimiters[] = {", "};
  pthread_t thread;
  cpu_set_t cpuset;
  struct bitmask * cpus;
  char * snames, * rnames, host[80];
  int bind_threads, match, color, local_rank;
  int ranks_per_node, ranks_per_socket;
  MPI_Comm local_comm;

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);

  bind_threads = 0;
  ptr = getenv("BIND_THREADS");
  if (ptr != NULL ) {
    if (strncasecmp(ptr,"yes",3) == 0) bind_threads = 1;
  }

  verbose = 0;
  ptr = getenv("BIND_VERBOSE");
  if (ptr != NULL) {
     if (strncasecmp(ptr,"yes",3) == 0) verbose = 1;
  }

  ptr = getenv("SYSTEM_CORES");
  if (ptr != NULL) system_cores = atoi(ptr);
  else             system_cores = 2;

  // make a communicator of all ranks on this node
  snames = (char *) malloc(nranks*sizeof(host));
  rnames = (char *) malloc(nranks*sizeof(host));
  gethostname(host, sizeof(host));

  for (i=0; i<sizeof(host); i++) {    
     if (host[i] == '.') {
         host[i] = '\0';
         break;
     }
  }

  for (i=0; i<nranks; i++)  {
    ptr = snames + i*sizeof(host);
    strncpy(ptr, host, sizeof(host));
  }

  MPI_Alltoall(snames, sizeof(host), MPI_BYTE,
               rnames, sizeof(host), MPI_BYTE, MPI_COMM_WORLD);
  color = 0;
  match = 0;
  for (i=0; i<nranks; i++) {
    ptr = rnames + i*sizeof(host);
    if (strcmp(host, ptr) == 0) {    
      match++;
      if (match == 1) color = i;
    }
  }

  MPI_Comm_split(MPI_COMM_WORLD, color, myrank, &local_comm);
  MPI_Comm_rank(local_comm, &local_rank);
  MPI_Comm_size(local_comm, &ranks_per_node);

  *ppn = ranks_per_node;
  *lrank = local_rank;

  // do not return until local_rank and ranks_per_node are set
  if (bind_threads == 0) {
    if (myrank==0) fprintf(stderr,"bindthreads: not binding because BIND_THREADS is not set to yes.\n");
    return;
  }
 
  ranks_per_socket = ranks_per_node / 2;
  socket = local_rank / ranks_per_socket;

  MPI_Allreduce(&ranks_per_node, &max_ranks_per_node, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  // max_cpus = max #cpus on the system
  max_cpus = sysconf(_SC_NPROCESSORS_ONLN);

  processor_version = get_pvr();

  if (processor_version != POWER9) {
    if (myrank==0) fprintf(stderr,"bindthreads: not binding because the processor is not power9.\n");
  }

  j = 0;
  // Power9 socket 0 is numa node 0
  cpus = numa_allocate_cpumask();
  rc = numa_node_to_cpus(0, cpus);
  if (rc >= 0) {
    for (i = 0; i < cpus->size; i++)
      if (numa_bitmask_isbitset(cpus, i)) socket0_cpus[j++] = i;
  }
  
  j = 0;
  // Power9 socket 1 is numa node 8
  cpus = numa_allocate_cpumask();
  rc = numa_node_to_cpus(8, cpus);
  if (rc >= 0) {
    for (i = 0; i < cpus->size; i++)
      if (numa_bitmask_isbitset(cpus, i)) socket1_cpus[j++] = i;
  }

  cpus_per_socket = j;

  socket1_base_cpu = socket1_cpus[0];

  k = 0; 
  for (i = 0; i < cpus_per_socket; i++) {
     if ( socket0_cpus[i] < (socket1_base_cpu - 2*system_cores) ) k++;
  }

  available_cpus_per_socket = k;

  // find max available OpenMP threads
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#else
  nthreads = 1;
#endif

 cpus_per_rank = (2*available_cpus_per_socket) / max_ranks_per_node;

  // find max available OpenMP threads
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#else
  nthreads = 1;
#endif

  use_envlist = 0;
  list_ptr = getenv("BIND_CPU_LIST");
  if (list_ptr != NULL) {
     k = 0;
     ptr = strtok(list_ptr, delimiters);
     while(ptr != NULL) {
        envlist[k] = atoi(ptr);
        ptr = strtok(NULL, delimiters);
        k++;
     }

     if (k != max_ranks_per_node*nthreads) {
        if (myrank == 0) fprintf(stderr, "bindthreads: BIND_CPU_LIST requires ranks_per_node*nthreads items ... using defaults\n");
        use_envlist = 0;
     }
     else  {
        badlist = 0;
        for (k=0; k<max_ranks_per_node*nthreads; k++) {
           if (envlist[k] < 0 || envlist[k] >= max_cpus) badlist = 1;
        }
        if (badlist) {
           if (myrank == 0) fprintf(stderr, "bindthreads: BIND_CPU_LIST must have 0<=cpu<max_cpus ... using defaults\n");
           use_envlist = 0;
        }
        else use_envlist = 1;
     }

     if (use_envlist == 1) {
        k = 0;
        for (j=0; j<nthreads; j++) {
           ndx = j + local_rank*nthreads;
           uselist[k] = envlist[ndx];
           k++;
        }
        goto setaffinity;
     }
  }
  else {
     ptr = getenv("BIND_STRIDE");
     if (ptr != NULL) cpus_per_rank = atoi(ptr);

     ptr = getenv("BIND_SKIP");
     if (ptr != NULL) {
        skip = atoi(ptr);
     }
     else skip = 0;

     inc = cpus_per_rank / nthreads;

     if (socket == 0) {
        my_base_ndx = skip + cpus_per_rank * local_rank;
        k = 0;
        for (i=0; i<nthreads; i++) {
           ndx = my_base_ndx + i*inc;
           uselist[k] = socket0_cpus[ndx];
           k++;
        }
     }
     else {
        my_base_ndx = skip + cpus_per_rank * (local_rank - ranks_per_socket);
        k = 0;
        for (i=0; i<nthreads; i++) {
           ndx = my_base_ndx + i*inc;
           uselist[k] = socket1_cpus[ndx];
           k++;
        }
     }
  }

setaffinity:

  // set affinity for OpenMP threads
#ifdef _OPENMP
#pragma omp parallel private(tid,thread,cpu,cpuset,rc)
{
  #pragma omp critical
  {
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif
    thread = pthread_self();
    cpu = uselist[tid];
    CPU_ZERO(&cpuset);
    CPU_SET(cpu, &cpuset);
    rc = pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
    if (rc == 0 && verbose) {
       fprintf(stderr, "bindthreads: binding MPI rank %d thread %d to cpu %d\n", myrank, tid, cpu);
    }
    if (rc != 0) {
      fprintf(stderr, "bindthreads: pthread_set_affinity_np failed for MPI rank %d, thread %d, cpu %d on host %s \n", myrank, tid, cpu, host);
    }
#ifdef _OPENMP
  }
}
#endif

  return;
}
