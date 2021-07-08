#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <pthread.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <errno.h>
#include <string.h>

#define BARRIER  10
#define EXCHANGE 11

#define SORT_ASCENDING_ORDER   1
#define SORT_DESCENDING_ORDER -1

#define POWER9 0x4e  // witherspoon
#define POWER8 0x4c  // minsky

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

void bindthreads(void);

void compute(long, int, double *, double *);

double lutexp(double);

void exch(int);

void print_help(void);

void sortx(double *, int, int *, int);

void TMrand(unsigned long seed, int npts, double * xrand);

static int myrank, myrow, mycol, npex, npey, exchcount;

static int local_rank, ranks_per_node;

static double * sendn,  * sende,  * sends,  * sendw,  * recvn,  * recve,  * recvs,  * recvw;

static MPI_Comm local_comm;


int main(int argc, char * argv[])
{
  int i, k, exchbytes, nranks, ch, help_flag = 0;
  int iter, maxiter, calib_iters, bin, numbins, method;
  long npts;
  double  ssum, data_volume, bw;
  double  t1, t2, t3, tmin, tmax, delay;
  double * xrand;
  unsigned long seed;
  int nrand = 1000;

  double * tcomm, * tcomp, * tstep;
  double compute_interval_msec;
  double elapsed1, elapsed2;

  double histo_bin_width, xhisto, hmin, prob;
  double tavg, ssq, samples, relative_variation;
  double dev, sigma, scb, sp4, skewness, kurtosis;
  int * histo;

  char host[80], * ptr, * snames, * rnames;
  int * sort_key;
  int num_nodes, color, key;
  double * compmax, compmin, compavg, compssq;
  double sigma_comp, samples_total, ssq_comp;
  double * aggregate_sigma_comp, * aggregate_compavg;
  MPI_Comm node_comm;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);
  bindthreads();

  if (myrank == 0) printf("checking program args ...\n");

  // set sensible default values then check program args
  numbins = 50;
  maxiter = 1000;
  compute_interval_msec = 50.0;
  exchbytes = 10000;
  method = EXCHANGE;

  while (1) {
     ch = getopt(argc, argv, "hm:t:i:n:x:");
     if (ch == -1) break;
     switch (ch) {
        case 'h':
           help_flag = 1;
           break;
        case 't':  // set the time for one compute interval in msec
           compute_interval_msec = (double) atoi(optarg);
           break;
        case 'm':  // choose barrier or exchange
           if (0 == strncasecmp(optarg, "barrier", 7)) method = BARRIER;
           else                                        method = EXCHANGE;
           break;
        case 'i':  // set the number of iterations of {compute, communicate}
           maxiter = atoi(optarg);
           break;
        case 'n':  // set the number of histogram bins
           numbins = atoi(optarg);
           break;
        case 'x':  // set the number of bytes for neighbor exchange
           exchbytes = atoi(optarg);
           break;
        default:
           break;
     }   
  }

  if (help_flag) {
     if (myrank == 0) print_help();
     MPI_Finalize();
     return 0;
  }

  if (myrank == 0) {
     if (method == EXCHANGE) {
        printf("using neighbor exchange : compute_interval_msec = %.2lf, maxiter = %d, numbins = %d, exchbytes = %d\n",
                compute_interval_msec, maxiter, numbins, exchbytes);
     }
     else {
        printf("using global barrier : compute_interval_msec = %.2lf, maxiter = %d, numbins = %d\n",
                compute_interval_msec, maxiter, numbins);
     }
  }

  npex = (int) (0.1 + sqrt((double) nranks));
  while (nranks % npex != 0) npex++;
  npey = nranks / npex;

  exchcount = exchbytes / sizeof(double);

  mycol = myrank % npex;
  myrow = myrank / npex;

  if (myrank == 0) printf("starting with nranks = %d\n", nranks);

  if (myrank == 0 && method == EXCHANGE) printf("using npex = %d, npey = %d\n", npex, npey);

  histo = (int *) malloc(numbins*sizeof(int));
  tcomm = (double *) malloc(maxiter*sizeof(double));
  tcomp = (double *) malloc(maxiter*sizeof(double));
  tstep = (double *) malloc(maxiter*sizeof(double));

  sendn = (double *) malloc(exchbytes);
  sends = (double *) malloc(exchbytes);
  sende = (double *) malloc(exchbytes);
  sendw = (double *) malloc(exchbytes);
  recvn = (double *) malloc(exchbytes);
  recvs = (double *) malloc(exchbytes);
  recve = (double *) malloc(exchbytes);
  recvw = (double *) malloc(exchbytes);

  for (i = 0; i< exchcount; i++) sendn[i] = (double) myrank;
  for (i = 0; i< exchcount; i++) sends[i] = (double) myrank;
  for (i = 0; i< exchcount; i++) sende[i] = (double) myrank;
  for (i = 0; i< exchcount; i++) sendw[i] = (double) myrank;

  for (i = 0; i< exchcount; i++) recvn[i] = 0.0;
  for (i = 0; i< exchcount; i++) recvs[i] = 0.0;
  for (i = 0; i< exchcount; i++) recve[i] = 0.0;
  for (i = 0; i< exchcount; i++) recvw[i] = 0.0;

  xrand = (double *) malloc(nrand*sizeof(double));

  seed = 13579UL;

  TMrand(seed, nrand, xrand);

  npts = (long) ( 1.0e7 * compute_interval_msec / 125.0 ); 
  ssum = 0.0;

  calib_iters = (int) ( 200.0 * 50.0 / compute_interval_msec );

  MPI_Barrier(MPI_COMM_WORLD);

  // make one call that is not timed
  compute(npts, nrand, xrand, &ssum);
  ssum += 1.0;

  if (method == EXCHANGE) exch(0);
  else                    MPI_Barrier(MPI_COMM_WORLD);

  if (myrank == 0) printf("calibrating compute time ..\n");

  MPI_Barrier(MPI_COMM_WORLD);

  // time calls to the compute routine with npts = 10^7 and adjust npts
  // compute for long enough to ramp up power
  t1 = MPI_Wtime();
  for (i = 0; i < calib_iters ; i++) {
    compute(npts, nrand, xrand, &ssum);
    ssum += 1.0;
  }
  t2 = MPI_Wtime();
  tmin = (t2 - t1)/((double) calib_iters);

  if (myrank == 0) printf("first-pass tmin = %.2lf msec\n", 1.0e3*tmin);

  // find the global minimum time and use that as the reference
  MPI_Allreduce(MPI_IN_PLACE, &tmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  // reset npts using measured compute time
  npts = (long) (((double) npts)*compute_interval_msec/(1.0e3*tmin));

  if (myrank == 0) printf("first iteration : npts = %ld\n", npts);

  // repeat this process one more time
  MPI_Barrier(MPI_COMM_WORLD);
  t1 = MPI_Wtime();
  for (i = 0; i < calib_iters ; i++) {
    compute(npts, nrand, xrand, &ssum);
    ssum += 1.0;
  }
  t2 = MPI_Wtime();
  tmin = (t2 - t1)/((double) calib_iters);

  // find the global minimum time and use that as the reference
  MPI_Allreduce(MPI_IN_PLACE, &tmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  // reset npts using measured compute time
  npts = (long) (((double) npts)*compute_interval_msec/(1.0e3*tmin));

  if (myrank == 0) printf("using npts = %ld, tmin = %.2lf msec\n", npts, 1.0e3*tmin);

  compute_interval_msec = 1.0e3*tmin;

  MPI_Barrier(MPI_COMM_WORLD);
  
  elapsed1 = MPI_Wtime();

  //==============================================
  // time a number of {compute, communicate} steps
  //==============================================
  for (iter = 0; iter < maxiter; iter++) {
    if ( (iter + 1) % 100 == 0) MPI_Barrier(MPI_COMM_WORLD);
    t1 = MPI_Wtime();
    compute(npts, nrand, xrand, &ssum);
    ssum += 1.0;
    t2 = MPI_Wtime();
    if (method == EXCHANGE) exch(iter);
    else                    MPI_Barrier(MPI_COMM_WORLD);
    t3 = MPI_Wtime();
    tcomp[iter] = t2 - t1;
    tcomm[iter] = t3 - t2;
    tstep[iter] = t3 - t1;
  }

  MPI_Allreduce(MPI_IN_PLACE, tcomm, maxiter, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  elapsed2 = MPI_Wtime();

  // max data volume per rank per step ... n,e,s,w both send and recv = 8x
  data_volume = 8.0*((double) exchbytes);

  tmax = 0.0; tmin = 1.0e30;
  for (iter = 0; iter < maxiter; iter++) {
    if (tcomm[iter] > tmax) tmax = tcomm[iter];
    if (tcomm[iter] < tmin) tmin = tcomm[iter];
  }

  if (myrank == 0) printf("max communication time in msec = %.3lf\n",  1.0e3*tmax);

  // histogram the data for max communication time per iteration
  histo_bin_width = (tmax - tmin)/((double) (numbins - 1));

  hmin = tmin - 0.5*histo_bin_width;

  for (bin = 0; bin < numbins; bin++) histo[bin] = 0;

  for (iter = 0; iter < maxiter; iter++) {
    if (histo_bin_width > 0.0)  bin = (int) ((tcomm[iter] - hmin)/histo_bin_width);
    else                        bin = 0;
    if ((bin >= 0) && (bin < numbins)) histo[bin]++;
  }

  if (myrank == 0) {
    printf("\n");
    printf("histogram of max communication times per step\n");
    printf("      msec       count\n");
    for (bin = 0; bin < numbins; bin++) {
      xhisto = tmin + histo_bin_width*((double) bin);
      printf("%10.3lf  %10d\n", xhisto*1.0e3, histo[bin]);
    }
  }

  // for per-node analysis, focus on the compute times
  compmax = (double *) malloc(maxiter*sizeof(double));

  // find the global minimum computation time over all iterations
  compmin = 1.0e30;
  for (iter = 0; iter < maxiter; iter++) {
    if (tcomp[iter] < compmin) compmin = tcomp[iter];
  }
  
  MPI_Allreduce(MPI_IN_PLACE, &compmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  // find the max compute time for each iteration on a per-node basis
  MPI_Allreduce(tcomp, compmax, maxiter, MPI_DOUBLE, MPI_MAX, local_comm);

  // find the global avg time for computation
  compavg = 0.0;
  for (iter = 0; iter < maxiter; iter++) compavg += tcomp[iter];

  MPI_Allreduce(MPI_IN_PLACE, &compavg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  samples_total = ((double) maxiter) * ((double) nranks);

  compavg = compavg / samples_total;

  // look at the max computation times per iteration for analysis by node
  compssq = 0.0;
  for (iter = 0; iter < maxiter; iter++) {
    dev = compmax[iter] - compavg;
    compssq += dev*dev;
  }

  sigma_comp = sqrt(compssq/((double) maxiter));

  num_nodes = nranks / ranks_per_node;

  color = local_rank;
  key = myrank;
  MPI_Comm_split(MPI_COMM_WORLD, color, key, &node_comm);

  if (local_rank == 0) {

     aggregate_sigma_comp = (double *) malloc(num_nodes*sizeof(double));
     aggregate_compavg    = (double *) malloc(num_nodes*sizeof(double));

     MPI_Allgather(&sigma_comp, 1, MPI_DOUBLE, aggregate_sigma_comp, 1, MPI_DOUBLE, node_comm);
     MPI_Allgather(&compavg,    1, MPI_DOUBLE, aggregate_compavg,    1, MPI_DOUBLE, node_comm);

     snames = (char *) malloc(num_nodes*sizeof(host));
     rnames = (char *) malloc(num_nodes*sizeof(host));
     sort_key = (int *) malloc(num_nodes*sizeof(int));

     gethostname(host, sizeof(host));

     for (i=0; i<sizeof(host); i++) {    
        if (host[i] == '.') {
           host[i] = '\0';
           break;
        }
     }

     for (i=0; i<num_nodes; i++)  {
       ptr = snames + i*sizeof(host);
       strncpy(ptr, host, sizeof(host));
     }

     MPI_Alltoall(snames, sizeof(host), MPI_BYTE, rnames, sizeof(host), MPI_BYTE, node_comm);

     // sort by node-local sigma for computation times
     // aggregate_sigma_comp is returned in sorted order, along with the sort key
     sortx(aggregate_sigma_comp, num_nodes, sort_key, SORT_DESCENDING_ORDER);

     if (myrank == 0) {
        printf("\n");
        printf(" percent variation = 100*sigma/mean for computation times by node:\n");
        printf("              host    percent variation\n");
        for (i=0; i< num_nodes; i++) {
           k = sort_key[i];
           ptr = rnames + k*sizeof(host);
           relative_variation = 100.0*aggregate_sigma_comp[i] / aggregate_compavg[k];
           printf("%18s    %10.3lf\n", ptr, relative_variation);
        }
        printf("\n");
     }
  }
  

  // analyze all step time data for all ranks
  tmin = 1.0e30; tmax = 0.0;
  for (iter = 0; iter < maxiter; iter++) {
     if (tstep[iter] > tmax) tmax = tstep[iter];
     if (tstep[iter] < tmin) tmin = tstep[iter];
  }

  if (method == BARRIER) {
     // in this case there is just one step time per iteration, the same for all ranks
     tavg = (elapsed2 - elapsed1) / ((double) maxiter);
     
     // compute moments for the time-step samples
     ssq = 0.0; scb = 0.0; sp4 = 0.0; 
     for (iter = 0; iter < maxiter; iter++) {
       dev = tstep[iter] - tavg;
       ssq += dev*dev;
       scb += dev*dev*dev;
       sp4 += dev*dev*dev*dev;
     }

     samples = ((double) maxiter);
  }
  else {
     // analyze all step time data for every rank
     tavg = 0.0;
     for (iter = 0; iter < maxiter; iter++) tavg += tstep[iter];

     // global avg for step times
     MPI_Allreduce(MPI_IN_PLACE, &tavg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   
     samples = ((double) maxiter) * ((double) nranks);
    
     tavg = tavg / samples;

     // compute moments for the time-step samples ... per rank
     ssq = 0.0; scb = 0.0; sp4 = 0.0; 
     for (iter = 0; iter < maxiter; iter++) {
       dev = tstep[iter] - tavg;
       ssq += dev*dev;
       scb += dev*dev*dev;
       sp4 += dev*dev*dev*dev;
     }

     // compute global statistics on the step times
     MPI_Allreduce(MPI_IN_PLACE, &ssq, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
     MPI_Allreduce(MPI_IN_PLACE, &scb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
     MPI_Allreduce(MPI_IN_PLACE, &sp4, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }

  sigma = sqrt(ssq/samples);

  relative_variation = sigma / tavg;

  skewness = scb / (samples * sigma * sigma * sigma);

  kurtosis = sp4 / (samples * sigma * sigma * sigma * sigma);

  // find the global min and max step times
  MPI_Allreduce(MPI_IN_PLACE, &tmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &tmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  histo_bin_width = (tmax - tmin)/((double) (numbins - 1));
  hmin = tmin - 0.5*histo_bin_width;

  for (bin = 0; bin < numbins; bin++) histo[bin] = 0;

  // each rank histograms its own data
  for (iter = 0; iter < maxiter; iter++) {
    if (histo_bin_width > 0.0)  bin = (int) ((tstep[iter] - hmin)/histo_bin_width);
    else                        bin = 0;
    if ((bin >= 0) && (bin < numbins)) histo[bin]++;
  }

  // for the exchange method, we should sum over all MPI ranks
  if (method == EXCHANGE) {
     MPI_Allreduce(MPI_IN_PLACE, histo, numbins, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  }

  if (myrank == 0) {
    printf("histogram of step times for all ranks\n");
    printf("      msec       count               density\n");
    for (bin = 0; bin < numbins; bin++) {
      xhisto = tmin + histo_bin_width*((double) bin);
      prob = 1.0e-3*((double) histo[bin]) / histo_bin_width;
      printf("%10.3lf  %10d  %20.4lf\n", xhisto*1.0e3, histo[bin], prob);
    }
    printf("\n");
  }

  if (myrank == 0) {
    bw = 1.0e-9*0.5*((double) ranks_per_node)*data_volume/(tavg - compmin);
    printf("min computational interval in msec = %.4lf\n", 1.0e3*compmin);
    printf("the average time per step in msec  = %.4lf\n", 1.0e3*tavg);
    printf("percent variation in step times    = %.3lf\n", 100.0*relative_variation);
    printf("skewness of step time distribution = %.3le\n", skewness);
    printf("kurtosis of step time distribution = %.3le\n", kurtosis);
    printf("elapsed time (sec) = %.3lf\n", elapsed2 - elapsed1);
    if (method == EXCHANGE) printf("effective exch bw per node = %.3le GB/sec\n", bw);
    printf("ssum =  %.6le\n\n", ssum);
  }

  MPI_Finalize();

  return 0;
}


// -----------------------------------
// routine for boundary exchange
// -----------------------------------
void exch(int iter)
{
  int ireq, tag, north, south, east, west;
  MPI_Request req[8];
  MPI_Status status[8];

  ireq = 0;

  tag = 99 + 20*iter;

  // exchange data with partner to the north
  if (myrow != (npey - 1)) {
    north = myrank + npex;
    MPI_Irecv(recvn, exchcount, MPI_DOUBLE, north, tag, MPI_COMM_WORLD, &req[ireq]);
    ireq = ireq + 1;
    MPI_Isend(sendn, exchcount, MPI_DOUBLE, north, tag, MPI_COMM_WORLD, &req[ireq]);
    ireq = ireq + 1;
  }

  // exchange data with partner to the east
  if (mycol != (npex - 1)) {
    east = myrank + 1;
    MPI_Irecv(recve, exchcount, MPI_DOUBLE, east, tag, MPI_COMM_WORLD, &req[ireq]);
    ireq = ireq + 1;
    MPI_Isend(sende, exchcount, MPI_DOUBLE, east, tag, MPI_COMM_WORLD, &req[ireq]);
    ireq = ireq + 1;
  }

  // exchange data with partner to the south
  if (myrow != 0) {
    south = myrank - npex;
    MPI_Irecv(recvs, exchcount, MPI_DOUBLE, south, tag, MPI_COMM_WORLD, &req[ireq]);
    ireq = ireq + 1;
    MPI_Isend(sends, exchcount, MPI_DOUBLE, south, tag, MPI_COMM_WORLD, &req[ireq]);
    ireq = ireq + 1;
  }

  // exchange data with partner to the west
  if (mycol != 0) {
    west = myrank - 1;
    MPI_Irecv(recvw, exchcount, MPI_DOUBLE, west, tag, MPI_COMM_WORLD, &req[ireq]);
    ireq = ireq + 1;
    MPI_Isend(sendw, exchcount, MPI_DOUBLE, west, tag, MPI_COMM_WORLD, &req[ireq]);
    ireq = ireq + 1;
  }

  MPI_Waitall(ireq, req, status);

  return;
  
}

// -----------------------------------
// routine for computation
// -----------------------------------
void compute(long n, int nrand, double * xrand, double * ssum)
{
  double tsum;
  long i, lrand, ndx;

  lrand = (long) nrand;
  
  tsum = 0.0;
  for (i = 0L; i < n; i++) {
    ndx = i % lrand;
    tsum = tsum + lutexp(xrand[ndx]);
  }

  *ssum = tsum;

}

// -----------------------------------
// help message
// -----------------------------------
void print_help(void)
{
   printf("Syntax: mpirun -np #ranks variability [-t compute_interval_msec] [-i maxiter] [-n numbins] [-x exchbytes] [-m method]\n");
}

// -----------------------------------
// thread/process binding
// -----------------------------------
#define MAX_CPUS_PER_NODE 192

//===========================================================================
// This routine will bind threads only if the env variable BIND_THREADS=yes.
// This version uses a simple cpulist : 0,1,2,..., max_cpus - 1.
// default : spread out MPI ranks and threads evenly
// options : BIND_STRIDE=number   ... stride per process not per thread
//           BIND_BASE=first logical cpu to occupy
//           BIND_CPU_LIST="cpu1,cpu2,cpu3,...,cpuN" a specific list
//===========================================================================

void bindthreads(void)
{
  char * ptr;
  int cpulist[MAX_CPUS_PER_NODE];
  int envlist[MAX_CPUS_PER_NODE];
  int uselist[MAX_CPUS_PER_NODE];
  int cpu, smt_width, max_cpus, num_cpus, inc, ndx;
  int i, j, k, processor_version, myrank, nranks, tid, skip;
  int max_ranks_per_node, base_cpu, cpus_per_rank, verbose, cpus_per_socket;
  int nthreads, rc, use_envlist, badlist;
  char * list_ptr;
  char delimiters[] = {", "};
  pthread_t thread;
  cpu_set_t cpuset;
  char * snames, * rnames, host[80];
  int bind_threads, match, color;

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);

//bind_threads = 1;
  bind_threads = 0;
 
  verbose = 0;
  ptr = getenv("BIND_VERBOSE");
  if (ptr != NULL) {
     if (strncasecmp(ptr,"yes",3) == 0) verbose = 1;
  }

  // max_cpus = max #cpus on the system
  max_cpus = sysconf(_SC_NPROCESSORS_ONLN);

#ifdef INTEL_HT
  // simple list for Intel CPUs with hyperthreading
  num_cpus = max_cpus;
  cpus_per_socket = num_cpus/2;
  k = 0;
  for (j = 0; j < 2; j++) {
     for (i = 0; i < cpus_per_socket; i++) {
        cpulist[k] = j*cpus_per_socket + i;
        cpulist[k+1] =  j*cpus_per_socket + i + num_cpus; 
        k += 2;
     }   
  }
#else
  processor_version = get_pvr();

  if      (processor_version == POWER9) smt_width = 4;
  else if (processor_version == POWER8) smt_width = 8;

  // reserve the last core from each socket if cpu isoloation is on
  ptr = getenv("LSF_CPU_ISOLATION");
  if (ptr != NULL) {
    if (0 == strncasecmp(ptr, "on", 2)) {
      num_cpus = max_cpus - 2*smt_width;
      j = 0;
      for (k=0; k<2; k++) {
        for (i=0; i<((max_cpus/2) - smt_width); i++) {
          cpulist[j] = k*(max_cpus/2) + i;
          j++;
        }
      }
    }
  }
  else {
    num_cpus = max_cpus;
    // make a simple list for all CPUs
    for (i=0; i<max_cpus; i++) cpulist[i] = i;
  }
#endif

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

  MPI_Allreduce(&ranks_per_node, &max_ranks_per_node, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  cpus_per_rank = num_cpus / max_ranks_per_node;

  nthreads = 1;

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
     }
  }
  else {
     ptr = getenv("BIND_STRIDE");
     if (ptr != NULL) cpus_per_rank = atoi(ptr);

     ptr = getenv("BIND_BASE");
     if (ptr != NULL) {
        skip = atoi(ptr);
     }
     else skip = 0;

     base_cpu = skip + cpus_per_rank * local_rank;
     inc = cpus_per_rank / nthreads;
     k = 0;
     for (i=0; i<nthreads; i++) {
        ndx = base_cpu + i*inc;
        uselist[k] = cpulist[ndx];
        k++;
     }
  }

  if (bind_threads == 0) return;

  tid = 0;
  thread = pthread_self();
  cpu = uselist[tid];
  CPU_ZERO(&cpuset);
  CPU_SET(cpu, &cpuset);
  rc = pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
  if (rc == 0 && verbose) {
     fprintf(stderr, "bindthreads: binding MPI rank %d thread %d to cpu %d\n", myrank, tid, cpu);
  }
  else if (rc != 0) {
    fprintf(stderr, "bindthreads: pthread_set_affinity_np failed for MPI rank %d, thread %d, cpu %d on host %s \n", myrank, tid, cpu, host);
  }

  return;
}

//===========================================================================
// incremental Shell sort with increment array: inc[k] = 1 + 3*2^k + 4^(k+1) 
//===========================================================================
void sortx(double * arr , int n, int * ind, int flag)
{
   int h, i, j, k, inc[20];
   int numinc, pwr2, pwr4;
   double val;

   if (n <= 1) {
      ind[0] = 0;
      return;
   }

   pwr2 = 1;
   pwr4 = 4;

   numinc = 0;
   h = 1;
   inc[numinc] = h;
   while (numinc < 20) {
      h = 1 + 3*pwr2 + pwr4;
      if (h > n) break;
      numinc++;
      inc[numinc] = h;
      pwr2 *= 2;
      pwr4 *= 4;
   }

   for (i=0; i<n; i++) ind[i] = i;

   if (flag > 0) { // sort in increasing order 
      for (; numinc >= 0; numinc--) {
         h = inc[numinc];
         for (i = h; i < n; i++) {
            val = arr[i];
            k   = ind[i];

            j = i;
   
            while ( (j >= h) && (arr[j-h] > val) ) {
               arr[j] = arr[j-h];
               ind[j] = ind[j-h];
               j = j - h;
            }

            arr[j] = val;
            ind[j] = k;
         }
      }
   }
   else { // sort in decreasing order 
      for (; numinc >= 0; numinc--) {
         h = inc[numinc];
         for (i = h; i < n; i++) {
            val = arr[i];
            k   = ind[i];

            j = i;
   
            while ( (j >= h) && (arr[j-h] < val) ) {
               arr[j] = arr[j-h];
               ind[j] = ind[j-h];
               j = j - h;
            }

            arr[j] = val;
            ind[j] = k;
         }
      }
   }
}
