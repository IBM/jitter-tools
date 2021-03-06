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
#include <sys/types.h>
#include <unistd.h>
#include <time.h>

#define BARRIER   10
#define EXCHANGE  11
#define ALLREDUCE 12

#define SORT_ASCENDING_ORDER   1
#define SORT_DESCENDING_ORDER -1

void get_local_communicator(int *, int *, MPI_Comm *);

void compute(long, int, double *, double *);

double lutexp(double);

void exch(int);

void print_help(void);

void sortx(double *, int, int *, int);

void TMrand(unsigned long seed, int npts, double * xrand);

static int myrank, myrow, mycol, npex, npey, exchcount;

static double * sendn,  * sende,  * sends,  * sendw,  * recvn,  * recve,  * recvs,  * recvw;


void bindthreads(int *, int *);


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
  int local_rank, ranks_per_node;
  MPI_Comm local_comm;

  double * tcomm, * tcomp, * tstep;
  double compute_interval_msec;
  double elapsed1, elapsed2, target_measurement_time;

  double histo_bin_width, xhisto, hmin, prob;
  double tavg, ssq, samples, relative_variation;
  double dev, sigma, scb, sp4, skewness, kurtosis;
  long * histo;

  char host[80], * ptr, * snames, * rnames;
  int * sort_key;
  int num_nodes, color, key;
  double * compmax, compmin, compavg, compssq;
  double sigma_comp, samples_total, ssq_comp;
  double * aggregate_sigma_comp, * aggregate_compavg;
  double mean, * alldev, * allavg;
  double minavg, maxavg;
  int minrank, maxrank;
  int jobid, rank, green_light, dump_data, tag = 99;
  float * compbuf, * stepbuf, * floatcomp, * floatstep;
  char outfile[80];
  time_t current_time;
  char * time_string;
  FILE * ofp;
  MPI_Status status;
  MPI_Comm node_comm;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);

  MPI_Barrier(MPI_COMM_WORLD);
  current_time = time(NULL);
  time_string = ctime(&current_time);
  if (myrank == 0) fprintf(stderr, "starting time : %s\n", time_string);

  bindthreads(&ranks_per_node, &local_rank);

  get_local_communicator(&ranks_per_node, &local_rank, &local_comm);

  num_nodes = nranks / ranks_per_node;

  // create a communicator for collecting data across nodes
  color = local_rank;
  key = myrank;
  MPI_Comm_split(MPI_COMM_WORLD, color, key, &node_comm);

  if (myrank == 0) printf("checking program args ...\n");

  // set sensible default values then check program args
  numbins = 50;
  compute_interval_msec = 5.0;
  target_measurement_time = 400.0;  // units of seconds
  exchbytes = 100;
  method = EXCHANGE;
  dump_data = 0;

  while (1) {
     ch = getopt(argc, argv, "hm:c:t:n:x:d");
     if (ch == -1) break;
     switch (ch) {
        case 'h':
           help_flag = 1;
           break;
        case 'c':  // set the compute interval in msec
           compute_interval_msec = (double) atoi(optarg);
           break;
        case 'm':  // choose the comunication method : barrier or exchange
           if (0 == strncasecmp(optarg, "barrier", 7))        method = BARRIER;
           else if (0 == strncasecmp(optarg, "allreduce", 9)) method = ALLREDUCE;
           else                                               method = EXCHANGE;
           break;
        case 't':  // set the target measurement time in sec
           target_measurement_time = atoi(optarg);
           break;
        case 'n':  // set the number of histogram bins
           numbins = atoi(optarg);
           break;
        case 'x':  // set the number of bytes for neighbor exchange
           exchbytes = atoi(optarg);
           break;
        case 'd':  // optionally dump all timing data
           dump_data = 1;
           break;
        default:
           break;
     }   
  }

  maxiter = (int) (1.0e3 * target_measurement_time / compute_interval_msec );

  if (help_flag) {
     if (myrank == 0) print_help();
     MPI_Finalize();
     return 0;
  }

  if (myrank == 0) {
     if (method == EXCHANGE) {
        printf("using neighbor exchange : compute_interval_msec = %.2lf, target measurement time = %.1lf sec, numbins = %d, exchbytes = %d\n",
                compute_interval_msec, target_measurement_time, numbins, exchbytes);
     }
     else {
        printf("using global collectives : compute_interval_msec = %.2lf, target measurement time = %.1lf sec, numbins = %d\n",
                compute_interval_msec, target_measurement_time, numbins);
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

  histo = (long *) malloc(numbins*sizeof(long));
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

  if (method == EXCHANGE)       exch(0);
  else if (method == ALLREDUCE) MPI_Allreduce(MPI_IN_PLACE, &ssum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  else                          MPI_Barrier(MPI_COMM_WORLD);

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

  // define a minimum time and use that as the reference
  // take the average-over-nodes of the node-local minimum time
  MPI_Allreduce(MPI_IN_PLACE, &tmin, 1, MPI_DOUBLE, MPI_MIN, local_comm);
  MPI_Allreduce(MPI_IN_PLACE, &tmin, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  tmin = tmin / ((double) nranks);

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

  // re-define tmin : average-over-nodes of the node-local minimum time
  MPI_Allreduce(MPI_IN_PLACE, &tmin, 1, MPI_DOUBLE, MPI_MIN, local_comm);
  MPI_Allreduce(MPI_IN_PLACE, &tmin, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  tmin = tmin / ((double) nranks);

  // reset npts using measured compute time
  npts = (long) (((double) npts)*compute_interval_msec/(1.0e3*tmin));

  if (myrank == 0) printf("using npts = %ld, tmin = %.2lf msec\n", npts, 1.0e3*tmin);

  compute_interval_msec = 1.0e3*tmin;

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Pcontrol(1);
  
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
    if (method == EXCHANGE)       exch(iter);
    else if (method == ALLREDUCE) MPI_Allreduce(MPI_IN_PLACE, &ssum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    else                          MPI_Barrier(MPI_COMM_WORLD);
    t3 = MPI_Wtime();
    tcomp[iter] = t2 - t1;
    tcomm[iter] = t3 - t2;
    tstep[iter] = t3 - t1;
  }

  MPI_Allreduce(MPI_IN_PLACE, tcomm, maxiter, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  elapsed2 = MPI_Wtime();

  MPI_Pcontrol(0);

  // compute mean and relative variation rank by rank
  mean = 0.0;
  for (iter = 0; iter < maxiter; iter++) mean += tcomp[iter];
  mean = mean / ((double) maxiter);

  ssq = 0.0;
  for (iter = 0; iter < maxiter; iter++) ssq += (tcomp[iter] - mean) * (tcomp[iter] - mean);

  sigma = sqrt(ssq/((double) maxiter));

  relative_variation = 100.0*sigma/mean;

  alldev = (double *) malloc(nranks*sizeof(double));
  allavg = (double *) malloc(nranks*sizeof(double));

  MPI_Gather(&relative_variation, 1, MPI_DOUBLE, alldev, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&mean,               1, MPI_DOUBLE, allavg, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

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

  for (bin = 0; bin < numbins; bin++) histo[bin] = 0L;

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
      printf("%10.3lf  %10ld\n", xhisto*1.0e3, histo[bin]);
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
        printf(" percent variation = 100*sigma/mean for the max computation times per step by node:\n");
        printf("              host         mean(msec)    percent variation\n");
        for (i=0; i< num_nodes; i++) {
           k = sort_key[i];
           ptr = rnames + k*sizeof(host);
           relative_variation = 100.0*aggregate_sigma_comp[i] / aggregate_compavg[k];
           printf("%18s    %10.3lf    %10.3lf\n", ptr, 1.0e3*aggregate_compavg[k], relative_variation);
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

  if (method != EXCHANGE) {
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
     MPI_Allreduce(MPI_IN_PLACE, histo, numbins, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  }

  if (myrank == 0) {
    printf("histogram of step times for all ranks\n");
    printf("      msec       count               density\n");
    for (bin = 0; bin < numbins; bin++) {
      xhisto = tmin + histo_bin_width*((double) bin);
      prob = 1.0e-3*((double) histo[bin]) / histo_bin_width;
      printf("%10.3lf  %10ld  %20.4lf\n", xhisto*1.0e3, histo[bin], prob);
    }
    printf("\n");
  }

  // optionally dump data in a single binary file, use floats to save space
  if (dump_data) {
    floatcomp = (float *) malloc(maxiter*sizeof(float));
    floatstep = (float *) malloc(maxiter*sizeof(float));
    for (iter=0; iter<maxiter; iter++) {
      floatcomp[iter] = (float) tcomp[iter];
      floatstep[iter] = (float) tstep[iter];
    }
    if (myrank == 0) {
      compbuf = (float *) malloc(maxiter*sizeof(float));
      stepbuf = (float *) malloc(maxiter*sizeof(float));
      ptr = getenv("LSB_JOBID");
      if (ptr == NULL) jobid = getpid();
      else             jobid = atoi(ptr);
      sprintf(outfile, "%d.timing_data", jobid);
      ofp = fopen(outfile, "w");
      if (ofp != NULL) {
        fwrite(floatcomp, sizeof(float), maxiter, ofp);
        fwrite(floatstep, sizeof(float), maxiter, ofp);
      }
      for (rank=1; rank<nranks; rank++) {
         MPI_Send(&green_light, 1, MPI_INT, rank, tag, MPI_COMM_WORLD);
         MPI_Recv(compbuf, maxiter, MPI_FLOAT, rank, tag, MPI_COMM_WORLD, &status);
         MPI_Recv(stepbuf, maxiter, MPI_FLOAT, rank, tag, MPI_COMM_WORLD, &status);
         if (ofp != NULL) {
           fwrite(compbuf, sizeof(float), maxiter, ofp);
           fwrite(stepbuf, sizeof(float), maxiter, ofp);
         }
      }
    }
    else {
      MPI_Recv(&green_light, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
      MPI_Send(floatcomp, maxiter, MPI_FLOAT, 0, tag, MPI_COMM_WORLD);
      MPI_Send(floatstep, maxiter, MPI_FLOAT, 0, tag, MPI_COMM_WORLD);
    }
  }

  if (myrank == 0) {
    bw = 1.0e-9*0.5*((double) ranks_per_node)*data_volume/(tavg - compmin);
    printf("min computational interval in msec = %.4lf\n", 1.0e3*compmin);
    printf("the average time per step in msec  = %.4lf\n", 1.0e3*tavg);
    printf("percent variation in step times    = %.3lf\n", 100.0*relative_variation);
    printf("skewness of step time distribution = %.3le\n", skewness);
    printf("kurtosis of step time distribution = %.3le\n", kurtosis);
    printf("elapsed time (sec)                 = %.3lf\n", elapsed2 - elapsed1);
    printf("measured time / target time        = %.3lf\n", (elapsed2 - elapsed1)/target_measurement_time); 
    printf("effective parallel efficiency      = %.3lf\n", target_measurement_time/(elapsed2 - elapsed1));
    if (method == EXCHANGE) printf("effective exch bw per node         = %.3le GB/sec\n", bw);
    printf("ssum =  %.6le\n\n", ssum);
  }

  if (myrank == 0) {
    maxavg = 0.0;
    minavg = 1.0e30;
    maxrank = 0;
    minrank = 0;
    for (i=0; i<nranks; i++) {
      if (allavg[i] > maxavg)  {
        maxavg = allavg[i];
        maxrank = i;
      }   
      if (allavg[i] < minavg)  {
        minavg = allavg[i];
        minrank = i;
      }   
    }   
    printf("\n");
    printf("max avg compute interval = %.3lf msec for rank %d\n", 1.0e3*maxavg, maxrank);
    printf("min avg compute interval = %.3lf msec for rank %d\n", 1.0e3*minavg, minrank);
    printf("\n");
    printf("mean compute interval (msec) and percent relative variation by rank:\n");
    for (i=0; i<nranks; i++) printf("%6d %12.3lf  %10.2lf\n", i, 1.0e3*allavg[i], alldev[i]);
    printf("\n");
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
   printf("Syntax: mpirun -np #ranks osnoise [-c compute_interval_msec] [-t target_measurement_time] [-n numbins] [-x exchbytes] [-m method] [-d]\n");
}


void get_local_communicator(int * ppn, int * prank, MPI_Comm * pcomm)
{
  int i, myrank, nranks;
  char * ptr, * snames, * rnames, host[80];
  int match, color;
  int local_rank, ranks_per_node;
  MPI_Comm local_comm;

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);

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
  *prank = local_rank;
  *pcomm = local_comm;

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
