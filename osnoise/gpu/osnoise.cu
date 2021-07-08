#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <pthread.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <errno.h>
#include <string.h>
#include <cuda_runtime.h>
#include <nvml.h>
#include "lutexp.h"

#define MAX_BLOCKS 256
#define THREADS_PER_BLOCK 256
#define WARP_SIZE 32
#define NUM_WARPS THREADS_PER_BLOCK/WARP_SIZE

#define BARRIER  10
#define EXCHANGE 11

#define SORT_ASCENDING_ORDER   1
#define SORT_DESCENDING_ORDER -1

#define CUDA_RC(rc) if( (rc) != cudaSuccess ) \
  {fprintf(stderr, "Error %s at %s line %d\n", cudaGetErrorString(cudaGetLastError()), __FILE__,__LINE__); MPI_Abort(MPI_COMM_WORLD, 1);}

void get_local_communicator(int *, int *, MPI_Comm *);

__global__ void compute(long, int, double *, double *, double *);

__device__ double lutexp(double, double *);

void exch(int);

void print_help(void);

void sortx(double *, int, int *, int);

extern "C" void TMrand(unsigned long seed, int nrand, double * xrand);

static int myrank, myrow, mycol, npex, npey, exchcount;

static double * sendn,  * sende,  * sends,  * sendw,  * recvn,  * recve,  * recvs,  * recvw;


/*=============================================================*/
/* Algorithm for exp(x):                                       */
/*                                                             */
/*       x*recln2 = p + r                                      */
/*       exp(x) = (2**p) * exp(r*ln2)                          */
/*                                                             */
/*       where p = integer + n-bit fraction                    */
/*       and r is a small remainder.                           */
/*                                                             */
/*       The integer part of p goes in the exponent field,     */
/*       the n-bit fraction is handled with a lookup table.    */
/*       The remainder term gets a Taylor expansion.           */
/*                                                             */
/* Use an 8-bit fraction, and a 5th order Taylor series.       */
/*                                                             */
/* p = (2**44 + 2**43) + (x*recln2) - (2**44 + 2**43)          */
/* gives 8-bit accuracy in the n-bit fraction                  */
/*=============================================================*/

__device__ double lutexp(double x, double * table)
{
   int hind;
   double result, hfactor;
   double recln2 = 1.44269504088896340736;
   double twop44_plus_twop43 = 2.6388279066624000e+13;

   struct uPair { unsigned lo ; unsigned hi ; };

   union { double d; struct uPair up; } X;
   union { double d; struct uPair up; } Result;

   double t, poly;
   double f0, f1, f2, f3, f4;
   double c2 = 1.0/2.0;
   double c3 = 1.0/6.0;
   double c4 = 1.0/24.0;
   double c5 = 1.0/120.0;

   /*--------------------------------------------------------*/
   /* multiply the input value by the reciprocal of ln(2)    */
   /* and shift by 2**44 + 2**43; save the exponent          */
   /*--------------------------------------------------------*/
   X.d = x*recln2 + twop44_plus_twop43;
   Result.up.hi = ( ( (X.up.lo >> 8) + 1023 )  << 20 ) & 0x7ff00000;
   Result.up.lo = 0;

   /*--------------------------------------------------*/
   /* compute the small remainder for the polynomial   */
   /* use the last 8 bits of the shifted X as an index */
   /*--------------------------------------------------*/
   t = x - (X.d - twop44_plus_twop43)*M_LN2;
   hind = X.up.lo & 0x000000ff;
   hfactor = Result.d * table[hind];

// /*---------------------------------------*/
// /* use a polynomial expansion for exp(t) */
// /*---------------------------------------*/
// poly = 1.0 + t*(c1 + t*(c2 + t*(c3 + t*(c4 + t*c5))));

   /*----------------------------------------------------*/
   /* for a single call, better to factor the polynomial */
   /*----------------------------------------------------*/
   f0 = 1.0 + t;
   f2 = t*t;
   f4 = f2*f2;
   f1 = c2 + t*c3;
   f3 = c4 + t*c5;

   poly = (f0 + f1*f2) + f3*f4;

   /*---------------------------------*/
   /* construct the result and return */
   /*---------------------------------*/
   result = hfactor* poly;

   /*-----------------------------------*/
   /* check input value for valid range */
   /*-----------------------------------*/
   if      (x < -709.0) return 0.0;
   else if (x >  709.0) return HUGE_VAL;
   else return result;
}

// TRACING =================================================

#include <time.h>

double timediff(
    struct timespec* t1,
    struct timespec* t0)
{
  struct timespec dt;

  if (t1->tv_nsec >= t0->tv_nsec) {
    dt.tv_sec  = t1->tv_sec  - t0->tv_sec;
    dt.tv_nsec = t1->tv_nsec - t0->tv_nsec;
  }
  else {
    dt.tv_sec  = t1->tv_sec  - t0->tv_sec  - 1;
    dt.tv_nsec = t1->tv_nsec - t0->tv_nsec + 1E+9;
  }

  return (double) dt.tv_sec + (double) dt.tv_nsec / 1.0E+9;
}

// =========================================================

int main(int argc, char * argv[])
{
  int i, k, exchbytes, nranks, ch, help_flag = 0;
  int iter, maxiter, calib_iters, bin, numbins, method;
  long npts;
  double  ssum, data_volume, bw;
  double  t1, t2, /*t3,*/ tmin, tmax;
  double * xrand, * table, * block_data;
  double * dev_xrand, * dev_table, * dev_block_data;
  double mean, * alldev, * allavg;
  double minavg, maxavg;
  int minrank, maxrank;
  unsigned long seed;
  int nrand = 1000, ntable = 256;
  int local_rank, ranks_per_node;
  MPI_Comm local_comm;

  // TRACING
  int tracing = 0;
  double threshold = INFINITY;

  double * tcomm, * tcomp, * tstep;
  double compute_interval_msec;
  double elapsed1, elapsed2, target_measurement_time;

  double histo_bin_width, xhisto, hmin, prob;
  double tavg, ssq, samples, relative_variation;
  double dev, sigma, scb, sp4, skewness, kurtosis;
  int * histo;

  int jobid, rank, green_light, dump_data = 0, tag = 99;
  float * compbuf, * stepbuf, * floatcomp, * floatstep;
  char outfile[80];
  MPI_Status status;

  char host[80], * ptr, * snames, * rnames;
  int * sort_key;
  int num_nodes, color, key;
  double * compmax, compmin, compavg, compssq;
  double sigma_comp, samples_total;
  double * aggregate_sigma_comp, * aggregate_compavg;
  time_t current_time;
  char * time_string;
  FILE * ofp;
  MPI_Comm node_comm;

  int adjust_freq = 0;
  nvmlDevice_t nvmldevice, * device;
  unsigned int device_count, graphicsClockMHz;
  int numDevices, myDevice, gpu;

  cudaEvent_t cudaBeg, cudaEnd;
  //float kernel_msec;
  int kernel_timing;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);
  get_local_communicator(&ranks_per_node, &local_rank, &local_comm);

  MPI_Barrier(MPI_COMM_WORLD);
  current_time = time(NULL);
  time_string = ctime(&current_time);
  if (myrank == 0) printf("starting time : %s\n", time_string);

  // each MPI rank selects a device
  CUDA_RC(cudaGetDeviceCount(&numDevices));
  myDevice = myrank % numDevices;
  CUDA_RC(cudaSetDevice(myDevice));

  if (myrank == 0) printf("checking program args ...\n");

  // set sensible default values then check program args
  numbins = 50;
  compute_interval_msec = 50.0;
  target_measurement_time = 50.0;  // units of seconds
  exchbytes = 100;
  method = EXCHANGE;
  kernel_timing = 0;

  while (1) {
     ch = getopt(argc, argv, "hdkf:m:c:t:n:x:j:");
     if (ch == -1) break;
     switch (ch) {
        case 'h':
           help_flag = 1;
           break;
        case 'c':  // set the compute interval in msec
           compute_interval_msec = (double) atoi(optarg);
           break;
        case 'm':  // choose the comunication method : barrier or exchange
           if (0 == strncasecmp(optarg, "barrier", 7)) method = BARRIER;
           else                                        method = EXCHANGE;
           break;
        case 't':  // set the target measurement time in sec
           target_measurement_time = atoi(optarg);
           break;
        case 'f':  // set the GPU frequency
           adjust_freq = 1;
           graphicsClockMHz = atoi(optarg);
           break;
        case 'n':  // set the number of histogram bins
           numbins = atoi(optarg);
           break;
        case 'x':  // set the number of bytes for neighbor exchange
           exchbytes = atoi(optarg);
           break;
        /*
        case 'k':  // optionally use kernel time for tcomp
           kernel_timing = 1;
           break;
        */
        case 'j':  // set jitter tracing threshold value in msec
           tracing = 1;
           threshold = atof(optarg) / 1000.0;
           break;
        case 'd':  // optionally dump all timing data
           dump_data = 1;
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

  maxiter = (int) (1.0e3 * target_measurement_time / compute_interval_msec );

  if (kernel_timing) {
     CUDA_RC(cudaEventCreate(&cudaBeg));
     CUDA_RC(cudaEventCreate(&cudaEnd));
     if (myrank == 0) printf("compute times are cuda kernel times, not counting memcpy\n");
  }
  else {
     if (myrank == 0) printf("compute times include cuda kernel plus memcpy times\n");
  }

  if (NVML_SUCCESS != nvmlInit()) {
     fprintf(stderr, "failed to initialize NVML ... exiting\n");
     MPI_Abort(MPI_COMM_WORLD, 1);
  }

  if (NVML_SUCCESS != nvmlDeviceGetCount(&device_count)) {
     fprintf(stderr, "nvmlDeviceGetCount failed ... exiting\n");
     MPI_Abort(MPI_COMM_WORLD, 1);
  }

  device = (nvmlDevice_t *) malloc(device_count*sizeof(nvmlDevice_t));

  for (i = 0; i < device_count; i++) {
     if (NVML_SUCCESS != nvmlDeviceGetHandleByIndex(i, &device[i])) {
        fprintf(stderr, "nvmlDeviceGetHandleByIndex failed ... exiting\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
     }
  }

  // for use in nvml queries
  nvmldevice = device[myDevice];

  if (adjust_freq) {
     if (NVML_SUCCESS != nvmlDeviceSetApplicationsClocks ( nvmldevice, 877, graphicsClockMHz )) {
        fprintf(stderr, "rank %d failed to set the GPU freq to %d for device %d ... continuing\n", myrank, graphicsClockMHz, myDevice);
     }
  }

  if (myrank == 0) {
     if (method == EXCHANGE) {
        printf("using neighbor exchange : compute_interval_msec = %.2lf, target measurement time = %.1lf sec, numbins = %d, exchbytes = %d\n",
                compute_interval_msec, target_measurement_time, numbins, exchbytes);
     }
     else {
        printf("using global barrier : compute_interval_msec = %.2lf, target measurement time = %.1lf sec, numbins = %d\n",
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

  table = (double *) malloc(ntable*sizeof(double));

  // make a copy of the exp look-up table
  for (i = 0; i < ntable; i++) table[i] = exp_table[i];

  // allocate GPU data
  CUDA_RC(cudaMalloc((void **)&dev_table, ntable*sizeof(double)));
  CUDA_RC(cudaMalloc((void **)&dev_xrand, nrand*sizeof(double)));
  CUDA_RC(cudaMalloc((void **)&dev_block_data, MAX_BLOCKS*sizeof(double)));

  block_data = (double *) malloc(MAX_BLOCKS*sizeof(double));

  // copy data to the GPU one time
  CUDA_RC(cudaMemcpy(dev_table, table, ntable*sizeof(double), cudaMemcpyHostToDevice));
  CUDA_RC(cudaMemcpy(dev_xrand, xrand, nrand*sizeof(double),  cudaMemcpyHostToDevice));

  // this initial guess is pretty good for Volta ... adjust if necessary
  npts = (long) ( 1.6e10 * compute_interval_msec / 125.0 );
  ssum = 0.0;

  if (myrank == 0) printf("initial npts = %ld\n", npts);

  calib_iters = (int) ( 200.0 * 50.0 / compute_interval_msec );

  MPI_Barrier(MPI_COMM_WORLD);

  // make one call that is not timed
  int threadsPerBlock = THREADS_PER_BLOCK;
  int numBlocks = (npts + threadsPerBlock - 1) / threadsPerBlock;
  if (numBlocks > MAX_BLOCKS) numBlocks = MAX_BLOCKS;
  compute<<<numBlocks, threadsPerBlock>>>(npts, nrand, dev_xrand, dev_table, dev_block_data);
  CUDA_RC(cudaMemcpy(block_data, dev_block_data, numBlocks*sizeof(double), cudaMemcpyDeviceToHost));
  for (int j = 0; j < numBlocks; j++)  ssum += block_data[j];

  ssum += 1.0;

  if (method == EXCHANGE) exch(0);
  else                    MPI_Barrier(MPI_COMM_WORLD);

  if (myrank == 0) printf("calibrating compute time ..\n");

  MPI_Barrier(MPI_COMM_WORLD);

  // time calls to the compute routine and adjust npts
  // compute for long enough to ramp up power
  t1 = MPI_Wtime();
  for (i = 0; i < calib_iters ; i++) {
    if (kernel_timing) CUDA_RC(cudaEventRecord(cudaBeg));
    compute<<<numBlocks, threadsPerBlock>>>(npts, nrand, dev_xrand, dev_table, dev_block_data);
    if (kernel_timing) CUDA_RC(cudaEventRecord(cudaEnd));
    CUDA_RC(cudaMemcpy(block_data, dev_block_data, numBlocks*sizeof(double), cudaMemcpyDeviceToHost));
    for (int j = 0; j < numBlocks; j++)  ssum += block_data[j];
    ssum += 1.0;
  }
  t2 = MPI_Wtime();
  tmin = (t2 - t1)/((double) calib_iters);

  if (myrank == 0) printf("initial tmin = %.2lf msec\n", 1.0e3*tmin);

  // find the global minimum time and use that as the reference
  MPI_Allreduce(MPI_IN_PLACE, &tmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  // reset npts using measured compute time
  npts = (long) (((double) npts)*compute_interval_msec/(1.0e3*tmin));

  if (myrank == 0) printf("first iteration : npts = %ld\n", npts);

  numBlocks = (npts + threadsPerBlock - 1) / threadsPerBlock;
  if (numBlocks > MAX_BLOCKS) numBlocks = MAX_BLOCKS;

  // repeat this process one more time
  MPI_Barrier(MPI_COMM_WORLD);
  t1 = MPI_Wtime();
  for (i = 0; i < calib_iters ; i++) {
    if (kernel_timing) CUDA_RC(cudaEventRecord(cudaBeg));
    compute<<<numBlocks, threadsPerBlock>>>(npts, nrand, dev_xrand, dev_table, dev_block_data);
    if (kernel_timing) CUDA_RC(cudaEventRecord(cudaEnd));
    CUDA_RC(cudaMemcpy(block_data, dev_block_data, numBlocks*sizeof(double), cudaMemcpyDeviceToHost));
    for (int j = 0; j < numBlocks; j++)  ssum += block_data[j];
    ssum += 1.0;
  }
  t2 = MPI_Wtime();
  tmin = (t2 - t1)/((double) calib_iters);

  // find the global minimum time and use that as the reference
  MPI_Allreduce(MPI_IN_PLACE, &tmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  // reset npts using measured compute time
  npts = (long) (((double) npts)*compute_interval_msec/(1.0e3*tmin));

  if (myrank == 0) printf("using npts = %ld, tmin = %.2lf msec\n", npts, 1.0e3*tmin);

  numBlocks = (npts + threadsPerBlock - 1) / threadsPerBlock;
  if (numBlocks > MAX_BLOCKS) numBlocks = MAX_BLOCKS;

  compute_interval_msec = 1.0e3*tmin;

  // TRACING
  struct timespec* T1 = (struct timespec*) malloc(maxiter * sizeof(struct timespec));
  struct timespec* T2 = (struct timespec*) malloc(maxiter * sizeof(struct timespec));
  struct timespec* T3 = (struct timespec*) malloc(maxiter * sizeof(struct timespec));

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Pcontrol(1);

  elapsed1 = MPI_Wtime();

  //==============================================
  // time a number of {compute, communicate} steps
  //==============================================

  for (iter = 0; iter < maxiter; iter++) {
    if ( (iter + 1) % 100 == 0) MPI_Barrier(MPI_COMM_WORLD);

clock_gettime(CLOCK_MONOTONIC, &T1[iter]);
    compute<<<numBlocks, threadsPerBlock>>>(npts, nrand, dev_xrand, dev_table, dev_block_data);
    CUDA_RC(cudaMemcpy(block_data, dev_block_data, numBlocks*sizeof(double), cudaMemcpyDeviceToHost));
clock_gettime(CLOCK_MONOTONIC, &T2[iter]);
    if (method == EXCHANGE)     exch(iter);
    else                        MPI_Barrier(MPI_COMM_WORLD);
clock_gettime(CLOCK_MONOTONIC, &T3[iter]);
  }

  /*
  for (iter = 0; iter < maxiter; iter++) {
    if ( (iter + 1) % 100 == 0) MPI_Barrier(MPI_COMM_WORLD);
    t1 = MPI_Wtime();
    if (kernel_timing) CUDA_RC(cudaEventRecord(cudaBeg));
    compute<<<numBlocks, threadsPerBlock>>>(npts, nrand, dev_xrand, dev_table, dev_block_data);
    if (kernel_timing) CUDA_RC(cudaEventRecord(cudaEnd));
    CUDA_RC(cudaMemcpy(block_data, dev_block_data, numBlocks*sizeof(double), cudaMemcpyDeviceToHost));
    if (kernel_timing) CUDA_RC(cudaEventSynchronize(cudaEnd));
    for (int j = 0; j < numBlocks; j++)  ssum += block_data[j];
    ssum += 1.0;
    t2 = MPI_Wtime();
    if (method == EXCHANGE) exch(iter);
    else                    MPI_Barrier(MPI_COMM_WORLD);
    t3 = MPI_Wtime();
    if (kernel_timing) {
       CUDA_RC(cudaEventElapsedTime(&kernel_msec, cudaBeg, cudaEnd));
       tcomp[iter] = 1.0e-3 * ((double) kernel_msec);
    }
    else tcomp[iter] = t2 - t1;
    tcomm[iter] = t3 - t2;
    tstep[iter] = t3 - t1;
  }
  */

  MPI_Allreduce(MPI_IN_PLACE, tcomm, maxiter, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  elapsed2 = MPI_Wtime();

  MPI_Pcontrol(0);

  for (iter = 0; iter < maxiter; iter++) {
    tcomp[iter] = timediff(&T2[iter], &T1[iter]);
    tcomm[iter] = timediff(&T3[iter], &T2[iter]);
    tstep[iter] = timediff(&T3[iter], &T1[iter]);
  }

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

  MPI_Gather(&relative_variation,     1, MPI_DOUBLE, alldev,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&mean,                   1, MPI_DOUBLE, allavg,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

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
    printf("  rank      host   gpu  <compute(msec)> %%variation\n");
    for (i=0; i<nranks; i++) {
       k = i / numDevices;
       ptr = rnames + k*sizeof(host);
       gpu = i % numDevices;
       printf("%6d %10s %4d %12.3lf  %10.2lf\n", i, ptr, gpu, 1.0e3*allavg[i], alldev[i]);
    }
    printf("\n");
  }

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

  //
  // TRACING
  // 

  if (tracing) {
    char hostname[64] = {0};
    gethostname(hostname, 64);
    int cpu;
    cpu_set_t set;
    sched_getaffinity(0, sizeof(set), &set);
    int j;
    for (j = 0; j < CPU_SETSIZE; ++j) {
      if (CPU_ISSET(j, &set)) {
        cpu = j;    
        break;
      }
    }

    if (myrank == 0) {
      struct timespec* T1_buf = (struct timespec*) malloc(maxiter * sizeof(struct timespec));
      struct timespec* T2_buf = (struct timespec*) malloc(maxiter * sizeof(struct timespec));
      struct timespec* T3_buf = (struct timespec*) malloc(maxiter * sizeof(struct timespec));
      memcpy(T1_buf, T1, maxiter * sizeof(struct timespec));
      memcpy(T2_buf, T2, maxiter * sizeof(struct timespec));
      memcpy(T3_buf, T3, maxiter * sizeof(struct timespec));

      char rank_host[64] = {0};
      memcpy(rank_host, hostname, 64);
      int rank_cpu;
      rank_cpu = cpu;

      char tracing_fn[256] = {0};
      sprintf(tracing_fn, "job%s.jtf", getenv("LSB_JOBID"));
      FILE* tracing_fp = fopen(tracing_fn, "w");

      for (rank = 0; rank < nranks; rank++) {
        if (rank) {
          MPI_Recv(rank_host, 64, MPI_BYTE, rank, tag, MPI_COMM_WORLD, &status);
          MPI_Recv(&rank_cpu, 1,  MPI_INT,  rank, tag, MPI_COMM_WORLD, &status);
          MPI_Recv(T1_buf, maxiter * sizeof(struct timespec), MPI_BYTE, rank, tag, MPI_COMM_WORLD, &status);
          MPI_Recv(T2_buf, maxiter * sizeof(struct timespec), MPI_BYTE, rank, tag, MPI_COMM_WORLD, &status);
          MPI_Recv(T3_buf, maxiter * sizeof(struct timespec), MPI_BYTE, rank, tag, MPI_COMM_WORLD, &status);
        }

        printf("[jtf] rank=%d host=%s cpu=%d maxiter=%d\n", rank, rank_host, rank_cpu, maxiter);

        for (iter = 0; iter < maxiter; iter++) {
          double dt = timediff(&T2_buf[iter], &T1_buf[iter]);
          if (dt >= threshold) {
            fprintf(tracing_fp, 
                "%-10d %s %-3d %10ld.%06ld %10ld.%06ld %16.6lf\n",
                rank,
                rank_host,
                rank_cpu,
                T1_buf[iter].tv_sec, T1_buf[iter].tv_nsec / 1000,
                T2_buf[iter].tv_sec, T2_buf[iter].tv_nsec / 1000,
                dt);
          } // threshold
        } // iter loop
      } // rank loop
    } // rank 0
    else {
      MPI_Send(hostname, 64, MPI_BYTE, 0, tag, MPI_COMM_WORLD);
      MPI_Send(&cpu,     1,  MPI_INT,  0, tag, MPI_COMM_WORLD);
      MPI_Send(T1, maxiter * sizeof(struct timespec), MPI_BYTE, 0, tag, MPI_COMM_WORLD);
      MPI_Send(T2, maxiter * sizeof(struct timespec), MPI_BYTE, 0, tag, MPI_COMM_WORLD);
      MPI_Send(T3, maxiter * sizeof(struct timespec), MPI_BYTE, 0, tag, MPI_COMM_WORLD);
    } // rank x
  } // tracing

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

// ----------------------------------------------------------
// routine for computation; data is returned in the out array
// which must be copied back to the host for final reduction
// ----------------------------------------------------------
__global__ void compute(long npts, int nrand, double * xrand, double * table, double * out)
{
  double tsum, y;
  long lrand, ndx;
  __shared__ double acc[NUM_WARPS];
  int lane   = threadIdx.x % WARP_SIZE;
  int warp   = threadIdx.x / WARP_SIZE;

  lrand = (long) nrand;

  tsum = 0.0;
  for (long i = blockDim.x * blockIdx.x + threadIdx.x;  i < npts; i += blockDim.x * gridDim.x) {
    ndx = i % lrand;
    y = lutexp(xrand[ndx], table);
    tsum += y;
  }

  // reduce over threads in a warp
  for (int shift = WARP_SIZE/2; shift > 0; shift /= 2) tsum += __shfl_down_sync(0xffffffff, tsum, shift);

  // save values for this thread block in shared memory
  if (lane == 0) acc[warp] = tsum;
  __syncthreads();

  // reduce once more to get a value per thread block
  if (warp == 0) {
     tsum = (lane < NUM_WARPS) ? acc[lane] : 0.0;
     for (int shift = NUM_WARPS/2; shift > 0; shift /= 2) tsum += __shfl_down_sync(0xffffffff, tsum, shift);
  }

  // save the per-block values for final reduction on the host
  if (threadIdx.x == 0) out[blockIdx.x] = tsum;

}

// -----------------------------------
// help message
// -----------------------------------
void print_help(void)
{
   printf("Syntax: mpirun -np #ranks osnoise [-d] [-c compute_interval_msec] [-t target_measurement_time] [-n numbins] [-x exchbytes] [-m method]\n");
   printf("  -d to dump data\n"
          "  -t in units of sec\n"
          "  -x in units of bytes\n");
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
