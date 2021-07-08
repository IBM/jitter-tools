#include <stdio.h>
#include <mpi.h>
#include <cuda_runtime.h>
#include "lutexp.h"
#include "lutexp.cu.h"

#define MAX_BLOCKS 256
#define THREADS_PER_BLOCK 256
#define WARP_SIZE 32
#define NUM_WARPS THREADS_PER_BLOCK/WARP_SIZE

#define CUDA_RC(rc) if( (rc) != cudaSuccess ) \
  {fprintf(stderr, "Error %s at %s line %d\n", cudaGetErrorString(cudaGetLastError()), __FILE__,__LINE__); MPI_Abort(MPI_COMM_WORLD, 1);}


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



// ----------------------------------------------------------
// routine for computation; data is returned in the out array
// which must be copied back to the host for final reduction
// ----------------------------------------------------------
__global__ void compute(long npts, int nrand, double * xrand, double * table)
{
  double tsum = 0.0;
  long lrand, ndx;
  lrand = (long) nrand;
  for (long i = blockDim.x * blockIdx.x + threadIdx.x; i < npts; i += blockDim.x * gridDim.x) {
    ndx = i % lrand;
    tsum = tsum + lutexp(xrand[ndx], table);
  }
  static volatile double keep = 0.0;
  keep += tsum; // keep compiler from discarding loop
  __syncthreads();
}


double * block_data;
double * dev_xrand, * dev_table, * dev_block_data;
int nrand;
int numBlocks;

void __config_lutexp_cuda(double *xrand, int _nrand)
{
  nrand = _nrand;

  // each MPI rank selects a device
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  int numDevices, myDevice;
  CUDA_RC(cudaGetDeviceCount(&numDevices));
  myDevice = myrank % numDevices;
  CUDA_RC(cudaSetDevice(myDevice));

  double *table = (double *) malloc(sizeof(exp_table));

  // make a copy of the exp look-up table
  for (int i = 0; i < sizeof(exp_table) / sizeof(exp_table[0]); i++) table[i] = exp_table[i];

  // allocate GPU data
  CUDA_RC(cudaMalloc((void **)&dev_table, sizeof(exp_table)));
  CUDA_RC(cudaMalloc((void **)&dev_xrand, nrand*sizeof(double)));
  CUDA_RC(cudaMalloc((void **)&dev_block_data, MAX_BLOCKS*sizeof(double)));

  block_data = (double *) malloc(MAX_BLOCKS*sizeof(double));

  // copy data to the GPU one time
  CUDA_RC(cudaMemcpy(dev_table, table, sizeof(exp_table), cudaMemcpyHostToDevice));
  CUDA_RC(cudaMemcpy(dev_xrand, xrand, nrand*sizeof(double),  cudaMemcpyHostToDevice));
}

void __compute_lutexp_cuda(uint64_t calib_iters /* == npts? */)
{
  int threadsPerBlock = THREADS_PER_BLOCK;
  int numBlocks = (calib_iters + threadsPerBlock - 1) / threadsPerBlock;
  if (numBlocks > MAX_BLOCKS) numBlocks = MAX_BLOCKS;
  compute<<<numBlocks, threadsPerBlock>>>(calib_iters, nrand, dev_xrand, dev_table);
  CUDA_RC(cudaMemcpy(block_data, dev_block_data, numBlocks*sizeof(double), cudaMemcpyDeviceToHost));
}
