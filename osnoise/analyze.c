#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <fcntl.h>

int main(int argc, char * argv[])
{
  FILE * fp1;
  int i, j, n, fd, ppn, nnodes;
  ssize_t rc;
  char host[80], filename[80];
  char * allhosts, * ptr;
  int cpulist[40];
  int iter, maxiter = 80000;  // 400 sec at 5 msec per sample
  double *** alldata;
  double mintime, maxtime, mean, samples;
  int rank, minhost, mincpu, maxhost, maxcpu;
  double ssq, sigma, relative_variation;
  double * compmax, * compmin;
  int bin, numbins = 50;
  double xhisto, hmin, histo_bin_width, prob;
  double avgmin, sumtime, efficiency;
  long * histo;

  ppn = 40;

  j = 0;
  for (i=0; i<80; i+=4) {
     cpulist[j] = i;
     j++;
  }
  for (i=88; i<168; i+=4) {
     cpulist[j] = i;
     j++;
  }

  printf("starting with ppn = %d\n", ppn);

  fp1 = fopen("hf", "r");
  if (fp1 == NULL) {
    printf("can't open hf ... exiting\n");
    exit(0);
  }

  n = 0;
  while (EOF != fscanf(fp1, "%s", host)) n++;

  nnodes = n;

  printf("got nnodes = %d\n", nnodes);

  rewind(fp1);

  histo = (long *) malloc(numbins*sizeof(long));

  allhosts = (char *) malloc(nnodes*sizeof(host));

  alldata = (double ***) malloc(nnodes*sizeof(double *));
  alldata[0] = (double **) malloc(nnodes*ppn*sizeof(double *));
  alldata[0][0] = (double *) malloc(nnodes*ppn*maxiter*sizeof(double));

  for (n=0; n<nnodes; n++) alldata[n] = alldata[0] + n*ppn;

  for (n=0; n<nnodes; n++) {
    for (i=0; i<ppn; i++) {
      alldata[n][i] = alldata[0][0] + i*maxiter + n*ppn*maxiter;
    }
  }
  

  // read in all data
  for (n=0; n<nnodes; n++) {
    fscanf(fp1, "%s", host);
    strcpy(allhosts + n*sizeof(host), host);
    printf("reading data for host %d = %s ...\n", n, host);
#pragma omp parallel private(i,j,filename,fd,rc)
    for (i=0; i<ppn; i++) {
      j = cpulist[i];
      sprintf(filename, "%s.%d", host, j);
      fd = open(filename, O_RDONLY);
      if (fd < 0) {
        printf("missing data file %s ... exiting\n", filename);
        exit(0);
      }
      rc = read(fd, &alldata[n][i][0], maxiter*sizeof(double));
      if (rc < 0) {
        printf("read failed for file %s ... exiting\n", filename);
        exit(0);
      }
      close(fd);
    }
  }

  // compute global min/max per rank
  mintime = 1.0e30;
  maxtime = 0.0;
  minhost = 0;
  mincpu = 0;
  maxhost = 0;
  maxcpu = 0;
  for (n=0; n<nnodes; n++) {
    for (i=0; i<ppn; i++) {
      for (iter=0; iter<maxiter; iter++) {
        mean += alldata[n][i][iter];
        if (alldata[n][i][iter] > maxtime)  {
           maxtime = alldata[n][i][iter];
           maxhost = n;   
           maxcpu = cpulist[i];
        }
        if (alldata[n][i][iter] < mintime)  {
           mintime = alldata[n][i][iter];
           minhost = n;   
           mincpu = cpulist[i];
        }
      }
    }
  }

  samples = ((double) nnodes) * ((double) ppn) * ((double) maxiter);

  mean = mean / samples;

  printf("\n");
  printf("global min time = %.3lf msec on host %s cpu %d \n", mintime, allhosts + minhost*sizeof(host), mincpu);
  printf("\n");
  printf("global max time = %.3lf msec on host %s cpu %d \n", maxtime, allhosts + maxhost*sizeof(host), maxcpu);
  printf("\n");
  printf("global avg time = %.3lf\n", mean);
  printf("\n");

  ssq = 0.0;
  for (n=0; n<nnodes; n++) {
    for (i=0; i<ppn; i++) {
      for (iter=0; iter<maxiter; iter++) {
        ssq += (alldata[n][i][iter] - mean) * (alldata[n][i][iter] - mean);
      }
    }
  }

  sigma = sqrt(ssq/samples);

  relative_variation = 100.0 * sigma / mean;

  printf("overall relative variation = %.2lf percent\n", relative_variation);

  // compute mean and relative variation by host
  samples = ((double) ppn) * ((double) maxiter); 

  // compmax is the max compute time in any given iteration
  compmax = (double *) malloc(maxiter*sizeof(double));

  // compmin is the min compute time per node
  compmin = (double *) malloc(nnodes*sizeof(double));

  for (iter=0; iter<maxiter; iter++) compmax[iter] = 0.0;
  for (n=0; n<nnodes; n++) compmin[n] = 1.0e30;
  
  printf("\n");
  printf("percent variation = 100*sigma/mean for the max computation times per step by node\n");
  printf("       host    mean(msec)   percent variation\n");
  for (n=0; n<nnodes; n++) {
    // the compmax used here is node-local
    mean = 0.0;
    for (iter=0; iter<maxiter; iter++) {
      compmax[iter] = 0.0;
      for (i=0; i<ppn; i++) {
        if (alldata[n][i][iter] > compmax[iter]) compmax[iter] = alldata[n][i][iter];
        if (alldata[n][i][iter] < compmin[n]) compmin[n] = alldata[n][i][iter];
        mean += alldata[n][i][iter];
      }
    }

    mean = mean / samples;

    ssq = 0.0;
    for (i=0; i<ppn; i++) {
      for (iter=0; iter<maxiter; iter++) {
        ssq += (alldata[n][i][iter] - compmax[iter]) * (alldata[n][i][iter] - compmax[iter]);
      }
    }

    sigma = sqrt(ssq/samples);

    relative_variation = 100.0 * sigma / mean;

    printf("%14s %10.2lf %10.2lf\n", allhosts + n*sizeof(host), mean, relative_variation);
  }

  avgmin = 0.0;
  for (n=0; n<nnodes; n++) avgmin += compmin[n];
   
  // avgmin is the min compute time per node, averaged over all nodes
  avgmin = avgmin / ((double) nnodes);

  // re-define compmax to be the max time in any rank for a given iteration
  for (iter=0; iter<maxiter; iter++) {
    compmax[iter] = 0.0;
    for (n=0; n<nnodes; n++) {
      for (i=0; i<ppn; i++) {
         if (alldata[n][i][iter] > compmax[iter]) compmax[iter] = alldata[n][i][iter];
      }
    }
  }

  // sumtime is the time expected for a parallel job
  sumtime = 0.0;
  for (iter=0; iter<maxiter; iter++) sumtime += compmax[iter];

  efficiency = ((double) maxiter) * avgmin / sumtime;

  printf("\n");
  printf("estimated overall efficiency = %.3lf\n", efficiency);

  // histogram all samples
  histo_bin_width = (maxtime - mintime)/((double) (numbins - 1));

  hmin = mintime - 0.5*histo_bin_width;

  for (bin = 0; bin < numbins; bin++) histo[bin] = 0L; 

  for (n=0; n<nnodes; n++) {
    for (i=0; i<ppn; i++) {
      for (iter=0; iter<maxiter; iter++) { 
        if (histo_bin_width > 0.0)  bin = (int) ((alldata[n][i][iter] - hmin)/histo_bin_width);
        else                        bin = 0;
        if ((bin >= 0) && (bin < numbins)) histo[bin]++;
      }
    }
  }

  printf("\n");
  printf("histogram of step times for all ranks\n");
  printf("      msec       count               density\n");
  for (bin = 0; bin < numbins; bin++) {
    xhisto = mintime + histo_bin_width*((double) bin);
    prob = 1.0e-3*((double) histo[bin]) / histo_bin_width;
    printf("%10.3lf  %10ld  %20.4lf\n", xhisto, histo[bin], prob);
  }   
  printf("\n");

  printf("summary data by rank:\n");
  printf("      host        cpu   mean(msec)  relative variation (percent)\n");
  for (n=0; n<nnodes; n++) {
    for (i=0; i<ppn; i++) {
      mean = 0.0;
      for (iter=0; iter<maxiter; iter++) mean += alldata[n][i][iter];
      mean = mean / ((double) maxiter);
      ssq = 0.0;
      for (iter=0; iter<maxiter; iter++) ssq += (alldata[n][i][iter] - mean) * (alldata[n][i][iter] - mean);
      sigma = sqrt(ssq/((double) maxiter));
      relative_variation = 100.0 * sigma / mean;
      printf("%14s %6d  %8.2lf %8.2lf\n", allhosts + n*sizeof(host), cpulist[i], mean, relative_variation);
    }
  }
  
  

  return 0;
}
