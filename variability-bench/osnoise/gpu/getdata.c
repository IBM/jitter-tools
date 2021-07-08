#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

int main(int argc, char * argv[])
{
  int ifd;
  FILE * ofp;
  int rank, nranks, myrank;
  int i, npts;
  ssize_t rc, filesize;
  size_t nbytes;
  float * val;
  char infile[80];
  char outfile[80];
  off_t offset, roff;
  int bin, nbins;
  long * histo;
  float x, xmin, xmax, width;
  struct stat statbuf;

  if (argc != 2) {
    printf("syntax : getdata input_file\n");
    return 0;
  }

  sprintf(infile, argv[1]);

  rc = stat(infile, &statbuf);
  if (rc != 0) {
     printf("stat failed for the input file ... exiting\n");
     exit(0);
  }

  filesize = statbuf.st_size;

  printf("enter the number of MPI ranks : ");
  scanf("%d", &nranks);

  npts = filesize / ( 2 * nranks * sizeof(float) );

  printf("using nsamples = %d\n", npts);

  val = (float *) malloc(npts*sizeof(float));

  nbytes = npts*sizeof(float);
  ifd = open(infile, O_RDONLY);
  if (ifd < 0) {
    printf("open failed for the input file ... exiting\n");
    exit(0);
  }

//myrank = 11396;
  printf("enter the MPI rank : ");
  scanf("%d", &myrank);

  nbins = 31;

  histo = (long *) malloc((nbins + 1) * sizeof(long));

  for (bin=0; bin<=nbins; bin++) histo[bin] = 0L;

//for (rank = 0; rank < 24576; rank++) {

//  if (rank == myrank) continue;
    rank = myrank;

    offset = ((long) 8) * ((long) rank) * ((long) npts);

    roff = lseek(ifd, offset, SEEK_SET);
    if (roff < 0) {
      printf("lseek failed for the input file ... exiting\n");
      exit(0);
    }

    rc = read(ifd, val, nbytes);
    if (rc != nbytes) {
      printf("read failed for the input file ... exiting\n");
      exit(0);
    }

    sprintf(outfile, "%d.comp", myrank);
    ofp = fopen(outfile, "w");

    xmin = 1.0e30;
    xmax = 0.0;
    for (i=0; i<npts; i++) {
      if (val[i] < xmin) xmin = val[i];
      if (val[i] > xmax) xmax = val[i];
    }
    printf("got xmin = %.3le, xmax = %.3le\n", xmin, xmax);

    width = (xmax - xmin) / ((double) nbins);

    for (i=0; i<npts; i++) {
      fprintf(ofp, "%6d %.6e\n", i, val[i]);
      bin = (val[i] - xmin) / width;
      if (bin >= 0  && bin <= nbins) histo[bin]++;
    }   


    fclose(ofp);

//}

  rc = close(ifd);
  if (rc < 0) {
    printf("close failed for the input file\n");
  }

/*
  x = xmin;
  printf("%.3e  %d\n", x, 0L); 
  for (bin=0; bin<nbins; bin++) {
    x = xmin + width * ((float) bin);
    printf("%.3e  %ld\n", x, histo[bin]);
    x += width;
    printf("%.3e  %ld\n", x, histo[bin]);
  }   
  x = xmax;
  printf("%.3e  %ld\n", x, 0L); 
*/

  for (bin=0; bin<=nbins; bin++) {
    x = xmin + width * ((float) bin) + 0.5f*width;
    printf("%.3e  %ld\n", x, histo[bin]);
  }   
  return 0;
}
