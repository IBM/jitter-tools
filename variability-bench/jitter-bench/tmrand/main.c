#include <stdio.h>

void TMrand(unsigned long, int, double *);

int main(int argc, char * argv[]) 
{
   const int npts = 30;
   unsigned long seed;
   double xrand[npts];

   seed = 1357911UL;

   TMrand(seed, npts, xrand);

   for (int i = 0; i < npts; i++)  printf("%.6lf\n", xrand[i]);

   return 0;
}
