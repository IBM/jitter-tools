#include "tinymt64.h"

static tinymt64_t random;

// C interface
void TMrand(unsigned long seed, int npts, double * xrand)
{
   tinymt64_init(&random, seed);

   for (int i = 0; i < npts; i++)  xrand[i] = tinymt64_generate_double(&random);
}

#pragma weak tmrand_=tmrand

// Fortran interface
void tmrand(unsigned long * pseed, int * pnpts, double * xrand)
{
   unsigned long seed = *pseed;
   int npts = *pnpts;

   tinymt64_init(&random, seed);

   for (int i = 0; i < npts; i++)  xrand[i] = tinymt64_generate_double(&random);
}
