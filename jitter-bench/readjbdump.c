#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int main(
    int argc,
    char** argv)
{
  // --------------------------------------------------------
  // parse and show command line arguments
  // --------------------------------------------------------

  if (argc <= 1) {
    printf("usage: %s <jbdump>\n"
           "  jbdump: path to jbdump file\n",
           argv[0]);
    exit(1);
  }

  char* bin_fn = argv[1];
  fprintf(stderr, "jbdump: %s\n", bin_fn);

  // --------------------------------------------------------
  // open file for reading
  // --------------------------------------------------------

  FILE* bin_f = fopen(bin_fn, "r");
  if (!bin_f) {
    fprintf(stderr, "error: fopen(%s): %d: %s\n",
        bin_fn, errno, strerror(errno));
    exit(1);
  }

  // --------------------------------------------------------
  // read, check, and show total ranks and phases
  // --------------------------------------------------------

  int ranks;  // total number of ranks
  int phases; // total number of phases
  double exp; // expected compute phase duration

  if ((fread(&ranks,  sizeof(int), 1, bin_f) != 1)) {
    fprintf(stderr, "error: fread(ranks): %d: %s\n",
        errno, strerror(errno));
    exit(1);
  }
  fprintf(stderr, "ranks:  %d\n", ranks);
  if (ranks < 0) {
    fprintf(stderr, "error: invalid ranks\n");
    exit(1);
  }

  if ((fread(&phases, sizeof(int), 1, bin_f) != 1)) {
    fprintf(stderr, "error: fread(phases): %d: %s\n",
        errno, strerror(errno));
    exit(1);
  }
  fprintf(stderr, "phases: %d\n", phases);
  if (phases < 0) {
    fprintf(stderr, "error: invalid phases\n");
    exit(1);
  }

  if ((fread(&exp, sizeof(double), 1, bin_f) != 1)) {
    fprintf(stderr, "error: fread(exp): %d: %s\n",
        errno, strerror(errno));
    exit(1);
  }
  fprintf(stderr, "exp:    %f\n", exp);
  if (exp < 0.0) {
    fprintf(stderr, "error: invalid exp\n");
    exit(1);
  }

  // --------------------------------------------------------
  // process jbdump
  // --------------------------------------------------------

  for (int r = 0; r < ranks; r++) {

    printf("\n");

    int tcp_s;
    int tcp_ns;

    int dac_s;
    int dac_ns;

    for (int p = 0; p < phases; p++) {

      int tap_s;
      int tap_ns;
      int tbp_s;
      int tbp_ns;

      if (p == 0) {
        fread(&tap_s,  sizeof(int), 1, bin_f);
        fread(&tap_ns, sizeof(int), 1, bin_f);
      } else {
        tap_s = tcp_s;
        tap_ns = tcp_ns;
      }

      fread(&tbp_s,  sizeof(int), 1, bin_f);
      fread(&tbp_ns, sizeof(int), 1, bin_f);

      if (p < phases - 1 || r == 0) {
        fread(&tcp_s,  sizeof(int), 1, bin_f);
        fread(&tcp_ns, sizeof(int), 1, bin_f);

        if (!r) {
          dac_s = tcp_s - tap_s;
          dac_ns = tcp_ns - tap_ns;
        }
      } else {
        tcp_s = tap_s + dac_s;
        tcp_ns = tap_ns + dac_ns;
      }

      int dt_s   = tbp_s  - tap_s;
      int dt_ns  = tbp_ns - tap_ns;
      double dt  = (double) dt_s + (double) dt_ns / 1.0E9;

      int dtc_s   = tcp_s  - tbp_s;
      int dtc_ns  = tcp_ns - tbp_ns;

      double comm  = (double) dtc_s + (double) dtc_ns / 1.0E9;
      double jit = dt - exp;

      // -----------------------------------------
      // print line
      // -----------------------------------------
      printf("%5d %10d "
             "%10d.%.09d "
             "%10d.%.09d "
             "%10d.%.09d "
             "%16.9lf %16.9lf %6.2f%% %16.9lf\n",
             r, p,
             tap_s, tap_ns,
             tbp_s, tbp_ns,
             tcp_s, tcp_ns,
             dt, jit, 100.0 * jit / exp, comm > 0 ? comm : 0);
    }
  }
}

