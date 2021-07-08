 /* Copyright IBM corp.
 * author: Alessandro Morari, Bryan Rosenburg, Nelson Mimura Gonzalez, Patrick Siegl
 *
 */


#define _GNU_SOURCE

#include <assert.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <sched.h>
#include <signal.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/time.h>
#include <sys/file.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>

#include <mpi.h>

#include "../jitter-trace/libjt.c"

static int              rank;
static int              np;
static int              COMPUTE_BEST_TRIAL              = 2000;
static int              per_rank_stats                  = 0;
static char             perror_str[2048];

double  lutexp(double);
void    TMrand(unsigned long seed, int npts, double * xrand);
int musl_execve(const char *payload_elf,
                const char **payload_argv);

#define nrand 1000
double          xrand[nrand];
unsigned long   seed  = 13579UL;

#define USE_MPI_WTIME 0
#define HIST_NUM_BINS 20
#define MAX_HOSTNAME_CHAR 255

typedef enum comp_kern_t {
  COMPUTE_LUTEXP      = 0,
  COMPUTE_LUTEXP_CUDA = 1,
  COMPUTE_SQRT        = 2,
  COMPUTE_ADD         = 3
} comp_kern_t;

#define COMM_NONE      0
#define COMM_BARRIER   1
#define COMM_ALLREDUCE 2


void errExit(const char *msg, ...)
{
    va_list args;
    va_start(args, msg);
    vfprintf(stderr, msg, args);
    va_end(args);

    MPI_Finalize();
    exit(EXIT_FAILURE);
}

#define perrorMsg(msg) do { snprintf(perror_str, 1024, "[rank %d] %s", rank, msg); perror(perror_str);\
                        } while (0)

#define SORT_ASCENDING_ORDER   1
#define SORT_DESCENDING_ORDER -1

typedef struct node_info_t {
    char hostname[MAX_HOSTNAME_CHAR];
    int cpu;
} node_info_t;

static node_info_t  node_info;
static node_info_t* node_info_array =  NULL;

struct  calibration{
    uint64_t iters;
    uint64_t num_phases;
    double compute_usec;
};


void usage(const char *name)
{
    if ( ! rank)
        printf("usage: %s -p <number of phases>\n"
               "  Pinning options (use both):\n"
               "    -l <#ranks> : pin #ranks per node\n"
               "    -s <#cpus>  : pin with a stride of #cpus\n"
               "  Run options:\n"
               "    -m lutexp|lutexp_cuda|sqrt|add : compute phase math kernel\n"
               "    -c <usec>       : compute phase microseconds\n"
               "    -i <iterations> : compute phase iterations\n"
               "        specify -c or -i\n"
               "    -r <comm>: communication between compute phases\n"
               "            'n' : none\n"
               "            'b' : barrier (default)\n"
               "            'a' : all-reduce\n"
               "    -C mono|real : clock to use for timestamps\n"
               "            mono: CLOCK_MONOTONIC (default)\n"
               "            real: CLOCK_REALTIME (time-of-day))\n"
               "    -a <binary> <args...> : provide a payload binary followed\n"
               "                            by its own arguments.\n"
               "                            Payload is executed as follows:\n"
               "                            ┌ %s calibration\n"
               "                            ├ execution of <binary> <args...>\n"
               "                            └ %s benchmark\n"
               "  Output options:\n"
               "    -n: print per-rank histogram\n"
               "    -d <usec>: print raw data for phases longer than <usec>\n"
               "    -e: print extra statistics\n"
               "    -b: dump raw data in binary format\n"
               "  Tracing options:\n"
               "    -t:     enable filtered tracing\n",
               name, name, name );
}

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
        if (h > n)
            break;
        numinc++;
        inc[numinc] = h;
        pwr2 *= 2;
        pwr4 *= 4;
    }

    for (i=0; i<n; i++)
        ind[i] = i;

    // sort in increasing order
    if (flag > 0) {
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
    // sort in decreasing order
    else {
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

int timestamp_clock = CLOCK_MONOTONIC;

static inline double time_sec(void)
{
#if (USE_MPI_WTIME)
    return MPI_Wtime();
#else
    struct timespec ts;
    if (clock_gettime(timestamp_clock, &ts)) {
        errExit("clock_gettime error");
    }
    return ((double)ts.tv_sec + (double)ts.tv_nsec / 1000000000.0);
#endif
}

int get_cpu(void)
{
    cpu_set_t set;
    if( sched_getaffinity(0, sizeof(set), &set) == -1 )
        errExit("sched_getaffinity error");
    int j;
    for (j = 0; j < CPU_SETSIZE; ++j){
        if (CPU_ISSET(j, &set))
            return j;
    }
    errExit("Error:cannot find cpu");
    return -1;
}

void set_cpu(int cpu)
{
    cpu_set_t set;
    CPU_ZERO(&set);
    CPU_SET(cpu, &set);
    printf("rank %d setting cpu %d\n", rank, cpu);
    if (sched_setaffinity(0, sizeof(set), &set) == -1)
        errExit("sched_setaffinity error");
}

typedef struct  {
  uint64_t count[HIST_NUM_BINS];
  uint64_t time[HIST_NUM_BINS];
  uint64_t samples;
  double bin_size;
  uint64_t aggregate;
  double min_value;
  double max_value;
  double local_min;
  double local_max;
  double mean_value;
  double sum_squares;
  double perc_var;
} histogram_t;

void print_hist(histogram_t * hist, const char * title){
  printf("%s\nsamples: %lu  aggregate: %lu usec\n",
      title, hist->samples, hist->aggregate);
  printf("min: %.2lf usec  avg: %.2lf usec  max: %.2lf usec\n",
      hist->local_min, hist->mean_value, hist->local_max);
  printf("variation (100*sigma/mean): %.3lf%%\n",hist->perc_var);
  printf("%19s \t%10s\t%12s\t%9s\n",
         "bin range", "count", "time (usec)", "time %");
  printf("------------------- \t----------\t------------\t---------\n");
  uint64_t i;
  for(i = 0; i < HIST_NUM_BINS; ++i) {
      printf("%9.2lf %9.2lf:\t%10lu\t%12lu\t%8.4lf%%\n",
             exp(log(hist->min_value) + (i * hist->bin_size)),
             exp(log(hist->min_value) + ((i+1) * hist->bin_size)),
             hist->count[i],
             hist->time[i],
             100.0 * (((double) hist->time[i]) / ((double) hist->aggregate)));
  }
}


void min_max(double *series, uint64_t series_size, double *_min, double *_max)
{
    double min = UINT64_MAX;
    double max = 0;
    uint64_t i;
    for (i = 0; i < series_size; ++i) {
        min = (min < series[i])?min:series[i];
        max = (max > series[i])?max:series[i];
    }
    *_min = min;
    *_max = max;
}

void min_max_global(double *series, uint64_t series_size,
                    double *_min, double *_max)
{
  double min, max;
  min_max(series, series_size, &min, &max);
  MPI_Allreduce(MPI_IN_PLACE, &min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  *_min = min;
  *_max = max;
}


void do_hist(double *series, uint64_t series_size,
             histogram_t *hist, double min, double max)
{
  uint64_t i;
  hist->aggregate = 0;
  hist->samples = series_size;
  hist->min_value = min;
  hist->max_value = max;
  min_max(series, series_size, &hist->local_min, &hist->local_max);
  hist->bin_size = (log(hist->max_value) - log(hist->min_value)) /
                                                          HIST_NUM_BINS;
  for(i = 0; i < HIST_NUM_BINS; ++i) {
    hist->count[i] = 0;
    hist->time[i] = 0;
  }

  double log_min_value = log(hist->min_value);
  for(i = 0; i < series_size; ++i) {
    int bin = (int) ((log(series[i]) - log_min_value) / hist->bin_size);
    if (bin == HIST_NUM_BINS) bin--; // kludge needed only for max_value
    hist->count[bin]++;
    hist->time[bin] += series[i];
    hist->aggregate += series[i];
  }

  hist->mean_value = ((double) hist->aggregate) / ((double) hist->samples);

  hist->sum_squares = 0;
  for (i = 0; i < series_size; ++i) {
    double dev = series[i] - hist->mean_value;
    hist->sum_squares += dev*dev;
  }

  //Perc variation: 100*sigma/mean
  hist->perc_var = 100*sqrt(hist->sum_squares/hist->mean_value)/hist->mean_value;
}

void do_hist_global(double *series, uint64_t series_size,
                    histogram_t *hist, double min, double max)
{
  do_hist(series, series_size, hist, min, max);
  hist->local_min = min;
  hist->local_max = max;

  MPI_Allreduce(MPI_IN_PLACE, &(hist->aggregate), 1 , MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, hist->count, HIST_NUM_BINS, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, hist->time, HIST_NUM_BINS, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &(hist->samples), 1 , MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &(hist->sum_squares), 1 , MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  // global average
  hist->mean_value = ((double) hist->aggregate) / ((double) hist->samples);

  //Global perc variation: 100*sigma/global-avg
  hist->perc_var = 100*sqrt(hist->sum_squares/(double)hist->samples)/hist->mean_value;
}

void print_rank_stats(char * title, histogram_t *chist) {
  printf("*************** %s Results ********************\n", title);
  print_hist(chist,"Compute and Jitter histogram");
  printf("\n");
}

void print_global_stats(char * title, histogram_t *chist, histogram_t *phist,
    struct calibration *calib, double runtime_elapsed) {
  printf("*************** %s Results ********************\n", title);
  double compute_usec = (double)(calib->num_phases*chist->min_value);
  printf("Runtime expected (compute only): %.4lf sec\n", compute_usec/1e6);
  printf("Runtime elapsed (jitter + compute + communication): %.4lf sec\n",
         runtime_elapsed);
  printf("Slowdown due to jitter and communication: %.4lfx\n",
         runtime_elapsed/(compute_usec/1e6));
  printf("Cpu time (all ranks): %.4lf sec composed as:\n",
         ((double) phist->aggregate)/1e6);
  printf("  %06.3lf%% of cpu time is compute\n",
         100*np*compute_usec/((double) phist->aggregate));
  printf("  %06.3lf%% of cpu time is jitter\n",
         100*(((double) chist->aggregate) - np*compute_usec) /
                                              ((double) phist->aggregate));
  printf("  %06.3lf%% of cpu time is communication\n",
         100*(((double) phist->aggregate) - ((double) chist->aggregate)) /
                                              ((double) phist->aggregate));
  print_hist(chist,"\nCompute and Jitter histogram");
  print_hist(phist,"\nPhase histogram (compute, jitter and communication)");
}

int try_lock(const char* lockfile)
{
    int fd = open(lockfile, O_CREAT|O_RDONLY,0666);
    if (fd == -1)
        perrorMsg("try_lock: open");

    int err = flock(fd, LOCK_EX|LOCK_NB);
    if (err == 0)
        return fd; //Aquired
    else if (errno != EWOULDBLOCK)
        perrorMsg("try_lock: flock"); //error

    close(fd);
    return 0; // Not acquired
}

void release_lock(int fd)
{
    int err = flock(fd, LOCK_UN);
    if (err)
        perrorMsg("release_lock: flock");
    close(fd);
}

// -----------------------------------------------------------------------------

#define USEC_PER_SEC    1000000
#define NSEC_PER_SEC    1000000000
#define NSEC_PER_MSEC   1000000
#define NSEC_PER_USEC   1000

void time_get(struct timespec* t)
{
    clock_gettime(timestamp_clock, t);
}

double dtime(struct timespec *t)
{
    return ((double) t->tv_sec) + (((double) t->tv_nsec) * 1e-9);
}

void time_diff(
        struct timespec* dt,
        struct timespec* t1,
        struct timespec* t0)
{
    dt->tv_sec  = t1->tv_sec  - t0->tv_sec;
    dt->tv_nsec = t1->tv_nsec - t0->tv_nsec;
    if (t1->tv_nsec < t0->tv_nsec) {
        dt->tv_sec  -= 1;
        dt->tv_nsec += NSEC_PER_SEC;
    }
}

void hline(int n, char ch)
{
    for (int i = 0; i < n; i++)
        putchar(ch);
    putchar('\n');
}

// -----------------------------------------------------------------------------

static inline void compute_lutexp(uint64_t iters)
{
    double tsum = 0.0;
    long lrand, ndx;
    uint64_t i;
    lrand = (long) nrand;
    for (i = 0L; i < iters; i++) {
        ndx = i % lrand;
        tsum = tsum + lutexp(xrand[ndx]);
    }
    static volatile double keep = 0.0;
    keep += tsum; // keep compiler from discarding loop
}

void __config_lutexp_cuda(double *xrand, int _nrand) __attribute__((weak));
void __compute_lutexp_cuda(uint64_t iters) __attribute__((weak));
static inline void compute_lutexp_cuda(uint64_t iters)
{
  __compute_lutexp_cuda(iters);
}

static inline void compute_sqrt(uint64_t iters)
{
    uint64_t i;
    double dummy = 0.0;

    for (i = 0L; i < iters; i++)
        dummy += sqrt((double) (i % 10L));
    static volatile double keep = 0.0;
    keep += dummy; // keep compiler from discarding loop
}

static inline void compute_add(uint64_t iters)
{
    uint64_t i;
    uint64_t dummy = 0;

    for (i = 0L; i < iters; i++)
        dummy += i;
    static volatile uint64_t keep = 0;
    keep += dummy; // keep compiler from discarding loop
}

// If "workload" is a compile-time constant, the switch will disappear.
static inline void compute(comp_kern_t workload, uint64_t iters)
{
    switch (workload) {
    case COMPUTE_LUTEXP:
        compute_lutexp(iters);
        break;
    case COMPUTE_LUTEXP_CUDA:
        compute_lutexp_cuda(iters);
        break;
    default:
    case COMPUTE_SQRT:
        compute_sqrt(iters);
        break;
    case COMPUTE_ADD:
        compute_add(iters);
        break;
    }
}

static inline void phase_loop(
    int communication, comp_kern_t workload,
    uint64_t phases, uint64_t iters, struct timespec *timestamp)
{
    uint64_t p;
    struct timespec *tp = timestamp;
    time_get(tp++);
    for (p = 0; p < phases; p++) {
        compute(workload, iters);
        time_get(tp++);
        if (communication != COMM_NONE) {
            if (communication == COMM_BARRIER) {
                MPI_Barrier(MPI_COMM_WORLD);
            } else {
                int r = 1;
                MPI_Allreduce(MPI_IN_PLACE, &r, 1, MPI_INT,
                        MPI_SUM, MPI_COMM_WORLD);
            }
            time_get(tp++);
        }
    }
}

// Generate multiple copies of the loop without conditionals.
static void do_phase_loop(
    int communication, comp_kern_t workload,
    uint64_t phases, uint64_t iters, struct timespec *timestamp)
{
    #define LOOP(barrier,work) \
        phase_loop(barrier, work, \
                   phases, iters, timestamp)
    switch (communication) {
    case COMM_NONE:
        switch (workload) {
        case COMPUTE_LUTEXP:      LOOP(COMM_NONE, COMPUTE_LUTEXP);      break;
        case COMPUTE_LUTEXP_CUDA: LOOP(COMM_NONE, COMPUTE_LUTEXP_CUDA); break;
        case COMPUTE_SQRT:        LOOP(COMM_NONE, COMPUTE_SQRT);        break;
        case COMPUTE_ADD:         LOOP(COMM_NONE, COMPUTE_ADD);         break;
        }
        break;
    case COMM_BARRIER:
        switch (workload) {
        case COMPUTE_LUTEXP:      LOOP(COMM_BARRIER, COMPUTE_LUTEXP);      break;
        case COMPUTE_LUTEXP_CUDA: LOOP(COMM_BARRIER, COMPUTE_LUTEXP_CUDA); break;
        case COMPUTE_SQRT:        LOOP(COMM_BARRIER, COMPUTE_SQRT);        break;
        case COMPUTE_ADD:         LOOP(COMM_BARRIER, COMPUTE_ADD);         break;
        }
        break;
    case COMM_ALLREDUCE:
        switch (workload) {
        case COMPUTE_LUTEXP:      LOOP(COMM_ALLREDUCE, COMPUTE_LUTEXP);      break;
        case COMPUTE_LUTEXP_CUDA: LOOP(COMM_ALLREDUCE, COMPUTE_LUTEXP_CUDA); break;
        case COMPUTE_SQRT:        LOOP(COMM_ALLREDUCE, COMPUTE_SQRT);        break;
        case COMPUTE_ADD:         LOOP(COMM_ALLREDUCE, COMPUTE_ADD);         break;
        }
        break;
    }
    #undef LOOP
}

double compute_best(int workload, int communication, uint64_t iters)
{
    int i;
    int timestamp_count = 1 + COMPUTE_BEST_TRIAL * (1 + !!communication);
    size_t timestamp_size = timestamp_count * sizeof(struct timespec);
    struct timespec* timestamp = malloc(timestamp_size);
    if (timestamp == NULL)
        errExit("error allocating timestamp array for calibration\n");
    mlock(timestamp, timestamp_size);
    memset(timestamp, 0, timestamp_size);
    do_phase_loop(communication, workload,
                  COMPUTE_BEST_TRIAL, iters, timestamp);
    struct timespec *tp = timestamp;
    double best = 1000000000.0;
    for (i = 0; i < COMPUTE_BEST_TRIAL; i++) {
        double duration = dtime(tp+1) - dtime(tp);
        if (duration < best) best = duration;
        tp += communication ? 2 : 1;
    }
    free(timestamp);
    return best * USEC_PER_SEC;
}

void set_node_info()
{
    gethostname(node_info.hostname, MAX_HOSTNAME_CHAR);
    node_info.cpu = get_cpu();
    MPI_Gather(&node_info, sizeof(node_info_t), MPI_BYTE, node_info_array, sizeof(node_info_t),
            MPI_BYTE, 0, MPI_COMM_WORLD);
}

void calibrate(int workload, int communication, struct calibration *calib)
{
    double calib_time = time_sec();
    TMrand(seed, nrand, xrand);
    if(__config_lutexp_cuda)
      __config_lutexp_cuda(xrand, nrand);

    MPI_Barrier(MPI_COMM_WORLD);
    if (calib->iters == 0) {
        double x1 = 1000, x2 = 10000, y1, y2, m, b;
        volatile double dummy = compute_best(workload, communication, x1); (void) dummy;// Warmup
        y1 = compute_best(workload, communication, x1);  // First point
        // Take min among ranks
        MPI_Allreduce(MPI_IN_PLACE, &y1, 1, MPI_DOUBLE,
                MPI_MIN, MPI_COMM_WORLD);
        y2 = compute_best(workload, communication, x2);  // Second point
        // Take min among ranks
        MPI_Allreduce(MPI_IN_PLACE, &y2, 1, MPI_DOUBLE,
                MPI_MIN, MPI_COMM_WORLD);
        // Slope
        m = (y2 - y1)/(x2 - x1);
        // Intercept
        b = y1 - m * x1;
        // Compute iterations with slope and intercept
        calib->iters = (uint64_t) ((calib->compute_usec - b) / m); // Iterations
        /*printf("x1 %lf  y1 %lf  x2 %lf  y2 %lf  m %lf  b %lf compute_usec %lf  iters %lu\n",
                x1, y1, x2, y2, m, b, calib->compute_usec, calib->iters);*/
    } else {
        // Iterations set by the argument
        calib->compute_usec = compute_best(workload, communication, calib->iters);
        // Take the min elapsed time among ranks
        MPI_Allreduce(MPI_IN_PLACE, &calib->compute_usec, 1, MPI_DOUBLE,
                MPI_MIN, MPI_COMM_WORLD);
    }

    if (rank == 0) {
        printf("******************** Calibration ******************\n");
        printf("Compute phase expected time: %.3lf msec\n", calib->compute_usec/1e3);
        printf("Compute phase iterations: %lu\n", (long)calib->iters);
        printf("Calibration time: %3lf sec\n",time_sec() - calib_time);
    }
}

int main(int argc, char *argv[])
{
    // initialize mpi
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    // parameters and defaults
    int cpu_stride     = -1;
    int ranks_per_node =  0;
    int communication  =  COMM_BARRIER;

    int  print_raw     =  0;        // print raw data
    long print_cutoff  =  LONG_MAX; // raw data phase threshold
    int  print_binary  =  0;        // print raw data in binary format
    int  tracing       =  0;        // enable filtered tracing
    int  print_extra   =  0;        // print extra statistics
    comp_kern_t workload = COMPUTE_SQRT;
    char *compute_math =  "sqrt";

    const char *payload_elf = NULL;
    const char **payload_argv = NULL;

    // calibration structure
    struct calibration calib;
    calib.num_phases       = 1000;                // number of compute+barrier phases
    calib.compute_usec     = 500.0;              // compute phase size in microseconds
    calib.iters      = 0;               //number of iterations not set (do calibration)

    // parse arguments
    int c, found_payload = 0;
    while ( ! found_payload
            && (c = getopt (argc, argv, "l:s:p:c:nfr:C:d:beti:j:w:b:m:ha:")) != -1){
        switch (c) {
        case 'l':
            ranks_per_node = atoi(optarg);
            break;
        case 's':
            cpu_stride = atoi(optarg);
            break;
        case 'p':
            calib.num_phases = atol(optarg);
            break;
        case 'c':
            calib.compute_usec = atol(optarg);
            break;
        case 'n':
            per_rank_stats = 1;
            break;
        case 'r':
            switch (optarg[0]) {
            case 'n': communication = COMM_NONE; break;
            case 'b': communication = COMM_BARRIER; break;
            case 'a': communication = COMM_ALLREDUCE; break;
            default: errExit("Unknown communication (-r) option\n");
            }
            break;
        case 'C':
            if (strcmp(optarg, "mono") == 0)
                timestamp_clock = CLOCK_MONOTONIC;
            else if (strcmp(optarg, "real") == 0)
                timestamp_clock = CLOCK_REALTIME;
            else
                errExit("Unknown -C (timestamp_clock) value\n");
            break;
        case 'd':
            print_raw = 1;
            print_cutoff = atol(optarg);
            break;
        case 'b':
            print_raw = 1;
            print_binary = 1;
            break;
        case 't':
            tracing = 1;
            break;
        case 'e':
            print_extra = 1;
            break;
        case 'i':
            // use fixed number of iterations
            calib.iters = atol(optarg);
           break;
        case 'm':
            compute_math = optarg;
            if (strcmp(optarg, "lutexp") == 0)
                workload = COMPUTE_LUTEXP;
            else if (strcmp(optarg, "lutexp_cuda") == 0
                     && __compute_lutexp_cuda
                     && __config_lutexp_cuda)
                workload = COMPUTE_LUTEXP_CUDA;
            else if (strcmp(optarg, "sqrt") == 0)
                workload = COMPUTE_SQRT;
            else if (strcmp(optarg, "add") == 0)
                workload = COMPUTE_ADD;
            else
                errExit("Unknown -m value '%s'\n", optarg);
            break;
        case 'a': // '--args' such as in gdb
            if(access(optarg, F_OK) != -1) {
                payload_elf = optarg;
                payload_argv = (const char**) malloc (argc * sizeof(const char*));
                int i;
                for(i = 0; i < argc; ++i) {
                  if(argv[i] == payload_elf) {
                    int new_argc = 0;
                    for(++i; i < argc; ++i, ++new_argc)
                      payload_argv[new_argc] = argv[i];
                    payload_argv[new_argc] = NULL;
                  }
                }
                // exit getopt parsing as all subsequent parameters are for payload
                found_payload = 1;
            }
            else
                errExit("ERR: provided binary '%s' can not be accessed!\n", optarg);
            break;
        case 'h':
        default:
            usage(argv[0]);
            MPI_Finalize();
            return 0;
        }
    }

    // get node information
    node_info_array = (node_info_t*) malloc(np*sizeof(node_info_t));
    set_node_info();

    if (rank == 0 && cpu_stride > ranks_per_node) {
        errExit("Total cpus or cpu stride error\n");
    }

    // -------------------------------------------------------------------------
    // process tracing
    // -------------------------------------------------------------------------

    char* tracing_dir        = NULL;
    char* tracing_job        = NULL;
    char  tracing_fpath[256] = {0};
    FILE* tracing_fp         = NULL;

    if (!rank && tracing) { // rank 0 and tracing enabled

        // read tracing directory
        tracing_dir = getenv("JTF_DIR");
        if (!tracing_dir) {
            printf("[jitter-bench] tracing: missing JTF_DIR\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        // read tracing job id
        tracing_job = getenv("JTF_JOB");
        if (!tracing_job) {
            printf("[jitter-bench] tracing: missing JTF_JOB\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        // prepare JTF filepath string
        sprintf(tracing_fpath, "%s/job%s.jtf", tracing_dir, tracing_job);

        // open JTF file
        tracing_fp = fopen(tracing_fpath, "w");
        if (!tracing_fp) {
            printf("[jitter-bench] tracing: could not create JTF file: %s\n"
                    "error %d: %s\n",
                    tracing_fpath, errno, strerror(errno));
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
    }

    // -------------------------------------------------------------------------
    // greetings
    // -------------------------------------------------------------------------

    if (rank == 0) {
        printf("************************** Jitter bench *************************\n"
                "Num-ranks: %d Num-phases: %lu Compute-math:%s\n",
                np, calib.num_phases, compute_math);
        if ( cpu_stride > 0)
            printf("Ranks-per-node: %d Cpu-stride: %d\n",ranks_per_node, cpu_stride);

        char* clock_str;
        switch(timestamp_clock) {
        case CLOCK_MONOTONIC: clock_str = "monotonic"; break;
        case CLOCK_REALTIME: clock_str = "realtime"; break;
        default: clock_str = "unknown"; break;
        }
        char* comm_str;
        switch(communication) {
        case COMM_NONE: comm_str = "none"; break;
        case COMM_BARRIER: comm_str = "barrier"; break;
        case COMM_ALLREDUCE: comm_str = "all-reduce"; break;
        default: comm_str = "unknown"; break;
        }
        printf("clock: %s  communication: %s\n", clock_str, comm_str);
        printf("print: raw: %d  binary: %d  cutoff: %ld  extra: %d\n",
                print_raw, print_binary, print_cutoff, print_extra);
        printf("tracing: %d dir: %s job: %s fpath: %s\n", 
                tracing, tracing_dir, tracing_job, tracing_fpath);
    }

    // pin cpu
    if(cpu_stride > 0)
        set_cpu((rank-(rank/ranks_per_node)*ranks_per_node)*cpu_stride);


    // calibration
    calibrate(workload, communication, &calib);


    // We need space for an initial timestamp plus one or two timestamps
    // per phase, depending on whether we're doing communication.
    register uint64_t phases = calib.num_phases;
    int timestamp_count = 1 + (phases * ((communication == COMM_NONE) ? 1 : 2));
    size_t timestamp_size = timestamp_count * sizeof(struct timespec);
    struct timespec* timestamp = malloc(timestamp_size);
    if (timestamp == NULL) {
        printf("rank %d: error allocating timestamp array: (%d) %s\n",
                rank, errno, strerror(errno));
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    mlock(timestamp, timestamp_size);

    // -------------------------------------------
    // benchmark loop
    // -------------------------------------------
    //
    // warmup and synchronize before starting loop
    do_phase_loop(communication, workload,
                  calib.num_phases / 4, calib.iters, timestamp);
    memset(timestamp, 0, timestamp_size);

    MPI_Barrier(MPI_COMM_WORLD);
    if(payload_elf)
        musl_execve(payload_elf, payload_argv);
    MPI_Pcontrol(1);
    double runtime_elapsed = time_sec();

    do_phase_loop(communication, workload,
                  calib.num_phases, calib.iters, timestamp);

    runtime_elapsed = time_sec() - runtime_elapsed;
    MPI_Pcontrol(0);
    MPI_Barrier(MPI_COMM_WORLD);

    // -------------------------------------------
    // update legacy statistics code
    // -------------------------------------------

    // prepare old data structures (keep previous statistics code working)
    double *cseries = calloc(calib.num_phases,sizeof(double)); // compute
    double *pseries = calloc(calib.num_phases,sizeof(double)); // phase (jitter+compute+comm)

    struct timespec* tp = timestamp; // auxiliary pointer to timestamp array
    struct timespec* tap = tp++; // auxiliary pointer to timestamp structure A
    for (uint64_t i = 0; i < phases; i++) {

      struct timespec* tbp = tp++;   // auxiliary pointer to timestamp structure B
      struct timespec  ds;           // time difference between A and B (structure)
      double           dt;           // time difference between A and B (double)
      struct timespec  bs;           // time difference between B and C (structure)
      double            b = 0.0;     // barrier (double)


      // calculate time difference and determine jitter
      time_diff(&ds, tbp, tap);
      dt = (double) ds.tv_sec + (double) ds.tv_nsec / (double) NSEC_PER_SEC;

      // calculate comm duration
      if (communication != COMM_NONE) {
        struct timespec* tcp = tp++; // auxiliary pointer to timestamp structure C
        time_diff(&bs, tcp, tbp);
        b = (double) bs.tv_sec + (double) bs.tv_nsec / (double) NSEC_PER_SEC;
      }

      // update previous structures
      cseries[i] = dt * USEC_PER_SEC;
      pseries[i] = (dt + b) * USEC_PER_SEC;

      tap = tp - 1;
    }

    // -------------------------------------------
    // back to previous statistics
    // -------------------------------------------

    // Global jitter histogram
    double cmin, cmax, pmin, pmax;
    histogram_t chist_glob, phist_glob;
    min_max_global(cseries, calib.num_phases, &cmin, &cmax);
    min_max_global(pseries, calib.num_phases, &pmin, &pmax);
    do_hist_global(cseries, calib.num_phases, &chist_glob, cmin, cmax);
    do_hist_global(pseries, calib.num_phases, &phist_glob, pmin, pmax);

    // expected compute phase duration (double)
    double e = cmin / (double) USEC_PER_SEC;

    // Print global stats
    if (rank == 0) {
      print_global_stats("Global", &chist_glob, &phist_glob, &calib, runtime_elapsed);
    }

    // Compute per-node jitter histograms
    if (per_rank_stats) {
      // Recompute local histograms (different min,max than global)
      histogram_t chist_local, phist_local;
      do_hist(pseries, calib.num_phases, &phist_local, pmin, pmax);
      do_hist(cseries, calib.num_phases, &chist_local, cmin, cmax);


      // Gather histograms from all ranks
      histogram_t *chist_array = NULL, *phist_array = NULL;
      double* sort_array = NULL;
      int* nodes_sorted = NULL;
      if (rank == 0) {
        chist_array = malloc(np*sizeof(histogram_t));
        phist_array = malloc(np*sizeof(histogram_t));
        sort_array = malloc(np * sizeof(double));
        nodes_sorted = malloc(np * sizeof(int));
      }

      //Gather all stats
      MPI_Gather(&chist_local, sizeof(histogram_t), MPI_BYTE, chist_array, sizeof(histogram_t),
          MPI_BYTE, 0, MPI_COMM_WORLD);
      MPI_Gather(&phist_local, sizeof(histogram_t), MPI_BYTE, phist_array, sizeof(histogram_t),
          MPI_BYTE, 0, MPI_COMM_WORLD);

      if (rank == 0) {
        int n;
        for ( n = 0; n < np; n++) {
          // sort based on aggregated compute + jitter time
          sort_array[n] = ((double) chist_array[n].aggregate);
        }

        sortx(sort_array, np, nodes_sorted, SORT_DESCENDING_ORDER);
        // Print sorted stats
        for (n = 0; n < np; ++n) {
          int r=nodes_sorted[n];
          char title[256];
          sprintf(title, "Rank %d [%s cpu %d]", r,
              node_info_array[r].hostname, node_info_array[r].cpu);
          print_rank_stats(title, &chist_array[r]);
        }
        free(chist_array);
        free(phist_array);
        free(sort_array);
        free(nodes_sorted);
      }
    }

    // prepare binary file if needed
    FILE* bin_f = NULL;
    if (print_binary) {
        char  bin_fn[512] = {0};
        char* env;

        if ((env = getenv("LSB_OUTPUTFILE")))   sprintf(bin_fn, "%s.jbdump", env);
        else if ((env = getenv("LSB_JOBID")))   sprintf(bin_fn, "job%s.out.jbdump", env);
        else                                    sprintf(bin_fn, "job.out.jbdump");

        int total_ranks = np;
        int total_phases = phases;
        double exp = e;
        bin_f = fopen(bin_fn, "w");
        fwrite(&total_ranks, sizeof(int), 1, bin_f);
        fwrite(&total_phases, sizeof(int), 1, bin_f);
        fwrite(&exp, sizeof(double), 1, bin_f);
    }

    // print raw data (per rank)
    if (print_raw || print_extra) {
        struct timespec* timestamp_all = malloc(np * timestamp_size);

        MPI_Gather(timestamp, timestamp_size, MPI_BYTE,
                timestamp_all, timestamp_size, MPI_BYTE, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            printf("\n");
            printf("EXTRA STATISTICS\n");
            hline(150, '=');

            double TOT_C = 0.0;
            double TOT_J = 0.0;
            double TOT_B = 0.0;

            int    TOT_discr_cnt = 0;
            double TOT_discr_sum = 0.0;

            for (int r = 0; r < np; r++) {
                int header_needed = 1;

                // rank specific statistics
                double   jt = 0.0; // total jitter
                double   wj = 0.0; // worst jitter
                uint64_t wi =  -1; // worst jitter index
                double   bt = 0.0; // total barrier
                double   wb = 0.0; // worst barrier
                uint64_t bi =  -1; // worst barrier index

                struct timespec null_ts = {0, 0};
                struct timespec* tcp = &null_ts; // null for no-barrier case
                struct timespec* tp = &timestamp_all[r * timestamp_count];
                struct timespec* tap = tp++;
                for (uint64_t i = 0; i < phases; i++) {
                    struct timespec* tbp = tp++;
                    struct timespec  ds; // time difference between A and B (structure)
                    double           dt; // time difference between A and B (double)
                    double           j;  // jitter (compute duration - expected)
                    struct timespec  bs; // time difference between B and C (structure)
                    double           b;  // time difference between B and C (double)

                    // calculate compute phase duration
                    time_diff(&ds, tbp, tap);
                    dt = (double) ds.tv_sec + (double) ds.tv_nsec / (double) NSEC_PER_SEC;
                    j  = dt - e;
                    if (j < 0.0) {
                        TOT_discr_cnt += 1;
                        TOT_discr_sum += fabs(j);
                        j = 0.0;
                    }
                    jt += j;
                    TOT_C += e;
                    TOT_J += j;

                    // update worst jitter
                    if (j > wj) {
                        wj = j;
                        wi = i;
                    }

                    if (communication != COMM_NONE) {
                        // calculate communication duration
                        tcp = tp++;
                        time_diff(&bs, tcp, tbp);
                        b = (double) bs.tv_sec + (double) bs.tv_nsec /
                            (double) NSEC_PER_SEC;
                        bt += b;
                        TOT_B += b;

                        // update worst barrier
                        if (b > wb) {
                            wb = b;
                            bi = i;
                        }
                    }

                    if (print_raw && print_binary) {
                        /*
                        fwrite(tap, sizeof(struct timespec), 1, bin_f);
                        fwrite(tbp, sizeof(struct timespec), 1, bin_f);
                        fwrite(tcp, sizeof(struct timespec), 1, bin_f);
                        */
                        int tap_s  = tap->tv_sec;
                        int tap_ns = tap->tv_nsec;
                        int tbp_s  = tbp->tv_sec;
                        int tbp_ns = tbp->tv_nsec;
                        fwrite(&tap_s,  sizeof(int), 1, bin_f);
                        fwrite(&tap_ns, sizeof(int), 1, bin_f);
                        fwrite(&tbp_s,  sizeof(int), 1, bin_f);
                        fwrite(&tbp_ns, sizeof(int), 1, bin_f);
                        if (r == rank && i == phases - 1) {
                                int tcp_s = tcp->tv_sec;
                                int tcp_ns = tcp->tv_nsec;
                                fwrite(&tcp_s,  sizeof(int), 1, bin_f);
                                fwrite(&tcp_ns, sizeof(int), 1, bin_f);
                        }
                    }

                    else if (print_raw && ((dt * 1e6) >= ((double) print_cutoff))) {
                        if (header_needed) {
                            printf("\n");
                            hline(144, '-');
                            printf("%5s %10s %20s %20s %20s "
                                   "%16s %16s %7s %16s\n",
                                   "rank", "phase", "ta", "tb", "tc",
                                   "duration", "jitter", "%", "barrier");
                                    hline(144, '-');
                            header_needed = 0;
                        }
                        printf("%5d %10lu %10ld.%.09ld %10ld.%.09ld %10ld.%.09ld "
                               "%16.9lf %16.9lf %6.2f%% ",
                               r,
                               i,
                               tap->tv_sec, tap->tv_nsec,
                               tbp->tv_sec, tbp->tv_nsec,
                               tcp->tv_sec, tcp->tv_nsec,
                               dt,
                               j,
                               100.0 * j / e);

                        if (communication != COMM_NONE) {
                            printf("%16.9lf", b);
                        }
                        printf("\n");

                        // JTF <start>
                        if (tracing) {
                            fprintf(tracing_fp, 
                                    "%-10d %s %-3d %10ld.%06ld %10ld.%06ld %16.6lf\n",
                                    r,
                                    node_info_array[r].hostname, 
                                    node_info_array[r].cpu,
                                    tap->tv_sec, tap->tv_nsec / 1000,
                                    tbp->tv_sec, tbp->tv_nsec / 1000,
                                    dt);
                        }
                        // JTF <end>
                    }

                    tap = tp - 1;
                }

                if (print_extra) {
                    // calculate and print additional statistics
                    struct timespec* ti = &timestamp_all[r * timestamp_count];
                    struct timespec* tf = ti + timestamp_count - 1;
                    struct timespec  tot;
                    time_diff(&tot, tf, ti);

                    // worst compute interval start and end
                    struct timespec* tjap = ti + (wi * ((communication == COMM_NONE) ? 1 : 2));
                    struct timespec* tjbp = tjap + 1;

                    printf("\n");
                    printf("(rank %5d) First time:     %6ld.%09ld\n",
                            r, ti->tv_sec, ti->tv_nsec);
                    printf("(rank %5d) Last time:      %6ld.%09ld\n",
                            r,tf->tv_sec, tf->tv_nsec);
                    printf("(rank %5d) Total time:     %16.9lf\n",
                            r, (double) tot.tv_sec + (double) tot.tv_nsec / (double) NSEC_PER_SEC);
                    printf("(rank %5d) Total jitter:   %16.9lf (%6.2lf%%)\n",
                            r, jt, 100.0 * jt / (phases * e));
                    printf("(rank %5d) Avg jitter:     %16.9lf (%6.2lf%%)\n",
                            r, jt / phases, 100.0 * jt / phases / e);
                    printf("(rank %5d) Worst jitter:   %16.9lf (%6.2lf%%) @ phase %-10lu (%6ld.%09ld to %6ld.%09ld)\n",
                            r, wj, 100.0 * wj / e, wi,
                            tjap->tv_sec, tjap->tv_nsec,
                            tjbp->tv_sec, tjbp->tv_nsec);
                    if (communication != COMM_NONE) {
                        // worst barrier start and end
                        struct timespec* tbbp = ti + 1 + (bi * 2);
                        struct timespec* tbcp = tbbp + 1;
                        printf("(rank %5d) Total barrier:  %16.9lf\n",
                                r, bt);
                        printf("(rank %5d) Avg barrier:    %16.9lf\n",
                                r, bt / phases);
                        printf("(rank %5d) Worst barrier:  %16.9lf %9s @ phase %-10lu (%6ld.%09ld to %6ld.%09ld)\n",
                                r, wb, "", bi,
                                tbbp->tv_sec, tbbp->tv_nsec,
                                tbcp->tv_sec, tbcp->tv_nsec);
                    }
                    printf("\n");
                }
            }

            double TOT = TOT_C + TOT_J + TOT_B;         // total x ranks
            double EXP = (double) calib.num_phases * e; // expected total duration

            printf("\n");
            printf("AGGREGATE STATISTICS\n");
            hline(80, '=');
            printf("\n");
            printf("RANKS:         %20d\n", np);
            printf("TOTAL:         %20.6lf\n", TOT);
            printf("TOTAL COMPUTE: %20.6lf (%6.2lf%%)\n", TOT_C, 100.0 * TOT_C / TOT);
            printf("TOTAL JITTER:  %20.6lf (%6.2lf%%)\n", TOT_J, 100.0 * TOT_J / TOT);
            printf("TOTAL BARRIER: %20.6lf (%6.2lf%%)\n", TOT_B, 100.0 * TOT_B / TOT);
            printf("SLOWDOWN:      %20.6lf (%6.2lf%%)\n", runtime_elapsed, 100.0 * (runtime_elapsed - EXP) / EXP);
            hline(80, '=');
            printf("TOT_discr_cnt: %20d\n",    TOT_discr_cnt);
            printf("TOT_discr_sum: %20.6lf\n", TOT_discr_sum);
            printf("\n");
        }
        free(timestamp_all);
    }

    // cleanup and finish
    free(timestamp);
    if(payload_argv)
      free(payload_argv);
    free(cseries);
    free(pseries);
    free(node_info_array);
    MPI_Finalize();
    return 0;
}

//////////////////////////////////////////////////////
//
// https://www.musl-libc.org/
// musl/src/process/system.c
//
// modified POSIX system() call
// to behave like an execve that is returning
//
// MIT License
//
//////////////////////////////////////////////////////

#include <spawn.h>
#include <pthread.h>
int musl_execve(const char *payload_elf,
                const char **payload_argv)
{
    pid_t pid;
    sigset_t old, reset;
    struct sigaction sa = { .sa_handler = SIG_IGN }, oldint, oldquit;
    int status = 0x7f00, ret;
    posix_spawnattr_t attr;

    pthread_testcancel();

    if (!payload_elf) return 1;

    sigaction(SIGINT, &sa, &oldint);
    sigaction(SIGQUIT, &sa, &oldquit);
    sigaddset(&sa.sa_mask, SIGCHLD);
    sigprocmask(SIG_BLOCK, &sa.sa_mask, &old);

    sigemptyset(&reset);
    if (oldint.sa_handler != SIG_IGN) sigaddset(&reset, SIGINT);
    if (oldquit.sa_handler != SIG_IGN) sigaddset(&reset, SIGQUIT);
    posix_spawnattr_init(&attr);
    posix_spawnattr_setsigmask(&attr, &old);
    posix_spawnattr_setsigdefault(&attr, &reset);
    posix_spawnattr_setflags(&attr, POSIX_SPAWN_SETSIGDEF|POSIX_SPAWN_SETSIGMASK);
    ret = posix_spawn(&pid, payload_elf, 0, &attr, (char * const*)payload_argv, __environ);
    posix_spawnattr_destroy(&attr);

    if (!ret) while (waitpid(pid, &status, 0)<0 && errno == EINTR);
    sigaction(SIGINT, &oldint, NULL);
    sigaction(SIGQUIT, &oldquit, NULL);
    sigprocmask(SIG_SETMASK, &old, NULL);

    if (ret) errno = ret;
    return status;
}
