#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define HOST  512 // Hostname length
#define LINE  512 // Line parsing buffer size
#define HLINE 80  // Length of horizontal line

//
//  JTF region
//

typedef struct {
  int rank;         // rank
  char host[HOST];  // hostname
  int cpu;          // region cpu
  double ta;        // region start timestamp
  double tb;        // region end timestamp
  double dt;        // region time delta
} reg_t;

//
//  Program variables
//

typedef struct {
  char* jtf;    // path to JTF file
  char* src;    // path to source trace (raw)
  char* dst;    // path to destination trace (filtered)

  FILE* jtf_f;
  FILE* src_f;
  FILE* dst_f;

  reg_t* regs;  // list of JTF regions
  int regc;     // number of JTF regions
  char* host;   // host to be processed
  int cpu;      // cpu to be processed
} var_t;

//
//  Function prototypes
//

void parse(var_t* var, int argc, char** argv);
void usage(char* exec);
void prepare_files(var_t* var);
void error_fopen(char* filepath);
void process_jtf(var_t* var);
void process_trace(var_t* var);
double advance(var_t* var, char* buffer);
double parse_time(char* timestamp);
void copy(var_t*, char* buffer);
void hline(FILE* fp);

//
//  Program entry point
//

int main(
    int argc,
    char** argv)
{
  var_t var;
  parse(&var, argc, argv);
  prepare_files(&var);
  process_jtf(&var);
  process_trace(&var);
  return 0;
}

//
//  Parse command line arguments
//

void parse(
    var_t* var,
    int argc,
    char** argv)
{
  if (argc <= 6)
    usage(argv[0]);

  var->jtf  = argv[1];
  var->regc = atoi(argv[2]);
  var->host = argv[3];
  var->cpu  = atoi(argv[4]);
  var->src  = argv[5];
  var->dst  = argv[6];
}

//
//  Print usage information and exit
//

void usage(
    char* exec)
{
  printf("usage: %s jtf regc host cpu src dst\n"
         "  jtf   JTF auxiliary file\n"
         "  regc  number of JTF regions (entries)\n"
         "  host  host to be processed"
         "  cpu   cpu to be processed\n"
         "  src   source trace\n"
         "  dst   destination trace\n",
         exec);
  exit(1);
}

//
//  Prepare trace files
//

void prepare_files(
    var_t* var)
{
  var->jtf_f = fopen(var->jtf, "r");
  if (!var->jtf_f)
    error_fopen(var->jtf);

  var->src_f = fopen(var->src, "r");
  if (!var->src_f)
    error_fopen(var->src);

  var->dst_f = fopen(var->dst, "w");
  if (!var->dst_f)
    error_fopen(var->dst);
}

//
//  Handle fopen errors
//

void error_fopen(
    char* filepath)
{
  printf("[tracefs-filter] error: fopen(%s): %d: %s\n",
      filepath, errno, strerror(errno));
  exit(1);
}

//
//  Process regions from JTF file
//

void process_jtf(
    var_t* var)
{
  // Allocate JTF list
  var->regs = malloc(var->regc * sizeof(reg_t));

  // Read JTF
  char buffer[LINE+1];
  for (int i = 0; i < var->regc; i++) {

    // Read line
    char* tmp = fgets(buffer, LINE, var->jtf_f);
    if (!tmp) {
      printf("[tracefs-filter] error: process_jtf(): fgets(): %d: %s\n",
          errno, strerror(errno));
      exit(1);
    }

    // Add JTF entry
    sscanf(buffer, "%d %s %d %lf %lf %lf\n",
        &var->regs[i].rank,
         var->regs[i].host,
        &var->regs[i].cpu,
        &var->regs[i].ta,
        &var->regs[i].tb,
        &var->regs[i].dt);
  }
}

//
//  Process trace file filtering entries based on JTF regions list
//

void process_trace(
    var_t* var)
{
  char buffer[LINE+1] = {0};
  double time = advance(var, buffer);

  // Iterate over JTF regions list
  for (int i = 0; i < var->regc; i++) {

    // Ignore entries from other hosts
    if (strcmp(var->regs[i].host, var->host))
      continue;

    // Ignore entries from other cpus
    if (var->regs[i].cpu != var->cpu)
      continue;

    // Write header
    hline(var->dst_f);
    fprintf(var->dst_f, "%-10d %s %-3d %16.6lf %16.6lf %16.6lf\n",
        var->regs[i].rank,
        var->regs[i].host,
        var->regs[i].cpu,
        var->regs[i].ta,
        var->regs[i].tb,
        var->regs[i].dt);
    hline(var->dst_f);
    fputs("\n", var->dst_f);

    // Advance until we are inside    
    while (time < var->regs[i].ta) {
      time = advance(var, buffer);
    }

    // Copy trace entries while inside
    while (var->regs[i].ta <= time && time <= var->regs[i].tb) {      
      copy(var, buffer);
      time = advance(var, buffer);
    }

    // Skip line before moving to next JTF
    fputc('\n', var->dst_f);
  }
}

//
//  Advance buffer and time (to a valid timestamp)
//

double advance(
    var_t* var,
    char* buffer)
{
  while (1) {

    // Read line
    char* tmp = fgets(buffer, LINE, var->src_f);
    if (!tmp) {
      if (!feof(var->src_f)) {
        printf("[tracefs-filter] error: advance(): fgets(): %d: %s\n",
            errno, strerror(errno));
        exit(1);
      }
      break;
    }

    // Ignore comments
    if (buffer[0] == '#')
      continue;

    // Read timestamp and return
    char timestamp[LINE];
    sscanf(buffer, "%*[^]]] %*s %[^:]", timestamp);
    double time = parse_time(timestamp);
    return time;
  }

  // Return a large number (probably reached end of file)
  return INFINITY;
}

//
//  Parse timestamp into floating point number
//

double parse_time(
    char* timestamp)
{
  char* s  = calloc(16, sizeof(char));
  char* us = calloc(16, sizeof(char));
  sscanf(timestamp, "%[0-9].%6[0-9]", s, us);

  if (!*s || !*us) {
    printf("[tracefs-filter] error: parse_time(): malformed timestamp: %s\n",
        timestamp);
    exit(1);
  }

  for (int i = 0; i < 6; i++) 
    if (!us[i]) 
      us[i] = '0';  

  return (double) atoll(s) + (double) atoll(us) / 1.0E+6;
}

//
//  Copy trace entry to destination
//

void copy(
    var_t* var,
    char* buffer)
{
  fputs(buffer, var->dst_f);
  if (ferror(var->dst_f)) {
    printf("[tracefs-filter] error: copy(): fputs(): %d: %s\n",
        errno, strerror(errno));
    exit(1);    
  }
}

//
//  Write horizontal line to file
//

void hline(
    FILE* fp)
{
  for (int i = 0; i < HLINE; i++)
    fputc('=', fp);
  fputc('\n', fp);
}

