#define _GNU_SOURCE

#include <ctype.h>     // isdigit
#include <dirent.h>    // opendir, closedir
#include <errno.h>     // errno
#include <signal.h>    // signal, __sighandler_t, SIGUSR1, SIGUSR2
#include <stdarg.h>    // va_list
#include <stdbool.h>   // bool, true, false
#include <stdio.h>     // printf, asprintf, sprintf
#include <stdlib.h>    // exit, system
#include <string.h>    // strerror, strcmp
#include <sys/types.h> // DIR
#include <unistd.h>    // sleep, gethostname

#include "pidfile.c"   // writepid

// -----------------------------------------------------------------------------

volatile bool g_jt_started = false;
volatile bool g_jt_stopped = false;

char* g_dir;
char* g_events;
char* g_buffer;
char* g_job;
char* g_cpus;
char* g_outfile;

char* g_host;
char* g_perfdat;
char* g_perfout;

#define HOSTNAME 1024 // length of hostname
#define CMD      1024 // length of command
#define SIZE     1024 // dat file size (bytes)
#define LEN      1024 // output file name size
#define SLEEP    60   // sleep for infinite loop (seconds)

// -----------------------------------------------------------------------------

void parse_args(int argc, char** argv);
void proc_dir(char* str);
void proc_events(char* str);
void proc_buffer(char* str);
void proc_job(char* str);
void proc_cpus(char* str);
void proc_hostname();
void proc_perfiles();
void proc_outfile();

void regsig(int signum, __sighandler_t action);
void loop();

void atfprintf(char* fmt, ...);

void handler_start(int signum);
void handler_stop(int signum);

// -----------------------------------------------------------------------------

void parse_args(int argc, char** argv)
{
    if (argc < 6) {
        if (argc > 1)
            if (strcmp(argv[1], "-h") != 0 &&
                strcmp(argv[1], "--help") != 0)
                printf("error: missing arguments\n");
        printf("usage: %s dir events buffer job cpus\n",
                argv[0]);
        printf("  dir: directory to write perf outputs\n");
        printf("  events: events to be captured "
                    "(special options: \"all\", \"reduced\")\n");
        printf("  buffer: buffer size "
                    "(e.g., \"2G\")\n");
        printf("  job: job id\n");
        printf("  cpus: filter cpus to be traced "
                    "(e.g., \"off\", \"0-83,88-171\")\n");
        exit(1);
    }

    proc_dir(argv[1]);
    proc_events(argv[2]);
    proc_buffer(argv[3]);
    proc_job(argv[4]);
    proc_cpus(argv[5]);

    proc_hostname();
    proc_perfiles();
    proc_outfile();

    printf("[jt] dir: \"%s\"\n", g_dir);
    printf("[jt] events: \"%s\"\n", g_events);
    printf("[jt] buffer: \"%s\"\n", g_buffer);
    printf("[jt] job: \"%s\"\n", g_job);
    printf("[jt] cpus: \"%s\"\n", g_cpus);
    printf("[jt] hostname: \"%s\"\n", g_host);
    printf("[jt] perfdat: \"%s\"\n", g_perfdat);
    printf("[jt] perfout: \"%s\"\n", g_perfout);
    printf("[jt] outfile: \"%s\"\n", g_outfile);
    printf("\n");
}

void proc_dir(char* str)
{
   DIR* dir = opendir(str);
   if (!dir) {
       printf("[jt] error: directory \"%s\": (%d) %s\n",
               str, errno, strerror(errno));
       exit(1);
   }

   closedir(dir);
   g_dir = str;
}

void proc_events(char* str)
{
    if (strcmp(str, "all") == 0) {
        g_events = "-e "
                   "alarmtimer:*,"
                   "block:*,"
                   "bpf:*,"
                   "cgroup:*,"
                   "cma:*,"
                   "compaction:*,"
                   "context_tracking:*,"
                   "cpuhp:*,"
                   "cxl:*,"
                   "devlink:*,"
                   "dma_fence:*,"
                   "drm:*,"
                   "ext4:*,"
                   "fib6:*,"
                   "fib:*,"
                   "filelock:*,"
                   "filemap:*,"
                   "huge_memory:*,"
                   "i2c:*,"
                   "iommu:*,"
                   "irq:*,"
                   "jbd2:*,"
                   "kmem:*,"
                   "mdio:*,"
                   "migrate:*,"
                   "module:*,"
                   "napi:*,"
                   "net:*,"
                   "nfs:*,"
                   "oom:*,"
                   "page_isolation:*,"
                   "pagemap:*,"
                   "power:*,"
                   "powerpc:*,"
                   "printk:*,"
                   "random:*,"
                   "rcu:*,"
                   "regmap:*,"
                   "rpm:*,"
                   "sched:*,"
                   "scsi:*,"
                   "signal:*,"
                   "skb:*,"
                   "sock:*,"
                   "task:*,"
                   "thermal:*,"
                   "thp:*,"
                   "timer:*,"
                   "udp:*,"
                   "vmscan:*,"
                   "workqueue:*,"
                   "writeback:*,"
                   "xdp:*,"
                   "syscalls:sys_enter_sched_getscheduler,"
                   "syscalls:sys_enter_getuid";
    }
    else if (strcmp(str, "reduced") == 0) {
        g_events = "-e sched:*,irq:*,workqueue:*,syscalls:sys_enter_sched_getscheduler";
    }
    else {
        asprintf(&g_events, "-e %s", str);
    }
}

void proc_buffer(char* str)
{
    if (strcmp(str, "") == 0) {
        printf("error: empty buffer option\n");
        exit(1);
    }
    asprintf(&g_buffer, "--mmap-pages %s", str);
}

bool isnum(char* str)
{
    for (unsigned int i = 0; i < strlen(str); i++) {
        if (!isdigit(str[i])) {
            return false;
        }
    }

    return true;
}

void proc_job(char* str)
{
    if (!isnum(str)) {
        printf("error: invalid job id: \"%s\"\n", str);
        exit(1);
    }

    g_job = str;
}

void proc_cpus(char* str)
{
    if (strcmp(str, "off") != 0) {
        asprintf(&g_cpus, "-C %s", str);
    }
    else {
        g_cpus = "";
    }
}

void proc_hostname()
{
    g_host = malloc(HOSTNAME * sizeof(char));
    memset(g_host, 0, HOSTNAME);
    if (gethostname(g_host, HOSTNAME)) {
        printf("error: gethostname(): (%d) %s\n",
                errno, strerror(errno));
    }
}

void proc_perfiles()
{
    //asprintf(&g_perfdat, "%s/job%s-perf-%s.dat",
    asprintf(&g_perfdat, "/tmp/job%s-perf-%s.dat",
            //g_dir, g_job, g_host);
            g_job, g_host);
    //asprintf(&g_perfout, "%s/job%s-perf-%s.out",
    asprintf(&g_perfout, "/tmp/job%s-perf-%s.out",
            //g_dir, g_job, g_host);
            g_job, g_host);
}

void proc_outfile()
{
    asprintf(&g_outfile, "%s/job%s-jt-%s.out", 
            g_dir, g_job, g_host);    
}

// -----------------------------------------------------------------------------

void regsig(int signum, __sighandler_t action)
{
    if (signal(signum, action) == SIG_ERR) {
        printf("[jt] error: signal(): (%d) %s\n",
                errno, strerror(errno));
        exit(1);
    }
}

void loop()
{
    while (true)
        sleep(SLEEP);
}

// -----------------------------------------------------------------------------

void atfprintf(char* fmt, ...)
{
    FILE* fp = fopen(g_outfile, "a");    
    va_list ap;
    va_start(ap, fmt);    
    vfprintf(fp, fmt, ap);
    va_end(ap);
    fclose(fp);
}

void handler_start(int signum __attribute__((unused)))
{
    atfprintf("[jt] handler_start\n");
    if (__atomic_test_and_set(&g_jt_started, __ATOMIC_ACQ_REL))
        return;
    g_jt_stopped = false;
    atfprintf("[jt] handler_start (after atomic)\n");

    char cmd[CMD] = {0};
    sprintf(cmd, "sudo /usr/bin/perf record %s %s %s -k CLOCK_MONOTONIC -a -o %s &> %s &",
            g_events, g_buffer, g_cpus, g_perfdat, g_perfout);    
    atfprintf("[jt] starting perf\n");
    atfprintf("%s\n", cmd);
    system(cmd);
}

void handler_stop(int signum __attribute__((unused)))
{
    atfprintf("[jt] handler_stop\n");
    if (__atomic_test_and_set(&g_jt_stopped, __ATOMIC_ACQ_REL))
        return;
    g_jt_started = false;
    atfprintf("[jt] handler_stop (after atomic)\n");

    char cmd[CMD] = {0};
    sprintf(cmd, "sudo killall -2 perf");
    atfprintf("[jt] stopping perf\n");
    atfprintf("%s\n", cmd);
    system(cmd);
}

// -----------------------------------------------------------------------------

int main(int argc, char** argv)
{
    printf("[jt] started\n");
    parse_args(argc, argv);    
    writepid(g_dir, g_job, g_host);
    regsig(SIGUSR1, handler_start);
    regsig(SIGUSR2, handler_stop);
    loop();
}

