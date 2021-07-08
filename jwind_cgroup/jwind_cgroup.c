/* Copyright IBM corp.
 * author: Bryan Rosenburg, Fabio Checconi
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <errno.h>
#include <fcntl.h>
#include <signal.h>
#include <time.h>
#include <sched.h>
#include <unistd.h>
#include <string.h>

#define JWIND_PRIO	51
#define HIGH_PRIO	50
#define LOW_PRIO	0

#define JWIND_POLICY	SCHED_FIFO
#define HIGH_POLICY	SCHED_FIFO
#define LOW_POLICY	SCHED_OTHER

#define PROCESS_MAX 1024

char *program;
int running = 1;

typedef struct {
    int pid;
    int active;
} process_t;
process_t process[PROCESS_MAX];
int process_count = 0;

void update_processes(char *filename)
{
    static uint64_t last_sum = 0;

    int file = open(filename, O_RDONLY, 0);
    if (file < 0) {
	fprintf(stderr, "%s: open(\"%s\") failed (%s).\n",
		program, filename, strerror(errno));
	exit(-1);
    }

    char buffer[PROCESS_MAX*8];
    int cnt = read(file, buffer, sizeof(buffer));
    if (cnt < 0) {
	fprintf(stderr, "%s: read(<%s>) failed (%s).\n",
		program, filename, strerror(errno));
	exit(-1);
    }
    if (cnt == sizeof(buffer)) {
	fprintf(stderr, "%s: tasks file \"%s\" too large.\n",
		program, filename);
	exit(-1);
    }

    close(file);

    // zero-fill last word and then sum
    int i;
    for (i = cnt; i < (((cnt + 7) / 8) * 8); i++) buffer[i] = '\0';
    uint64_t sum = 0;
    for (i = 0; i < ((cnt + 7) / 8); i++) {
	sum += ((uint64_t *) buffer)[i];
    }

    // assume no changes if the sum hasn't changed
    if (sum == last_sum) return;

    last_sum = sum;

    process_count = 0;
    int my_pid = getpid();

    char *p = buffer;
    while (p < (buffer + cnt)) {
	if (process_count >= PROCESS_MAX) {
	    fprintf(stderr, "%s: too many processes in tasks file \"%s\".\n",
		    program, filename);
	    exit(-1);
	}

	int pid = strtol(p, &p, 10);
	if (*p++ != '\n') {
	    fprintf(stderr, "%s: error parsing tasks file \"%s\".\n",
		    program, filename);
	    exit(-1);
	}

	if (pid != my_pid) {
	    process[process_count].pid = pid;
	    process[process_count].active = 1;
	    process_count++;
	}
    }
}

void precise_sleep(timer_t *timerp, struct itimerspec *timep,
		   uint64_t duration_us, sigset_t *sigsetp)
{
    timep->it_value.tv_nsec += (duration_us * 1000);
    if (timep->it_value.tv_nsec >= 1000000000) {
	timep->it_value.tv_sec += (timep->it_value.tv_nsec / 1000000000);
	timep->it_value.tv_nsec = (timep->it_value.tv_nsec % 1000000000);
    }
    timer_settime(timerp, TIMER_ABSTIME, timep, NULL);

    int sig;
    int rc = sigwait(sigsetp, &sig);
    if (rc != 0) {
	fprintf(stderr, "%s: sigwait failed (%s).\n",
		program, strerror(errno));
	exit(-1);
    }
}

void set_process_priorities(int policy, int prio)
{
    struct sched_param sched;
    sched.sched_priority = prio;
    int i;
    for (i = 0; i < process_count; i++) {
	if (process[i].active) {
	    int rc = sched_setscheduler(process[i].pid, policy, &sched);
	    if (rc != 0) {
		if (errno == ESRCH) {
		    // process has gone away; just mark it inactive
		    process[i].active = 0;
		} else {
		    fprintf(stderr, "%s: sched_setscheduler failed (%s).\n",
			    program, strerror(errno));
		    exit(-1);
		}
	    }
	}
    }
}

void stop_running(int sig)
{
    running = 0;
}

int set_sysfs_value(char *filename, int value)
{
    char value_string[16];
    int rc;
    int file = open(filename, O_RDWR, 0);
    if (file < 0) {
	fprintf(stderr, "%s: open(\"%s\") failed (%s)\n",
		program, filename, strerror(errno));
	exit(-1);
    }
    rc = pread(file, value_string, sizeof(value_string), 0);
    if (rc < 0) {
	fprintf(stderr, "%s: pread(<%s>) failed (%s)\n",
		program, filename, strerror(errno));
	exit(-1);
    }
    if (rc == sizeof(value_string)) {
	fprintf(stderr, "%s: pread(<%s>) return value too large\n",
		program, filename);
	exit(-1);
    }
    int old_value = strtol(value_string, NULL, 0);
    snprintf(value_string, sizeof(value_string) - 1, "%d", value);
    value_string[sizeof(value_string) - 1] = '\0'; // just in case
    rc = pwrite(file, value_string, strlen(value_string), 0);
    if (rc != strlen(value_string)) {
	fprintf(stderr, "%s: pwrite(<%s>) failed (%s)\n",
		program, filename, strerror(errno));
	exit(-1);
    }
    return old_value;
}

void usage(FILE *out)
{
    fprintf(out,
	    "Usage:\n"
	    "    %s -c <cgroup_name> -p <period usec> -w <window usec>\n",
	    program);
}

int main(int argc, char *argv[])
{
    int rc;

    program = argv[0];

    char *cgroup = NULL;
    uint64_t period_us = 0;
    uint64_t window_us = 0;

    int opt;
    while ((opt = getopt(argc, argv, "c:p:w:h")) != -1) {
	switch (opt) {
	case 'c':
	    cgroup = optarg;
	    break;
	case 'p':
	    period_us = strtoul(optarg, NULL, 0);
	    break;
	case 'w':
	    window_us = strtoul(optarg, NULL, 0);
	    break;
	case 'h':
	default:
	    usage(stdout);
	    exit(0);
	}
    }

    if ((cgroup == NULL) ||
	(period_us == 0) || (window_us == 0) ||
	(period_us < window_us))
    {
	usage(stderr);
	fprintf(stderr,
		"\n"
		"    cgroup_name (-c) must be provided, period (-p) and\n"
		"    window (-w) must be specified, and window must not\n"
		"    exceed period.\n");
	exit(-1);
    }

    char tasks_file_name[128];
    snprintf(tasks_file_name, 127,
	     "/sys/fs/cgroup/cpuset/%s/tasks", cgroup);
    tasks_file_name[127] = '\0'; // just in case

    char rt_runtime_file_name[128];
    snprintf(rt_runtime_file_name, 127,
	     "/sys/fs/cgroup/cpu/%s/cpu.rt_runtime_us", cgroup);
    rt_runtime_file_name[127] = '\0'; // just in case

    // Remove any real-time scheduling constraints
    int old_global_runtime =
	set_sysfs_value("/proc/sys/kernel/sched_rt_runtime_us", -1);
    int old_parent_runtime =
	set_sysfs_value("/sys/fs/cgroup/cpu/cpu.rt_runtime_us", -1);
    int old_cgroup_runtime =
	set_sysfs_value(rt_runtime_file_name, -1);

    // Block SIGALRM
    sigset_t sigset_alrm;
    sigemptyset(&sigset_alrm);
    sigaddset(&sigset_alrm, SIGALRM);
    sigprocmask(SIG_BLOCK, &sigset_alrm, NULL);

    // Create timer
    timer_t timer;
    timer_create(CLOCK_REALTIME, NULL, &timer);

    // Set time as if we were at the start of the previous noise window
    struct timespec now;
    clock_gettime(CLOCK_REALTIME, &now);
    uint64_t now_us = (now.tv_sec * 1000000ul) + (now.tv_nsec / 1000ul);
    uint64_t noise_window_us = ((now_us / period_us) * period_us);
    struct itimerspec time;
    time.it_value.tv_sec = noise_window_us / 1000000;
    time.it_value.tv_nsec = (noise_window_us % 1000000) * 1000;
    time.it_interval.tv_sec = 0;
    time.it_interval.tv_nsec = 0;

    // Catch termination signals
    signal(SIGINT, stop_running);
    signal(SIGQUIT, stop_running);

    // Set our own priority
    struct sched_param sched;
    sched.sched_priority = JWIND_PRIO;
    rc = sched_setscheduler(0, JWIND_POLICY, &sched);
    if (rc != 0) {
	fprintf(stderr, "%s: sched_setscheduler failed (%s).\n",
		program, strerror(errno));
	exit(-1);
    }

    int trace = open("/sys/kernel/debug/tracing/trace_marker", O_WRONLY, 0);
    if (trace < 0) {
	fprintf(stderr, "%s: open(\".../trace_marker\") failed (%s).\n",
		program, strerror(errno));
	exit(-1);
    }

    while (running) {
	(void) write(trace, "noise window\n", 13);
	update_processes(tasks_file_name);

	// Sleep to end of noise window
	precise_sleep(&timer, &time, window_us, &sigset_alrm);

	set_process_priorities(HIGH_POLICY, HIGH_PRIO);

	// Sleep to start of noise window
	precise_sleep(&timer, &time, period_us - window_us, &sigset_alrm);

	set_process_priorities(LOW_POLICY, LOW_PRIO);
    }

    // Restore real-time scheduling constraints
    (void) set_sysfs_value(rt_runtime_file_name, old_cgroup_runtime);
    (void) set_sysfs_value("/sys/fs/cgroup/cpu/cpu.rt_runtime_us",
			   old_parent_runtime);
    (void) set_sysfs_value("/proc/sys/kernel/sched_rt_runtime_us",
			   old_global_runtime);

    return 0;
}
