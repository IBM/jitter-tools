#include <errno.h>     // errno
#include <signal.h>    // kill, SIGUSR1, SIGUSR2
#include <stdbool.h>   // bool, true, false
#include <stdio.h>     // fprintf, stderr
#include <stdlib.h>    // getenv, exit
#include <string.h>    // strerror
#include <sys/types.h> // pid_t
#include <unistd.h>    // gethostname

#include "pidfile.c"   // readpid

#define HOSTNAME 1024  // length of hostname

bool  g_jt_on = false;
pid_t g_jt_pid;

void jt_start()
{
    //fprintf(stderr, "[libjt] init\n");
    const char* jt_on_s;
    if (!(jt_on_s = getenv("JT_ON"))) {
        //fprintf(stderr, "[libjt] warning: missing JT_ON\n");
        return;
    }
    g_jt_on = atoi(jt_on_s);
    if (!g_jt_on) {
        return;
    }

    const char* jt_dir_s;
    if (!(jt_dir_s = getenv("JT_DIR"))) {
        //fprintf(stderr, "[libjt] warning: missing JT_DIR\n");
        g_jt_on = false;
        return;
    }

    const char* jt_job_s;
    if (!(jt_job_s = getenv("JT_JOB"))) {
        //fprintf(stderr, "[libjt] warning: missing JT_JOB\n");
        g_jt_on = false;
        return;
    }

    char host[HOSTNAME] = {0};
    if (gethostname(host, HOSTNAME)) {
        //fprintf(stderr, "[libjt] warning: gethostname() (%d) %s\n",
        //        errno, strerror(errno));
        g_jt_on = false;
        return;
    }

    g_jt_pid = readpid(jt_dir_s, jt_job_s, host);

    //
    //  OLD START
    //

    if (!g_jt_on || !g_jt_pid) {
        return;
    }

    //fprintf(stderr, "[libjt] start\n");
    if (kill(g_jt_pid, SIGUSR1)) {
        //fprintf(stderr, "[libjt] warning: start(): (%d) %s\n",
        //        errno, strerror(errno));
    }

    // wait for perf to actually start
    sleep(10);
}

void jt_stop()
{
    if (!g_jt_on || !g_jt_pid) {
        return;
    }

    //fprintf(stderr, "[libjt] stop\n");
    if (kill(g_jt_pid, SIGUSR2)) {
        //fprintf(stderr, "[libjt] warning: stop(): (%d) %s\n",
        //        errno, strerror(errno));
    }

    // no need for sleep here; calling program/script should determine it
}

