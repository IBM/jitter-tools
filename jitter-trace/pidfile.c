#define LEN 1024

void writepid(const char* dir, const char* job, const char* hostname)
{
    char pid_fname[LEN] = {0};
    sprintf(pid_fname, "%s/job%s-jt-%s.pid", dir, job, hostname); 

    FILE* pid_f = fopen(pid_fname, "w");
    fprintf(pid_f, "%d\n", getpid());
    fclose(pid_f);
}

pid_t readpid(const char* dir, const char* job, const char* hostname)
{
    char pid_fname[LEN] = {0};
    sprintf(pid_fname, "%s/job%s-jt-%s.pid", dir, job, hostname);

    pid_t pid = 0;
    FILE* pid_f = fopen(pid_fname, "r");
    if (pid_f) {
        fscanf(pid_f, "%d", &pid);
        fclose(pid_f);
    }
    return pid;
}

