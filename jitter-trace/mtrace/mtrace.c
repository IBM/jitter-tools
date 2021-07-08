#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dlfcn.h>
#include <math.h>
#include <errno.h>

#include <mpi.h>

/* Override functions to trace */
int (*__Allgather)(const void *, int, MPI_Datatype, void *, int, MPI_Datatype, MPI_Comm);
int (*__Allreduce)(const void *, void *, int, MPI_Datatype, MPI_Op, MPI_Comm);
int (*__Alltoall)(const void *, int, MPI_Datatype, void *, int, MPI_Datatype, MPI_Comm);
int (*__Barrier)(MPI_Comm);
int (*__Bcast)(void *, int, MPI_Datatype, int, MPI_Comm);
int (*__Exscan)(const void *, void *, int, MPI_Datatype, MPI_Op, MPI_Comm);
int (*__Gather)(const void *, int, MPI_Datatype, void *, int, MPI_Datatype, int, MPI_Comm);
int (*__Reduce)(const void *, void *, int, MPI_Datatype, MPI_Op, int, MPI_Comm);
int (*__Reduce_scatter)(const void *, void *, const int *, MPI_Datatype, MPI_Op, MPI_Comm);
int (*__Scan)(const void *, void *, int, MPI_Datatype, MPI_Op, MPI_Comm);
int (*__Scatter)(const void *, int, MPI_Datatype, void *, int, MPI_Datatype, int, MPI_Comm);
int (*__Wait)(MPI_Request *, MPI_Status *);
int (*__Waitall)(int, MPI_Request *, MPI_Status *);
int (*__Waitany)(int, MPI_Request *, int *, MPI_Status *);

int (*__PAllgather)(const void *, int, MPI_Datatype, void *, int, MPI_Datatype, MPI_Comm);
int (*__PAllreduce)(const void *, void *, int, MPI_Datatype, MPI_Op, MPI_Comm);
int (*__PAlltoall)(const void *, int, MPI_Datatype, void *, int, MPI_Datatype, MPI_Comm);
int (*__PBarrier)(MPI_Comm);
int (*__PBcast)(void *, int, MPI_Datatype, int, MPI_Comm);
int (*__PExscan)(const void *, void *, int, MPI_Datatype, MPI_Op, MPI_Comm);
int (*__PGather)(const void *, int, MPI_Datatype, void *, int, MPI_Datatype, int, MPI_Comm);
int (*__PReduce)(const void *, void *, int, MPI_Datatype, MPI_Op, int, MPI_Comm);
int (*__PReduce_scatter)(const void *, void *, const int *, MPI_Datatype, MPI_Op, MPI_Comm);
int (*__PScan)(const void *, void *, int, MPI_Datatype, MPI_Op, MPI_Comm);
int (*__PScatter)(const void *, int, MPI_Datatype, void *, int, MPI_Datatype, int, MPI_Comm);
int (*__PWait)(MPI_Request *, MPI_Status *);
int (*__PWaitall)(int, MPI_Request *, MPI_Status *);
int (*__PWaitany)(int, MPI_Request *, int *, MPI_Status *);

/* Type of events to trace */
enum {
        __event_noise_constant = 0,
        __event_noise_propagates = 1,
        __event_other = 2,
        __event_MAX,
};

const char *__event_names[] = {
        "noise_constant",
        "noise_propagates",
        "noise_other",
};

typedef struct mevent mevent_t;
struct mevent {
        double timestamp;
        float duration;
        int id;
};

#define MTRACE_PAGE_MAX 8192
typedef struct mtrace_page mtrace_page_t;
struct mtrace_page {
        long last;
        mevent_t events[MTRACE_PAGE_MAX];
        mtrace_page_t *next;
};

struct mtrace {
        double start_time;
        mevent_t *last;
        mtrace_page_t *pages;
        mtrace_page_t *last_page;
} __mtrace;


static inline void mtrace_register_event(int event_id)
{
	mtrace_page_t *page = __mtrace.last_page;
	if (!page || page->last == MTRACE_PAGE_MAX) {
		page = malloc(sizeof(mtrace_page_t));
		page->last = 0;
		page->next = NULL;

		if (__mtrace.pages)
			__mtrace.last_page->next = page;
		else
			__mtrace.pages = page;

		__mtrace.last_page = page;
        }

	mevent_t *e = &page->events[page->last++];
	e->timestamp = MPI_Wtime();
	e->id = event_id;

	__mtrace.last = e;
}

static int mtrace_log(const char *filename, const char *str, ...)
{
	char fname[32];
	int jobid, rank;
	FILE *f;
        char *ptr;
        va_list args;
        int ret;

        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        ptr = getenv("LSB_JOBID");
        jobid = (ptr) ? atoi(ptr) : getpid();
	snprintf(fname, 32, "mtrace.%d.%d.log", jobid, rank);
	f = (filename) ? fopen(filename, "a") : fopen(fname, "a");
	if (!f)
		return errno;

        va_start(args, str);
	ret = vfprintf(f, str, args);
        va_end(args);
	fclose(f);
        return ret;
}

/*
 * Iterate over all the events in the trace
 * @e: the event cursor
 */
#define for_each_event(e)                                               \
        mtrace_page_t *__p;                                             \
        int __i;                                                        \
        for (__p = __mtrace.pages; __p; __p = __p->next)                \
                for (__i = 1, e = __p->events; __i <= __p->last;        \
                     e = __p->events + __i++)

/*
 * It's better to calculate the duration of an event after the tracing has
 * been completed, to reduce cpu overhead during the tracing itself
 */
static inline void calculate_events_duration()
{
        double old_timestamp = __mtrace.start_time;
        mevent_t *e;

        for_each_event(e) {
                e->duration = e->timestamp - old_timestamp;
                old_timestamp = e->timestamp;
        }
}

/*
 * The equivalent iteration time is the compute_phase time that is affected by
 * the same amount of jitter of the entire communication pattern registered.
 *
 * It's calculated with a weighted average and it's the interpolation of the
 * interval lengths over the curve 1/x^(-5/8)
 */
static void print_equivalent_iteration_time()
{
        int rank, np;
        int i;
        mevent_t *e;
        double num[__event_MAX];
        double weight[__event_MAX];
        unsigned long long tot[__event_MAX];
        char fname[16];

        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &np);
	snprintf(fname, 16, "eq.%d.log", np);
        memset(num, 0, __event_MAX * sizeof(double));
        memset(weight, 0, __event_MAX * sizeof(double));
        memset(tot, 0, __event_MAX * sizeof(unsigned long long));

        for_each_event(e) {
                weight[e->id] += pow(e->duration, -5/8);
                num[e->id] += e->duration * pow(e->duration, -5/8);
                tot[e->id]++;
        }
        
        for (i = 0; i < __event_MAX; i++) {
                if (rank) {
                        MPI_Reduce(&weight[i], MPI_IN_PLACE, 1, MPI_DOUBLE,
                                   MPI_SUM, 0, MPI_COMM_WORLD);
                        MPI_Reduce(&num[i], MPI_IN_PLACE, 1, MPI_DOUBLE,
                                   MPI_SUM, 0, MPI_COMM_WORLD);
                        MPI_Reduce(&tot[i], MPI_IN_PLACE, 1, MPI_DOUBLE,
                                   MPI_SUM, 0, MPI_COMM_WORLD);
                } else {
                        MPI_Reduce(MPI_IN_PLACE, &weight[i], 1, MPI_DOUBLE,
                                   MPI_SUM, 0, MPI_COMM_WORLD);
                        MPI_Reduce(MPI_IN_PLACE, &num[i], 1, MPI_DOUBLE,
                                   MPI_SUM, 0, MPI_COMM_WORLD);
                        MPI_Reduce(MPI_IN_PLACE, &tot[i], 1, MPI_DOUBLE,
                                   MPI_SUM, 0, MPI_COMM_WORLD);

                        mtrace_log(fname, "%d (%s):\t%.8lf %010llu\n",
                                   !!i, __event_names[i], num[i] / weight[i],
                                   tot[i] / np);
                }
        }
        if (!rank)
                mtrace_log(fname, "Total time spent by application: %lf\n",
                           MPI_Wtime() - __mtrace.start_time);
}

#ifdef full_log
static inline void mtrace_dump()
{
        mevent_t *e;
        for_each_event(e)
	        mtrace_log(NULL, "%.4lf %.9lf %s\n", e->timestamp, e->duration,
                           __event_names[e->id]);
}
#endif

static inline void __mtrace_initialize()
{
        __Allgather             = dlsym(RTLD_NEXT, "MPI_Allgather");
        __Allreduce             = dlsym(RTLD_NEXT, "MPI_Allreduce");
        __Alltoall              = dlsym(RTLD_NEXT, "MPI_Alltoall");
        __Barrier               = dlsym(RTLD_NEXT, "MPI_Barrier");
        __Bcast                 = dlsym(RTLD_NEXT, "MPI_Bcast");
        __Exscan                = dlsym(RTLD_NEXT, "MPI_Exscan");
        __Gather                = dlsym(RTLD_NEXT, "MPI_Gather");
        __Reduce                = dlsym(RTLD_NEXT, "MPI_Reduce");
        __Reduce_scatter        = dlsym(RTLD_NEXT, "MPI_Reduce_scatter");
        __Scatter               = dlsym(RTLD_NEXT, "MPI_Scatter");
        __Scan                  = dlsym(RTLD_NEXT, "MPI_Scan");
        __Wait                  = dlsym(RTLD_NEXT, "MPI_Wait");
        __Waitall               = dlsym(RTLD_NEXT, "MPI_Waitall");
        __Waitany               = dlsym(RTLD_NEXT, "MPI_Waitany");

        __PAllgather            = dlsym(RTLD_NEXT, "PMPI_Allgather");
        __PAllreduce            = dlsym(RTLD_NEXT, "PMPI_Allreduce");
        __PAlltoall             = dlsym(RTLD_NEXT, "PMPI_Alltoall");
        __PBarrier              = dlsym(RTLD_NEXT, "PMPI_Barrier");
        __PBcast                = dlsym(RTLD_NEXT, "PMPI_Bcast");
        __PExscan               = dlsym(RTLD_NEXT, "PMPI_Exscan");
        __PGather               = dlsym(RTLD_NEXT, "PMPI_Gather");
        __PReduce               = dlsym(RTLD_NEXT, "PMPI_Reduce");
        __PReduce_scatter       = dlsym(RTLD_NEXT, "PMPI_Reduce_scatter");
        __PScatter              = dlsym(RTLD_NEXT, "PMPI_Scatter");
        __PScan                 = dlsym(RTLD_NEXT, "PMPI_Scan");
        __PWait                 = dlsym(RTLD_NEXT, "PMPI_Wait");
        __PWaitall              = dlsym(RTLD_NEXT, "PMPI_Waitall");
        __PWaitany              = dlsym(RTLD_NEXT, "PMPI_Waitany");

        __mtrace.start_time = MPI_Wtime();
}

static inline void __mtrace_finalize()
{
        calculate_events_duration();
        print_equivalent_iteration_time();
#ifdef full_log
        mtrace_dump();
#endif
        mtrace_page_t *p, *n;
        for (p = __mtrace.pages; p; p = n) {
                n = p->next;
                free(p);
        }
	__mtrace.pages = __mtrace.last_page = NULL;
}

/* Override functions */
inline int PMPI_Init(int *argc, char ***argv)
{
        int (*__Init)(int *, char ***) = dlsym(RTLD_NEXT, "PMPI_Init");
        int ret =  __Init(argc, argv);
        __mtrace_initialize();
        return ret;
}
inline int MPI_Init(int *argc, char ***argv)
{
        int (*__Init)(int *, char ***) = dlsym(RTLD_NEXT, "MPI_Init");
        int ret =  __Init(argc, argv);
        __mtrace_initialize();
        return ret;
}
inline int MPI_Allgather(const void *send, int scount, MPI_Datatype st,
                         void *recv, int rcount, MPI_Datatype rt,
                         MPI_Comm comm)
{
        mtrace_register_event(__event_noise_propagates);
        return __Allgather(send, scount, st, recv, rcount, rt, comm);
}
inline int MPI_Allreduce(const void *send, void *recv, int count,
                         MPI_Datatype dt, MPI_Op op, MPI_Comm comm)
{
        mtrace_register_event(__event_noise_propagates);
        return __Allreduce(send, recv, count, dt, op, comm);
}
inline int MPI_Alltoall(const void *send, int scount, MPI_Datatype st,
                        void *recv, int rcount, MPI_Datatype rt, MPI_Comm comm)
{
        mtrace_register_event(__event_noise_propagates);
        return __Alltoall(send, scount, st, recv, rcount, rt, comm);
}
inline int MPI_Barrier(MPI_Comm comm)
{
        mtrace_register_event(__event_noise_propagates);
        return __Barrier(comm);
}
inline int MPI_Bcast(void *buf, int count, MPI_Datatype datatype, int root,
                     MPI_Comm comm)
{
        mtrace_register_event(__event_noise_propagates);
        return __Bcast(buf, count, datatype, root, comm);
}
inline int MPI_Exscan(const void *send, void *recv, int count,
                      MPI_Datatype dt, MPI_Op op, MPI_Comm comm)
{
        mtrace_register_event(__event_noise_propagates);
        return __Exscan(send, recv, count, dt, op, comm);
}
inline int MPI_Gather(const void *send, int scount, MPI_Datatype st,
                      void *recv, int rcount, MPI_Datatype rt, int root,
                      MPI_Comm comm)
{
        mtrace_register_event(__event_noise_constant);
        return __Gather(send, scount, st, recv, rcount, rt, root, comm);
}
inline int MPI_Reduce(const void *send, void *recv, int count,
                      MPI_Datatype dt, MPI_Op op, int root, MPI_Comm comm)
{
        mtrace_register_event(__event_noise_propagates);
        return __Reduce(send, recv, count, dt, op, root, comm);
}
inline int MPI_Reduce_scatter(const void *send, void *recv, const int count[],
                              MPI_Datatype dt, MPI_Op op, MPI_Comm comm)
{
        mtrace_register_event(__event_noise_propagates);
        return __Reduce_scatter(send, recv, count, dt, op, comm);
}
inline int MPI_Scatter(const void *send, int scount, MPI_Datatype st,
                       void *recv, int rcount, MPI_Datatype rt, int root,
                       MPI_Comm comm)
{
        mtrace_register_event(__event_noise_constant);
        return __Scatter(send, scount, st, recv, rcount, rt, root, comm);
}
inline int MPI_Scan(const void *send, void *recv, int count,
                    MPI_Datatype dt, MPI_Op op, MPI_Comm comm)
{
        mtrace_register_event(__event_noise_propagates);
        return __Scan(send, recv, count, dt, op, comm);
}
inline int MPI_Wait(MPI_Request *request, MPI_Status *status)
{
        mtrace_register_event(__event_other);
        return __Wait(request, status);
}
inline int MPI_Waitall(int count, MPI_Request array_of_requests[],
                       MPI_Status array_of_statuses[])
{
        mtrace_register_event(__event_other);
        return __Waitall(count, array_of_requests, array_of_statuses);
}
inline int MPI_Waitany(int count, MPI_Request array_of_requests[],
                       int *index, MPI_Status array_of_statuses[])
{
        mtrace_register_event(__event_other);
        return __Waitany(count, array_of_requests, index, array_of_statuses);
}
inline int PMPI_Allreduce(const void *send, void *recv, int count,
                          MPI_Datatype dt, MPI_Op op, MPI_Comm comm)
{
        mtrace_register_event(__event_noise_propagates);
        return __PAllreduce(send, recv, count, dt, op, comm);
}
inline int PMPI_Alltoall(const void *send, int scount, MPI_Datatype st,
                         void *recv, int rcount, MPI_Datatype rt,
                         MPI_Comm comm)
{
        mtrace_register_event(__event_noise_propagates);
        return __PAlltoall(send, scount, st, recv, rcount, rt, comm);
}
inline int PMPI_Barrier(MPI_Comm comm)
{
        mtrace_register_event(__event_noise_propagates);
        return __PBarrier(comm);
}
inline int PMPI_Bcast(void *buf, int count, MPI_Datatype datatype, int root,
                      MPI_Comm comm)
{
        mtrace_register_event(__event_noise_propagates);
        return __PBcast(buf, count, datatype, root, comm);
}
inline int PMPI_Exscan(const void *send, void *recv, int count,
                       MPI_Datatype dt, MPI_Op op, MPI_Comm comm)
{
        mtrace_register_event(__event_noise_propagates);
        return __PExscan(send, recv, count, dt, op, comm);
}
inline int PMPI_Gather(const void *send, int scount, MPI_Datatype st,
                       void *recv, int rcount, MPI_Datatype rt, int root,
                       MPI_Comm comm)
{
        mtrace_register_event(__event_noise_constant);
        return __PGather(send, scount, st, recv, rcount, rt, root, comm);
}
inline int PMPI_Reduce(const void *send, void *recv, int count,
                       MPI_Datatype dt, MPI_Op op, int root, MPI_Comm comm)
{
        mtrace_register_event(__event_noise_propagates);
        return __PReduce(send, recv, count, dt, op, root, comm);
}
inline int PMPI_Reduce_scatter(const void *send, void *recv, const int count[],
                               MPI_Datatype dt, MPI_Op op, MPI_Comm comm)
{
        mtrace_register_event(__event_noise_propagates);
        return __PReduce_scatter(send, recv, count, dt, op, comm);
}
inline int PMPI_Scatter(const void *send, int scount, MPI_Datatype st,
                        void *recv, int rcount, MPI_Datatype rt, int root,
                        MPI_Comm comm)
{
        mtrace_register_event(__event_noise_constant);
        return __PScatter(send, scount, st, recv, rcount, rt, root, comm);
}
inline int PMPI_Scan(const void *send, void *recv, int count,
                     MPI_Datatype dt, MPI_Op op, MPI_Comm comm)
{
        mtrace_register_event(__event_noise_propagates);
        return __PScan(send, recv, count, dt, op, comm);
}
inline int PMPI_Wait(MPI_Request *request, MPI_Status *status)
{
        mtrace_register_event(__event_other);
        return __PWait(request, status);
}
inline int PMPI_Waitall(int count, MPI_Request array_of_requests[],
                        MPI_Status array_of_statuses[])
{
        mtrace_register_event(__event_other);
        return __PWaitall(count, array_of_requests, array_of_statuses);
}
inline int PMPI_Waitany(int count, MPI_Request array_of_requests[],
                        int *index, MPI_Status array_of_statuses[])
{
        mtrace_register_event(__event_other);
        return __PWaitany(count, array_of_requests, index, array_of_statuses);
}

inline int PMPI_Finalize()
{
        int (*finalize)(void) = dlsym(RTLD_NEXT, "PMPI_Finalize");
        __mtrace_finalize();
        return finalize();
} 
inline int MPI_Finalize()
{
        int (*finalize)(void) = dlsym(RTLD_NEXT, "MPI_Finalize");
        __mtrace_finalize();
        return finalize();
}
