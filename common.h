#ifndef _COMMON_H

#define DIRECTED            0
#define UNDIRECTED          1
#define UNDEFINED           INT_MAX
#define UNDEFINED_EDGE      INT_MIN

#ifdef BF_DEBUG
#define dprintf(fmt, ...) printf(fmt, __VA_ARGS__)
#else
#define dprintf(...) \
        do {         \
        } while (0);
#endif

static inline double shmem_wtime(void) {
#ifdef USE_CLOCK
        clock_t uptime = clock();
        return (((double)uptime) / ((double)CLOCKS_PER_SEC));
#else
#include <sys/time.h>
        double wtime = 0;
        struct timeval tv;
        gettimeofday(&tv, NULL);
        wtime = tv.tv_sec;
        wtime += (double)tv.tv_usec / 1000000.0;
        return wtime;
#endif
}

struct node {
        struct edge* edges;
        int64_t nr_edges;
        int64_t nr_alloc_edges;
        int64_t internal_edges_only;
};

struct edge {
        int64_t dest;
        int64_t weight;
};

#endif
