/**
 * An implementation of the Bellman-Ford algorithm using SHMEM for 
 * parallelization
 *
 * optimizations used in algorithm:
 *  1. early end detection
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include <sched.h>

#include <shmem.h>

#include "common.h"

int64_t *distances;       /* node distances/predecessors from source */
int64_t *active_vertices; /* active vertices for a PE */
int64_t *pi;              /* predecessor list */
int64_t *active_vertex;   /* active vertex flag per PE */

int64_t nr_nodes         = 0; /* nr of vertices */
int64_t base_vertex      = 0; /* the base vertex for this PE */
int64_t orig_fair_offset = 0; /* the uniform partition size */
int64_t was_extended     = 0; /* if the partition size was extended */
int64_t nr_edges         = 0; /* the number of edges */
int64_t one_val          = 1;
int64_t neg_one_val      = -1;

struct node *nodes;
int          nr_alloc_edges;

char            filename[255];
size_t          shmem_memory_usage; // in bytes
size_t          local_memory_usage; // in bytes
static long     pWrk[SHMEM_REDUCE_MIN_WRKDATA_SIZE];
static long     pSync[SHMEM_REDUCE_SYNC_SIZE];
static int      pWrk2[SHMEM_REDUCE_MIN_WRKDATA_SIZE];

int64_t read_file(int64_t my_pe, int64_t nr_pes, int64_t *nr_nodes,
                  int64_t *fair_offset, int mode);

void relax(int64_t i, int64_t j, int64_t weight, int64_t nr_pes, int64_t my_pe)
{
    int64_t distance_j = 0, orig_j;
    int64_t loc        = j % (orig_fair_offset);
    int64_t pe         = j / (orig_fair_offset);
    int64_t pi_j       = 0, orig_pi;
    int64_t compared_distance;

    shmem_get64(&distance_j, &distances[loc], 1, pe);
    orig_j = distance_j;
    shmem_get64(&orig_pi, &pi[loc], 1, pe);
    compared_distance = distances[i] + weight;

    /* relax function */
    if (distance_j > compared_distance) {
        distance_j = compared_distance; 

        /* Purely local operation */
        if (pe == my_pe && nodes[j - base_vertex].internal_edges_only) {
            distances[loc]       = distance_j;
            pi[loc]              = base_vertex + i;
            *active_vertex       = 1;
            active_vertices[loc] = 1;
        }
        else {
            /* only try up to 10 times to ensure no deadlocks */
            for (int d_tries = 0; d_tries < 10; d_tries++) {
                int64_t old_j;
                old_j =
                    shmem_long_cswap(&distances[loc], orig_j, distance_j, pe);

                if (old_j == orig_j) {
                    pi_j = base_vertex + i;
                    for (int tries = 0; tries < 100; tries++) {
                        int64_t old_p;

                        old_p = shmem_long_cswap(&pi[loc], orig_pi, pi_j, pe);
                        if (old_p == orig_pi) {
                            break;
                        } else {
                            // FIXME: could have an
                            // issue here due to race
                            // condition between
                            // atomics
                            orig_pi = old_p;
                        }
                    }
                    if (my_pe == pe) {
                        active_vertices[loc] = 1;
                        *active_vertex       = 1;
                    }
                    else {
                        shmem_long_put(&active_vertices[loc], &one_val, 1, pe);
                        shmem_long_put(active_vertex, &one_val, 1, pe);
                    }
                    // TODO: replace with shmem_pe_quiet()
                    shmem_quiet();
                    break;
                }
                else {
                    /* someone else set this before us */
                    if (old_j < distance_j) {
                        break;
                    }
                    orig_j = old_j;
                }
            }
        }
    }
}

int bellman_ford_synchronous(int64_t nr_pes, int64_t my_pe, int64_t fair_offset,
                             int64_t nr_nodes, int *converged,
                             uint64_t *nr_traversed_edges)
{
    for (int k = 1; k < nr_nodes - 1; k++) {
        int        was_active    = 0;
        static int all_converged = 0;

        if (my_pe == 0) {
            converged[0] = 0;
        }

        if (*active_vertex == 1) {
            *active_vertex = 0;
            was_active     = 1;

            for (int64_t i = 0; i < fair_offset; i++) {
                if (active_vertices[i] == 1) {
                    active_vertices[i] = 0;
                    for (int64_t j = 0; j < nodes[i].nr_edges; j++) {
                        relax(i, nodes[i].edges[j].dest,
                              nodes[i].edges[j].weight, nr_pes, my_pe);
                        *nr_traversed_edges = *nr_traversed_edges + 1;
                    }
                }
            }
        }
        if (was_active == 0) {
            converged[0] = 1;
        }

        shmem_int_sum_to_all(&all_converged, converged, 1, 0, 0, nr_pes, pWrk2,
                             pSync);
        if (all_converged >= nr_pes) {
            break;
        }
    }
    return 0;
}

int main(int argc, char **argv)
{
    int64_t   nr_nodes;
    int64_t   my_pe, nr_pes;
    int64_t   fair_offset = 0;
    int       root     = 0;
    int       root_loc = 0;
    int       root_pe  = 0;
    int       mode     = 0;
    double    start, end;
    int *     converged;
    uint64_t *nr_traversed_edges;
    uint64_t *total_edges;
    int       nr_iterations = 0;
    double *  avg_times;
    double *  avg_teps;

    if (argc < 5) {
        fprintf(stderr, "usage:\n");
        fprintf(stderr,
                "\t%s <root> <directed | undirected> <graph file> <# of "
                "iterations>\n",
                argv[0]);
        return -1;
    }

    root = atoi(argv[1]);
    if (strcmp(argv[2], "undirected") == 0) {
        mode = UNDIRECTED;
    }
    strcpy(filename, argv[3]);
    nr_iterations      = atoi(argv[4]);
    local_memory_usage = shmem_memory_usage = 0;

    shmem_init();
    my_pe  = shmem_my_pe();
    nr_pes = shmem_n_pes();

    start = shmem_wtime();
    read_file(my_pe, nr_pes, &nr_nodes, &fair_offset, mode);
    end = shmem_wtime();
    if (my_pe == 0) {
        printf("Graph load time: %g s\n", end - start);
    }

    avg_times = (double *)malloc(sizeof(double) * nr_iterations);
    avg_teps  = (double *)malloc(sizeof(double) * nr_iterations);
    /* ferrol: i am not counting local timing allocations as part of timing
       measurements */
    active_vertex = (int64_t *)shmem_malloc(sizeof(int64_t));
    shmem_memory_usage += sizeof(int64_t);

    converged = (int *)shmem_malloc(sizeof(int));
    shmem_memory_usage += sizeof(int);
    converged[0] = 0;

    nr_traversed_edges = (uint64_t *)shmem_malloc(sizeof(uint64_t));
    total_edges        = (uint64_t *)shmem_malloc(sizeof(uint64_t));

    root_pe  = root / orig_fair_offset;
    root_loc = root % orig_fair_offset;

    if (root_pe >= nr_pes) {
        if (was_extended) {
            root_pe  = nr_pes - 1;
            root_loc = orig_fair_offset + was_extended - 1;
        }
    }

    /* finished init */
    for (int ii = 0; ii < nr_iterations; ii++) {
        /** BEGINNING OF AN ITERATION **/
        *nr_traversed_edges = 0;
        for (int64_t i = 0; i < fair_offset; i++) {
            distances[i]       = UNDEFINED;
            pi[i]              = -1;
            active_vertices[i] = 1;
        }
        converged[0] = 0;

        *active_vertex = 1;

        if (my_pe == root_pe) {
            distances[root_loc] = 0;
        }

        shmem_barrier_all();
        start = shmem_wtime();

        bellman_ford_synchronous(nr_pes, my_pe, fair_offset, nr_nodes,
                                 converged, nr_traversed_edges);
        end = shmem_wtime();

#if BF_DEBUG == 1
        for (int i = 0; i < nr_pes; i++) {
            if (my_pe == i) {
                for (int64_t j = 0; j < fair_offset; j++) {
                    int64_t distance = distances[j];
                    int64_t pi_j     = pi[j];

                    dprintf("distance[%ld] = %ld\n", j + base_vertex, distance);
                    dprintf("pi[%ld] = %ld\n", j + base_vertex, pi_j);
                }
            }
            shmem_barrier_all();
        }
#endif

        double total_time = end - start;
        double teps;
        shmem_long_sum_to_all((long *) total_edges, (long *) nr_traversed_edges, 1, 0, 0, nr_pes,
                              pWrk, pSync);
        teps = *total_edges / total_time;

        avg_times[ii] = total_time;
        avg_teps[ii]  = teps;

        shmem_barrier_all();
        if (my_pe == 0) {
            printf("iteration %d\n", ii);
            printf("\tGraph SSSP time: %g s\n", total_time);
            printf("\tGraph TEPS: %g\n", teps);
            printf("\tTraversed edges: %ld\n", *total_edges);
        }
    }

    if (my_pe == 0) {
        double a_time = 0.0, a_teps = 0.0, b_teps = 0.0, b_time = 0.0;
        printf("memory usage:\n");
        printf("\tsymmetric heap usage: %lu bytes\n", shmem_memory_usage);
        printf("\tlocal memory usage: %lu bytes\n", local_memory_usage);

        for (int i = 0; i < nr_iterations; i++) {
            a_time += 1.0 / avg_times[i];
            a_teps += 1.0 / avg_teps[i];
            b_teps += avg_teps[i];
            b_time += avg_times[i];
        }
        a_time = nr_iterations / a_time;
        a_teps = nr_iterations / a_teps;
        b_teps = b_teps / nr_iterations;
        b_time = b_time / nr_iterations;
        printf("harmonic mean (time): %g s\n", a_time);
        printf("harmonic mean (TEPS): %g\n", a_teps);
        printf("Time mean: %g s\n", b_time);
        printf("TEPS mean: %g\n", b_teps);
    }
    shmem_finalize();

    return 0;
}
