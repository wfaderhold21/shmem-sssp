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

extern int64_t *distances;       
extern int64_t *active_vertices; 
extern int64_t *pi;              
extern int64_t base_vertex;      
extern int64_t orig_fair_offset; 
extern int64_t nr_edges;         
extern struct node *nodes;

extern char            filename[255];
extern size_t          shmem_memory_usage; 
extern size_t          local_memory_usage; 


/**
 * Simple function to read the test file 
 *
 * FIXME: This needs to be improved for large scales
 */ 
int64_t read_file(int64_t my_pe, 
                  int64_t nr_pes, 
                  int64_t * nr_nodes, 
                  int64_t * fair_offset, 
                  int mode) 
{
    FILE * fp;
    int64_t ret; 

    fp = fopen(filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "could not open %s\n", filename);
        return -1;
    }
    
    ret = fscanf(fp, "%ld", nr_nodes);
    if (ret == EOF) {
        fprintf(stderr, "ERROR: empty file?\n");
        return -1;
    }

    if (*nr_nodes < nr_pes) {
        *nr_nodes = nr_pes;
    } else {
        if ((*nr_nodes % nr_pes) > 0) {
            *nr_nodes = (*nr_nodes + (nr_pes - (*nr_nodes % nr_pes)));
            if (my_pe == 0) {
                printf("*** setting nr_nodes to %ld\n", *nr_nodes);
            }
        }
    }
    orig_fair_offset = *fair_offset = *nr_nodes / nr_pes;
    base_vertex = *fair_offset * my_pe;

    distances = (int64_t *) shmem_malloc(sizeof(int64_t)*(*fair_offset));
    if (distances == NULL) {
        fprintf(stderr, "Could not alloc distances\n");
        abort();
    }
    shmem_memory_usage += (sizeof(int64_t) * (*fair_offset));
    active_vertices = (int64_t *) shmem_malloc(sizeof(int64_t) * (*fair_offset));
    if (active_vertices == NULL) {
        fprintf(stderr, "Could not alloc active_vertices\n");
        abort();
    }
    shmem_memory_usage += (sizeof(int64_t) * (*fair_offset));
    
    pi = (int64_t *)shmem_malloc(sizeof(int64_t)*(*fair_offset));
    shmem_memory_usage += (sizeof(int64_t) * (*fair_offset));

    nodes = (struct node *) malloc(sizeof(struct node) * (*fair_offset));
    local_memory_usage += sizeof(struct node) * (*fair_offset);
    for (int64_t i = 0; i < *fair_offset; i++) {
        nodes[i].nr_edges = 0;
        nodes[i].nr_alloc_edges = 1;
        nodes[i].internal_edges_only = 1;
        nodes[i].edges = (struct edge *)malloc(sizeof(struct edge) * nodes[i].nr_alloc_edges);
        local_memory_usage += (sizeof(struct edge) * nodes[i].nr_alloc_edges);
        active_vertices[i] = 1;
    }

    for (;;) {
        int64_t a, b, d;
        ret = fscanf(fp, "%ld %ld %ld", &a, &b, &d);
        if(ret == EOF) {
            break;
        }
        nr_edges++;
        if (a == b) { /* self loop pruning */
            continue;
        }
        if (a >= base_vertex && a < (base_vertex + (*fair_offset))) {
            int index = a - base_vertex;
            int a_nr_edges = nodes[index].nr_edges;

            nodes[index].edges[a_nr_edges].dest = b;
            nodes[index].edges[a_nr_edges].weight = d;
            nodes[index].nr_edges++;

            if(nodes[index].nr_edges >= nodes[index].nr_alloc_edges) {
                local_memory_usage += (sizeof(struct edge) * nodes[index].nr_alloc_edges);  
                nodes[index].nr_alloc_edges *= 2;
                nodes[index].edges = realloc(nodes[index].edges, sizeof(struct edge) * nodes[index].nr_alloc_edges);
            }
        } else { // a does not belong to us
            if (b >= base_vertex && b < (base_vertex + (*fair_offset))) {
                //this is our vertex
                nodes[b-base_vertex].internal_edges_only = 0;
            }
        }
        if (mode == UNDIRECTED) {
            /* FIXME: i'm not sure this works */
            if (b >= base_vertex && b < (base_vertex + (*fair_offset))) {
            
                int index = b - base_vertex;
                int b_nr_edges = nodes[index].nr_edges;
                
                nodes[index].edges[b_nr_edges].dest = a;
                nodes[index].edges[b_nr_edges].weight = d;
                nodes[index].nr_edges++;

                if (nodes[index].nr_edges >= nodes[index].nr_alloc_edges) {
                    local_memory_usage += sizeof(struct edge) * nodes[index].nr_alloc_edges;
                    nodes[index].nr_alloc_edges *= 2;
                    nodes[index].edges = realloc(nodes[index].edges, sizeof(struct edge) * nodes[index].nr_alloc_edges);
                }
            }
        } 
    }
    fclose(fp);

    return 0;
}

