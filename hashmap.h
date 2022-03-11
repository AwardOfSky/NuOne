#ifndef HASHMAP_H
#define HASHMAP_H

#include <stdint.h>
#include <math.h>

#define MAX_TREE_NODES (1 << 30)
#define CHILDREN_ARR_RATIO 2.0
#define SCALAR_INIT_SIZE 1024
#define CHILD_INIT_SIZE 1024
#define MAX_SCALAR_SIZE 1073741824
#define MAX_CHILD_SIZE 1073741824
#define N_CHILDREN_ARR 64
#define INIT_HASHMAP_SIZE 65535


#define ALLOC_PTRS(T, PTRS, IND, ARR, SIZE, TYPE, MSG) do {                             \
    if (T->IND + slots > T->SIZE) {                                                     \
        T->ARR++;                                                                       \
        T->IND = slots;                                                                 \
                                                                                        \
        if (T->ARR < N_CHILDREN_ARR) {                                                  \
            T->SIZE = min(T->SIZE * CHILDREN_ARR_RATIO, MAX_CHILD_SIZE);                \
            T->PTRS[T->ARR] = (TYPE*)malloc(t->cur_child_arr_size * sizeof(TYPE));      \
        } else {                                                                        \
            printf(PREERR"Ran out of arrays to allocate pointers to childs ("MSG")\n"); \
            T->ARR = 0;                                                                 \
        }                                                                               \
        return &T->PTRS[T->ARR][0];                                                     \
    }                                                                                   \
                                                                                        \
    TYPE*res_ptr = &T->PTRS[t->ARR][t->IND];                                            \
    T->IND += slots;                                                                    \
    return res_ptr;                                                                     \
} while(0)


#define FREE_PTRS(T, PTRS, ARR) do {    \
    for(int i = 0; i <= T->ARR; ++i) {  \
        free(T->PTRS[i]);               \
    }                                   \
} while(0)


// 32 bytes
typedef struct b_node {
    struct b_node *next;
    dag_node node;
} b_node;

// 64 bytes
typedef struct HashTable {
    b_node *table;
    uint32_t size;
    dag_node ***child_ptrs;
    float **scalar_ptrs;
    struct cache_node *candidates;
    uint32_t cand_index;
    uint32_t cand_size;
    uint32_t n_nodes;
    int scalar_arr;
    int scalar_ind;
    int cur_scalar_arr_size;
    int child_arr;
    int child_ind;
    int cur_child_arr_size;
} HashTable;

#include "cache.h"

HashTable *create_hashtable(uint32_t n);
void print_table_stats(b_node *table, int n, int do_extra_stats);
void free_scalar_ptrs(HashTable *t);
void free_child_ptrs(HashTable *t);
void free_hashtable(HashTable *t);
void print_dag_table(b_node *table, int n);
void delete_dag_table(b_node *table, uint32_t tab_size);
uint32_t hash_dag_node(uint32_t primitive, int arity, c_node children);
dag_node **alloc_child_ptrs(HashTable *t, int slots);
float *alloc_scalar_ptrs(HashTable *t, int slots);
uint32_t hash_dag_node(uint32_t primitive, int arity, c_node children);
void print_candidate_list(HashTable *t, int print_cache_cands);
dag_node *add_dag_node(HashTable *t, uint32_t primitive, c_node children);
dag_node *add_dag_node_index(HashTable *t, uint32_t primitive, c_node children, int ind);


#endif