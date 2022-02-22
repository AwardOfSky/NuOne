#ifndef HASHMAP_H
#define HASHMAP_H

#include <stdint.h>
#include <math.h>

#define MAX_TREE_NODES (1 << 30)
#define ENABLE_CACHE 0
#define CHILDREN_ARR_RATIO 2.0
#define SCALAR_INIT_SIZE 1024
#define CHILD_INIT_SIZE 1024
#define MAX_SCALAR_SIZE 1073741824
#define MAX_CHILD_SIZE 1073741824
#define N_CHILDREN_ARR 64

//temp cache stuff
#define OA_NODES 4
#define RADIX_SORT_THRESH 64
#define UPDATE_CACHE_N_GEN 5
#define FACTOR_TO_GROW_CACHE 4
#define ENABLE_CACHE 0
#define MAX_CACHE_SIZE 16384
#define MIN_CACHE_SIZE 1
#define MAX_CACHE_CANDIDATES 65536
#define MIN_CACHE_REQ 0
#define IMPORTANCE(FREQ, NODES) (FREQ * NODES)
#define INIT_CACHE_SIZE 1024
#define INIT_CAND_SIZE 8192


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
    //cache_node *candidates;
    uint32_t *candidates;
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


HashTable *create_hashtable(uint32_t n);
void free_scalar_ptrs(HashTable *t);
void free_child_ptrs(HashTable *t);
void free_hashtable(HashTable *t);
void print_dag_table(b_node *table, int n);
void delete_dag_table(b_node *table, uint32_t tab_size);
uint32_t hash_dag_node(uint32_t primitive, int arity, c_node children);
dag_node **alloc_child_ptrs(HashTable *t, int slots);
float *alloc_scalar_ptrs(HashTable *t, int slots);
uint32_t hash_dag_node(uint32_t primitive, int arity, c_node children);
dag_node *add_dag_node(HashTable *t, uint32_t primitive, c_node children);
dag_node *add_dag_node_index(HashTable *t, uint32_t primitive, c_node children, int ind);

#endif