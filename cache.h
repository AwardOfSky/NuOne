#ifndef CACHE_H
#define CACHE_H


#define ENABLE_CACHE 1
#define TENSOR_CACHE 1

// NOTE: none of the factors should be lower than 1
// a higher shrink factor means that is harder (or less likely) to shrink the cache,
// which in turn promotes bigger cache sizes
// a higher growth factor means that is easier (or more likely) to grow the cache,
// which also promotes bigger cache sizes
// @deprecated ^ 

#define OA_NODES 4
#define RADIX_SORT_THRESH 128
#define FACTOR_TO_GROW_CACHE 2
#define FACTOR_TO_SHRINK_CACHE 3
#define FACTOR_TO_GROW_CANDLIST 2

// Both the max cache memory (in bytes) and the max cache size can be set to -1
// to ignore the setting. Setting both macros to -1 is not advisable
// as it will allow the cache to grow without boundaries!
#define MAX_CACHE_MEMORY (1024 * 1024 * 1024)
#define MAX_CACHE_SIZE 16384

#define MIN_CACHE_SIZE 2
#define MAX_CACHE_CANDIDATES 65536
#define MIN_CACHE_REQ 0
#define IMPORTANCE(FREQ, NODES) (FREQ * NODES)
#define INIT_CACHE_SIZE 1024
#define INIT_CAND_SIZE 4
#define GROW_CACHE_N_GEN 1
#define SHRINK_CACHE_N_GEN 5
#define SUB_CAND_IMP_DEP_LIM 2



// 16 bytes
typedef struct cache_node {
    dag_node *dag;
    int index;
    int importance;
} cache_node;

// 32 bytes
typedef struct Cache {
    fit_t *data;
    fit_t *vars;
    char *valid;
    uint32_t occupied;
    uint32_t size;
    uint32_t accesses;
    uint32_t validated;
} Cache;

Cache cache;

int get_cand_imp_to_sub(HashTable *t, dag_node *dag, int32_t dep_lim);
void subtract_cand_imp(HashTable *t, dag_node *dag, int dep_lim, int *imp, int *freq);
void handle_candidates(dag_node *dag, HashTable *t, int ind);
void init_cache(size_t n);
void free_cache();

DECLARE_REALLOC(cache_node);
DECLARE_REALLOC(fit_t);
// float *realloc_fit_t(float *data, const uint32_t new_size, uint32_t *sptr);

//fit_t *realloc_fit_t(fit_t *data, const uint32_t new_size, uint32_t *sptr);

void realloc_cache(uint32_t new_size);
void print_cache(uint32_t *format, int print_invalid);
int cmp_candidates(const void *a, const void *b);
void build_cache(HashTable *t, int gen);

#endif
