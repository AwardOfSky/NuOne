#ifndef GP_H
#define GP_H

#include <stdint.h>
#include <string.h>
#include "utils.h"

#define DEBUG 0

#define TSET 2
#define PSET_LEN (PSET_END - PSET_START + 1)
#define TSET_LEN (TSET_END - PSET_START + 1)
#define FSET_LEN (PSET_END - TSET_END)
#define MIN_SCALAR -5.0
#define MAX_SCALAR 5.0
#define DELTA_SCALAR (MAX_SCALAR - MIN_SCALAR)
#define VAR_START 2
#define VAR_END  3
#define MAX_PRIM_NAME_SIZE 256
#define DOMAIN_DELTA(E,I) (E->MAX_DOMAIN[I] - E->MIN_DOMAIN[I])
#define MAX_CHILDS 16

#define GET_CHILD(NODE, I) ((uintptr_t)NODE[I])
#define GET_SCALAR(NODE, I) (*((float *)(NODE) + I))
#define GET_SCALAR_INT(NODE, I) (*((uint32_t *)NODE + I))
#define IS_TERM(PRIM, ARITY) (!ARITY || PRIM == Scalar)
#define IS_FUNC(PRIM, ARITY) (ARITY && PRIM != Scalar)

// TODO: apply to other functions
#define PRIM_INDEX_TYPE uint32_t
#define FIT_TYPE float


#define DAG_EQ(DAG, ARITY, PRIM, CHILDREN)                                                              \
    (DAG->primitive == PRIM &&                                                                          \
        ((PRIM != Scalar && !memcmp(DAG->children.n, CHILDREN.n, ARITY * sizeof(*CHILDREN.n))) ||       \
        (PRIM == Scalar && !memcmp(DAG->children.f, CHILDREN.f, ARITY * sizeof(*CHILDREN.f)))))         \


#define IMPLEMENT_STATS(NODE, STATS, DEP, FUNC, TERM) do {      \
    if (TERM) {                                                 \
        STATS->depth = FUNC(DEP, STATS->depth);                 \
        stats->n_terms++;                                       \
    } else {                                                    \
        stats->n_prims++;                                       \
    }                                                           \
} while(0)


#define GENERATE_SCALAR(TABLE, ARITY, CHILDREN) do {                \
    ARITY = primitive_set[Scalar].arity;                            \
    float *scalar_ptr = alloc_scalar_ptrs(TABLE, ARITY);            \
    CHILDREN.f = scalar_ptr;                                        \
    for (int i = 0; i < ARITY; ++i) {                               \
        *scalar_ptr++ = (rand_float() * DELTA_SCALAR) + MIN_SCALAR; \
    }                                                               \
}while(0);

#define GENERATE_TERMINAL(TABLE, DUMMY) generate_program(TABLE, Full, 0, 0, 0.5, &DUMMY);


#define STEP_DOMAIN(RUN, CUR_VARS) do {                                                                 \
    int index = 0;                                                                                      \
    CUR_VARS[0] += RUN->STEP[0];                                                                        \
    while((CUR_VARS[index] > RUN->MAX_DOMAIN[index] + RUN->STEP[index] / 2.0) && (index + 1 < DIMS)) {  \
        CUR_VARS[index] = RUN->MIN_DOMAIN[index];                                                       \
        ++index;                                                                                        \
        CUR_VARS[index] += RUN->STEP[index];                                                            \
    }                                                                                                   \
} while(0)


#define PDEBUG 0
#if PDEBUG
    #define printfd(...) printf(__VA_ARGS__);
#else
    #define printfd(...)
#endif



// 16 bytes
typedef struct Prim{
    char *name;
    int arity;
    int id;
} Prim;

// 8 bytes
typedef union c_node {
    float *f;
    struct dag_node **n;
} c_node;

// 24 bytes
typedef struct dag_node {
    c_node children;
    uint16_t primitive; // also change in prim after
    short arity;
    int frequency;
    int n_offspring;
    int candid;
} dag_node;

// 32 bytes
typedef struct tree {
    dag_node *dag;
    float fitness;
    int n_prims;
    int n_terms;
    int depth;
    uint32_t count;
    short swap; // to facilitate the copy of old index tables
    short mode;
} tree;


enum Gen_Methods {Grow = 0, Full = 1, Ramped = 2};
enum Cnt_Methods {Primitives = 0, Terminals = 1, All_Prims = 2};
enum Ind_Arrays {Average = 0, StandardDev = 1, Best = 2};

enum PSet_IDs {Scalar = 1,
                X = 2,
                Y = 3,
                Add = 4, 
                Sub = 5,
                Div = 6,
                Mul = 7,
                Cos = 8,
                Sin = 9,
                Tan = 10,
                If = 11,
                PSET_START = 1,
                TSET_END = 3,
                PSET_END = 11
                };

extern const Prim primitive_set[];


#include "hashmap.h"
#include "engine.h"
#include "genetics.h"


void print_dag_node(dag_node *node);
void fancy_dag_print(dag_node *t);
void fancy_dag_print_w(dag_node *t, int dep, char *pre);
tree *create_tree();
tree *generate_tree(HashTable *t, int method, int min_depth, int max_depth, float term_prob);
dag_node *generate_program(HashTable *table, int method, int max_depth, int min_depth, float term_prob, tree *stats);
dag_node *generate_program_w(HashTable *table, int method, int depth, int limit_dep, float term_prob, tree *stats);
tree *str_to_tree(HashTable *table, const char* str);
dag_node *str_to_dag(HashTable *table, const char **strp, int depth, tree *stats);
Prim *get_prim_index(const char *prim_str);
PRIM_INDEX_TYPE *get_list_same_arity(int arity, int *n, int exclude, int ignore_specific, int search_set);
int get_prim_same_arity(int arity, int exclude, int ignore_specific, int search_set);
char *get_dag_expr(dag_node *t);
void print_tree(tree *t, int print_stats, int do_fancy);
void print_population(tree **population, int n, int generation, int do_fancy);
dag_node *get_dag_node_cnt(dag_node *node, uint32_t *count, const int mode);
int get_dag_node_dep(dag_node *node, int dep);
dag_node *copy_dag_node(dag_node *node, HashTable *to_table);
tree *copy_tree(tree *from, HashTable *to_table);
tree **generate_population(HashTable *t, int method, int min_depth, int max_depth, int n, float term_prob);

int cmp_tree(const void *a, const void *b);
tree **sort_k_min_trees(tree **array, int n, int size);
tree **get_k_min_trees(tree **arr, int n, int size);
tree **sel_k_min_trees(tree **arr, int n, int size);

void print_domain(FIT_TYPE *data, uint32_t n, const uint32_t *format);
FIT_TYPE icompute_node_g(dag_node *t, const FIT_TYPE* cur_vars);
FIT_TYPE icompute_node_g_d(dag_node *t, const FIT_TYPE* cur_vars);
int calc_pop_fit(Engine *run, tree **population);
FIT_TYPE no_tree_fit(tree *t);
FIT_TYPE calc_tree_fit(Engine *run, tree *t);
FIT_TYPE *icompute_domain_g(Engine *run, dag_node *t);
FIT_TYPE icalculate_fitness_g(Engine *run, FIT_TYPE *values);


#endif