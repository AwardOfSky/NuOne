#ifndef GENETICS_H
#define GENETICS_H

#define KOZA_RULE 0.9
#define SUBTRACT_CNT(TERM, MODE) (MODE == All_Prims || (!TERM && MODE == Primitives) || (TERM && MODE == Terminals))
#define MIN_MUT_DEP 1
#define MAX_MUT_DEP 5
#define DEFINED_MUTATIONS 4

#define COPY_NODE(ARITY, PRIM, NODE, TABLE, CHILDREN, RECURSE) do {     \
    if (ARITY > 0) {                                                    \
        if (PRIM != Scalar) {                                           \
            dag_node **child_ptr = alloc_child_ptrs(TABLE, ARITY);      \
            CHILDREN.n = child_ptr;                                     \
            for (int i = 0; i < ARITY; ++i) {                           \
                *child_ptr++ = RECURSE;                                 \
            }                                                           \
        } else {                                                        \
            CHILDREN.f = alloc_scalar_ptrs(TABLE, ARITY);               \
            memcpy(CHILDREN.f, NODE->CHILDREN.f, ARITY * sizeof(float));\
        }                                                               \
    }                                                                   \
} while(0)

dag_node *rnd_node_depn(dag_node *root, int max_dep, int n);
tree *subtree_crossover(tree *p1, tree *p2, HashTable *to_table, const int max_dep);
tree *subtree_crossover_d(tree *p1, tree *p2, HashTable *to_table, int c1, int c2, int m1, int m2, const int max_dep);
dag_node *copy_subtree_crossover(dag_node *p1, tree *p2, HashTable *to_table, tree *stats, int dep, const int max_dep);
tree *subtree_mutation(tree *t, HashTable *to_table, const int max_dep);
tree *subtree_mutation_d(tree *t, HashTable *to_table, int c, int m, const int max_dep);
dag_node *copy_subtree_mutation(dag_node *node, HashTable *to_table, tree *stats, int dep, const int max_dep);
tree *delete_mutation(tree *t, HashTable *to_table);
tree *delete_mutation_d(tree *t, HashTable *to_table, uint32_t c);
dag_node *copy_delete_mutation(dag_node *node, HashTable *to_table, tree *stats, int dep);
tree *insert_mutation(tree *t, HashTable *to_table, int max_dep);
tree *insert_mutation_d(tree *t, HashTable *to_table, uint32_t c, int m, int max_dep);
dag_node *copy_insert_mutation(dag_node *node, HashTable *to_table, tree *stats, int dep, int max_dep);
tree *point_mutation(tree *t, HashTable *to_table);
tree *point_mutation_d(tree *t, HashTable *to_table, uint32_t c, int m);
dag_node *copy_point_mutation(dag_node *node, HashTable *to_table, tree *stats, int dep);
tree* mutation(tree* parent, HashTable *t, int max_dep);
tree *tournament(Engine *run, tree **population);

#endif