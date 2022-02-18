#include <stdio.h>
#include <stdlib.h>
#include "gp.h"

tree *subtree_crossover(tree *p1, tree *p2, HashTable *to_table, const int max_dep) {
    int mode_p1 = (rand_float() < KOZA_RULE) ? Primitives : Terminals;
    int mode_p2 = (rand_float() < KOZA_RULE) ? Primitives : Terminals;
    uint32_t count_p1 = (mode_p1 == Primitives && p1->n_prims > 0) ? ((rand() % p1->n_prims) + 1) : ((rand() % p1->n_terms) + 1);
    uint32_t count_p2 = (mode_p2 == Primitives && p2->n_prims > 0) ? ((rand() % p2->n_prims) + 1) : ((rand() % p2->n_terms) + 1);
    
    //printf("cross (c1, c2, m1, m2, mdep) %d %d %d %d %d\n", count_p1, count_p2, mode_p1, mode_p2, UPPER_MAX_DEP);
    return subtree_crossover_d(p1, p2, to_table, count_p1, count_p2, mode_p1, mode_p2, max_dep);
}


tree *subtree_crossover_d(tree *p1, tree *p2, HashTable *to_table, int c1, int c2, int m1, int m2, const int max_dep) {
    tree *res_tree = create_tree();
 
    if (p1->n_prims + p1->n_terms > 0 && p2->n_prims + p2->n_terms > 0) {
        int mode_p2 = m2;
        uint32_t count_p2 = c2;
        tree *subtree_in_p2 = create_tree();
        subtree_in_p2->dag = get_dag_node_cnt(p2->dag, &count_p2, mode_p2);

        if (subtree_in_p2->dag != NULL) {
            int subtree_dep = get_dag_node_dep(subtree_in_p2->dag, 0);
            subtree_in_p2->depth = subtree_dep;
            //printf("Depth of subtree in p2: %d\n", subtree_dep);
            //fancy_dag_print(subtree_in_p2->dag);

            //fancy_dag_print(subtree_in_p2); //printf("This i sp2\n");
            res_tree->mode = m1;
            res_tree->count = c1;
            res_tree->dag = copy_subtree_crossover(p1->dag, subtree_in_p2, to_table, res_tree, 0, max_dep);
        } else {
            printf(PREERR"Failed crossover operation: failed to get subtree in second parent.\n");
            free(res_tree);
            return copy_tree(p1, to_table);
        }
    } else {
        printf(PREERR"Failed crossover operation: one of the trees is empty.\n");
    }

    return res_tree;
}


// return random node at depth n
// by default, the function does not return terminals unless if returning nodes at last depth (max_dep)
dag_node *rnd_node_depn(dag_node *root, int max_dep, int n) {
    dag_node *res = root;
    dag_node **prim_childs;

    if (MAX_CHILDS <= STACK_ALLOC) {
        prim_childs = (dag_node **)alloca(sizeof(*prim_childs) * MAX_CHILDS);
    } else {
        prim_childs = (dag_node **)malloc(sizeof(*prim_childs) * MAX_CHILDS);
    }

    for(int i = 0; i < n; ++i) {
        
        int n_prim_childs = 0;
        for(int j = 0; j < res->arity; ++j) {
            dag_node *c = res->children.n[j];

            if (IS_FUNC(c->primitive, c->arity) || i == max_dep - 1) {
                prim_childs[n_prim_childs++] = c;
            }
        }
        
        if (n_prim_childs > 0) {
            res = prim_childs[rand() % n_prim_childs];
        }
    }

    if (MAX_CHILDS > STACK_ALLOC) {
        free(prim_childs);
    }

    return res;
}


dag_node *copy_subtree_crossover(dag_node *p1, tree *p2, HashTable *to_table, tree *stats, int dep, const int max_dep) {
    c_node children = {NULL};
    int term = IS_TERM(p1->primitive, p1->arity);
    stats->count -= SUBTRACT_CNT(term, stats->mode);

    if (!stats->count) {
        stats->count = UINT_MAX;
        stats->swap++;

        dag_node *p2_dag = p2->dag;

        int dep_to_take = p2->depth + dep - max_dep;
        if (dep_to_take > 0) {
            p2_dag = rnd_node_depn(p2_dag, p2->depth, dep_to_take);
        }
        
        p1 = p2_dag;
    }
    uint32_t primitive = p1->primitive;
    int arity = p1->arity;

    int old_swap = stats->swap;
    COPY_NODE(arity, primitive, p1, to_table, children,
                copy_subtree_crossover(p1->children.n[i], p2, to_table, stats, dep + 1, max_dep));
    IMPLEMENT_STATS(p1, stats, dep, max, IS_TERM(primitive, arity));

    int ind = (stats->swap == old_swap) ? p1->candid : -1;
    return add_dag_node_index(to_table, primitive, children, ind);   
}


tree *subtree_mutation(tree *t, HashTable *to_table, const int max_dep) {
    return subtree_mutation_d(t, to_table, (rand() % (t->n_prims + t->n_terms)) + 1, All_Prims, max_dep);
}


tree *subtree_mutation_d(tree *t, HashTable *to_table, int c, int m, const int max_dep) {
    tree *res_tree = create_tree();
    res_tree->count = c;
    res_tree->mode = m;

    res_tree->dag = copy_subtree_mutation(t->dag, to_table, res_tree, 0, max_dep);
    return res_tree;
}


dag_node *copy_subtree_mutation(dag_node *node, HashTable *to_table, tree *stats, int dep, const int max_dep) {
    uint32_t primitive = node->primitive;
    int arity = node->arity;

    int term = IS_TERM(primitive, arity);
    stats->count -= SUBTRACT_CNT(term, stats->mode);

    dag_node *res_ptr = NULL;
    if (stats->count != 0) {
        c_node children = {NULL};
        int old_swap = stats->swap;

        COPY_NODE(arity, primitive, node, to_table, children,
                copy_subtree_mutation(node->children.n[i], to_table, stats, dep + 1, max_dep));
        int ind = (stats->swap == old_swap) ? node->candid : -1;

        res_ptr = add_dag_node_index(to_table, primitive, children, ind);
        
        IMPLEMENT_STATS(node, stats, dep, max, term);
    } else {
        stats->count = UINT_MAX;
        stats->swap++;

        tree *st = create_tree();

        //make sure tree doesnt pass overall depth limit
        int max_gen_dep = clip(0, rand() % (MAX_MUT_DEP + 1) + MIN_MUT_DEP, max_dep - dep);
        res_ptr = generate_program(to_table, Grow, max_gen_dep, min(MIN_MUT_DEP, max_gen_dep), 0.5, st);

        if (res_ptr == NULL) {
            printf(PREERR"Failed to generate child (subtree).\n");
        }

        stats->depth = max(stats->depth, dep + st->depth);
        stats->n_prims += st->n_prims;
        stats->n_terms += st->n_terms;
        free(st);
    }

    return res_ptr;
}

