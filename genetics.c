#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "gp.h"

tree *subtree_crossover(tree *p1, tree *p2, HashTable *to_table, const int max_dep) {
    int mode_p1 = (rand_float() < KOZA_RULE) ? Primitives : Terminals;
    int mode_p2 = (rand_float() < KOZA_RULE) ? Primitives : Terminals;
    uint32_t count_p1 = (mode_p1 == Primitives && p1->n_prims > 0) ? ((rand() % p1->n_prims) + 1) : ((rand() % p1->n_terms) + 1);
    uint32_t count_p2 = (mode_p2 == Primitives && p2->n_prims > 0) ? ((rand() % p2->n_prims) + 1) : ((rand() % p2->n_terms) + 1);
    
    printfd("cross (cp1, cp2, m1, m2, mdep) %d %d %d %d %d\n", count_p1, count_p2, mode_p1, mode_p2, max_dep);
    return subtree_crossover_d(p1, p2, to_table, count_p1, count_p2, mode_p1, mode_p2, max_dep);
}


tree *subtree_crossover_d(tree *p1, tree *p2, HashTable *to_table, int c1, int c2, int m1, int m2, const int max_dep) {
    tree *res_tree = create_tree();
 
    if (p1->n_prims + p1->n_terms > 0 && p2->n_prims + p2->n_terms > 0) {
        int mode_p2 = m2;
        uint32_t count_p2 = c2;
        tree *subtree_in_p2 = create_tree();

        printfd("Tree in p2: %s\n", get_dag_expr(p2->dag));
        printfd("(cp2, m2): (%d, %d)\n", count_p2, mode_p2);
        subtree_in_p2->dag = get_dag_node_cnt(p2->dag, &count_p2, mode_p2);

        if (subtree_in_p2->dag != NULL) {
            int subtree_dep = get_dag_node_dep(subtree_in_p2->dag, 0);
            subtree_in_p2->depth = subtree_dep;
            
            printfd("Depth of subtree in p2: %d\n", subtree_dep);
            printfd("Subtree in p2: %s\n", get_dag_expr(subtree_in_p2->dag));
            
            //fancy_dag_print(subtree_in_p2->dag);

            res_tree->mode = m1;
            res_tree->count = c1;
            res_tree->dag = copy_subtree_crossover(p1->dag, subtree_in_p2, to_table, res_tree, 0, max(max_dep, p1->depth));
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
        //printf("dep to take, dep, max_d, %d %d %d\n", dep_to_take, dep, max_dep);
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

    res_tree->dag = copy_subtree_mutation(t->dag, to_table, res_tree, 0, max(max_dep, t->depth));
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


tree *delete_mutation(tree *t, HashTable *to_table) {
    return delete_mutation_d(t, to_table, (rand() % t->n_prims) + 1);
}


tree *delete_mutation_d(tree *t, HashTable *to_table, uint32_t c) {
    tree *res_tree = create_tree();

    if (t->n_prims > 0) {
        res_tree->count = c;
        res_tree->dag = copy_delete_mutation(t->dag, to_table, res_tree, 0);
    } else {
        printf(PREWARN"Cannot perform delete mutation: tree \"%s\" only has one node.\n", get_dag_expr(t->dag));
        res_tree = copy_tree(t, to_table);
    }

    return res_tree;
}


dag_node *copy_delete_mutation(dag_node *node, HashTable *to_table, tree *stats, int dep) {
    int primitive = node->primitive;
    int arity = node->arity;

    c_node children = {NULL};

    int func = IS_FUNC(primitive, arity);
    stats->count -= func;

    if (!stats->count) {
        //printf("Chosen promotion prim: %s\n", primitive_set[primitive].name);
        stats->count = UINT_MAX;
        stats->swap++;

        int random_prim = (rand() % arity);
        node = node->children.n[random_prim];
        primitive = node->primitive;
        arity = node->arity;
    }

    int old_swap = stats->swap;
    COPY_NODE(arity, primitive, node, to_table, children,
                copy_delete_mutation(node->children.n[i], to_table, stats, dep + 1));

    IMPLEMENT_STATS(node, stats, dep, max, !func);
    int ind = (stats->swap == old_swap) ? node->candid : -1;

    return add_dag_node_index(to_table, primitive, children, ind);
}


tree *insert_mutation(tree *t, HashTable *to_table, const int max_dep) {
    return insert_mutation_d(t, to_table, (rand() % (t->n_prims + t->n_terms)) + 1, All_Prims, max_dep);
}


tree *insert_mutation_d(tree *t, HashTable *to_table, uint32_t c, int m, const int max_dep) {
    tree *res_tree = create_tree();
    int p = get_prim_same_arity(0, 1, -1, 1);

    if (p != -1) {
        res_tree->count = c;
        res_tree->mode = m;
        res_tree->dag = copy_insert_mutation(t->dag, to_table, res_tree, 0, max(max_dep, t->depth));
    } else {
        printf(PREERR"Cannot perform insert mutation: there is no function primitive with more than 0 arguments\n");
        res_tree = copy_tree(t, to_table);
    }

    return res_tree;
}


dag_node *copy_insert_mutation(dag_node *node, HashTable *to_table, tree *stats, int dep, const int max_dep) {
    int primitive = node->primitive;
    int arity = node->arity;
    c_node children = {NULL};
    dag_node *res_ptr = NULL;

    int term = IS_TERM(primitive, arity);
    stats->count -= SUBTRACT_CNT(term, stats->mode);
    
    if ((stats->count > 0) || (dep == max_dep)) {
        int old_swap = stats->swap;
        if (dep == max_dep && !term) {
            int random_prim = (rand() % arity);
            node = node->children.n[random_prim];
            primitive = node->primitive;
            arity = node->arity;
            stats->depth = dep;
        }
        COPY_NODE(arity, primitive, node, to_table, children,
                copy_insert_mutation(node->children.n[i], to_table, stats, dep + 1, max_dep));
        IMPLEMENT_STATS(node, stats, dep, max, term);
        int ind = (stats->swap == old_swap) ? node->candid : -1;

        res_ptr = add_dag_node_index(to_table, primitive, children, ind);
    } else {
        stats->count = UINT_MAX;
        stats->swap++;
        
        // chose a function prim that does not have 0 arguments
        primitive = get_prim_same_arity(0, 1, -1, 1);
        int new_arity = primitive_set[primitive].arity;
        int node_slot = rand() % new_arity;

        dag_node **child_ptr = alloc_child_ptrs(to_table, new_arity);
        children.n = child_ptr;
        for(int i = 0; i < new_arity; ++i) {
            if (i == node_slot) {
                *child_ptr++ = copy_insert_mutation(node, to_table, stats, dep + 1, max_dep);
            } else {
                tree dummy_tree = {0};
                dag_node *term_node = GENERATE_TERMINAL(to_table, dummy_tree);
                //fancy_dag_print(term_node);
                *child_ptr++ = term_node;
            }
        }

        stats->n_prims++;
        stats->n_terms += new_arity - 1;
        res_ptr = add_dag_node(to_table, primitive, children);
    }

    return res_ptr;
}


tree *point_mutation(tree *t, HashTable *to_table) {
    return point_mutation_d(t, to_table, (rand() % (t->n_prims + t->n_terms)) + 1, All_Prims);
}


tree *point_mutation_d(tree *t, HashTable *to_table, uint32_t c, int m) {
    tree *res_tree = create_tree();
    res_tree->mode = m;
    res_tree->count = c;

    res_tree->dag = copy_point_mutation(t->dag, to_table, res_tree, 0);
    return res_tree;
}


dag_node *copy_point_mutation(dag_node *node, HashTable *to_table, tree *stats, int dep) {
    int primitive = node->primitive;
    int arity = node->arity;
    c_node children = {NULL};    

    int term = IS_TERM(primitive, arity);
    stats->count -= SUBTRACT_CNT(term, stats->mode);

    int old_swap = stats->swap;
    if (stats->count == 0) {
        stats->count = UINT_MAX;
        stats->swap++;

        if (term) {
            if (primitive == Scalar) {
                GENERATE_SCALAR(to_table, arity, children);
            } else {
                int ptemp = get_prim_same_arity(primitive_set[Scalar].arity, 1, primitive, 0);
                printf("got Prim: %d\n", ptemp);
                if (ptemp != -1) {
                    primitive = ptemp;
                }
            }
            arity = 0; // COPY_NODE without effect

        } else {
            int ptemp = get_prim_same_arity(-1, 1, primitive, 1);

            if (ptemp != -1) {
                primitive = ptemp;
                int new_arity = primitive_set[primitive].arity;

                if (new_arity <= arity) {
                    arity = new_arity;
                } else {
                    dag_node **child_ptr = alloc_child_ptrs(to_table, new_arity);
                    children.n = child_ptr;
                    for (int i = 0; i < arity; ++i) {
                        *child_ptr++ = copy_point_mutation(node->children.n[i], to_table, stats, dep + 1);
                    }
                    tree dummy_tree = {0};
                    for (int i = arity; i < new_arity; ++i) {
                        *child_ptr++ = GENERATE_TERMINAL(to_table, dummy_tree);
                    }
                    stats->n_terms += new_arity - arity;
                    arity = 0; // COPY_NODE without effect
                }
            }
        }
    }

    COPY_NODE(arity, primitive, node, to_table, children,
                copy_point_mutation(node->children.n[i], to_table, stats, dep + 1));

    IMPLEMENT_STATS(node, stats, dep, max, term);
    int ind = (stats->swap == old_swap) ? node->candid : -1;

    return add_dag_node_index(to_table, primitive, children, ind);
}


tree* mutation(tree* parent, HashTable *t, const int max_dep) {
    int mutation_selector = rand() % DEFINED_MUTATIONS;
    switch(mutation_selector) {
        default:
        case 0: return subtree_mutation(parent, t, max_dep);
        case 1: return point_mutation(parent, t);
        case 2: return delete_mutation(parent, t);
        case 3: return insert_mutation(parent, t, max_dep);
    }
}


tree *tournament(Engine *run, tree **population) {
    #define TS_BEST_IN_POP(IND) do {                        \
        float best_fit = FLT_MAX;                           \
        for(int i = 0; i < size; ++i) {                     \
            float fitness = population[IND]->fitness;       \
            if(best_fit > fitness) {                        \
                best_fit = fitness;                         \
                best_ind = IND;                             \
            }                                               \
        }                                                   \
    } while(0)


    int n = run->pop_size;
    int size = run->tournament_size;
    if(size < 0) {
        printf(PREWARN"Tournament size below 0: %d, constraining size to 1.\n", size);
        size = 1;
    }
    
    int best_ind = -1;
    if (size < n) {
        int *indices = range_random_sample(n, size);
        //PRINT_ARR(indices, size, "%d");
        TS_BEST_IN_POP(indices[i]);
        free(indices);
    } else {
        printf(PREWARN"Tournament size above n(%d): %d, constraining size to n.\n", n, size);
        size = clip(0, size, n);
        TS_BEST_IN_POP(i);
    }

    if (best_ind <= -1) {
        printf(PREERR"Tournament: Bad best fit index: %d\n", best_ind);
    }
    
    return population[best_ind];
}