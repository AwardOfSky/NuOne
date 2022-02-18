#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "gp.h"


const Prim primitive_set[] = {
    {"", 0, 0}, // placeholder for 0 ID
    {"scalar", 3, Scalar},
    {"x", 0, X},
    {"y", 0, Y},
    {"add", 2, Add},
    {"sub", 2, Sub},
    {"div", 2, Div},
    {"mul", 2, Mul}
};


void print_dag_node(dag_node *node) {
    int arity = node->arity;
    uint32_t primitive = node->primitive;

    if (primitive != 0xFFFF) {
        printf("Primitive: %s, Frequency: %d, Arity: %d, Nodes: %d, Candid: %d, Mem: %p, ",
                primitive_set[primitive].name, node->frequency, arity, node->n_offspring, node->candid, node);
    } else {
        printf("Primitive: [PLACEHOLDER], Frequency: %d, Arity: %d, Nodes: %d, Candid: %d, Mem: %p, ",
                node->frequency, arity, node->n_offspring, node->candid, node);
    }

    printf("Children: [");
    for(int k = 0; k < arity; ++k) {
        if (primitive != Scalar) {
            printf("%p", node->children.n[k]);
        } else {
            printf("%.3f", GET_SCALAR(node->children.f, k));
        }
        if (k < arity - 1) {
            printf(", ");
        }
    }
    printf("]\n");
}



void fancy_dag_print(dag_node *t) {
    char *temp_arr = (char *)calloc(64, sizeof(char));
    fancy_dag_print_w(t, 0, temp_arr);
    free(temp_arr);
}


void fancy_dag_print_w(dag_node *t, int dep, char *pre) {
    uint32_t primitive = t->primitive;
    int arity = t->arity;
    
    if (primitive != 0xFFFF) {
        printf("%s%s\n", pre, primitive_set[primitive].name);
    } else {
        printf("%s%s\n", pre, PREERR"This is a placeholder node, should not appear in a tree!");
    }
    
    pre[dep << 1] = '.';
    pre[(dep << 1) + 1] = '.';
    if (arity > 0) {
        if (primitive != Scalar) {
            for(int i = 0; i < arity; ++i) {
                //printf("ai jesus 0.5\n");
                fancy_dag_print_w(t->children.n[i], dep + 1, pre);
            }
        } else {
            for(int i = 0; i < arity; ++i) {
                printf("%s%.3f\n", pre, t->children.f[i]);
            }
        }
    }
    pre[dep << 1] = '\0';
}


tree *create_tree() {
    tree *t = (tree *)malloc(sizeof(tree));
    t->dag = NULL;
    t->depth = 0;
    t->fitness = FLT_MAX;
    t->n_terms = 0;
    t->n_prims = 0;
    t->swap = 0;
    t->count = UINT_MAX;
    return t;
}


tree *generate_tree(HashTable *t, int method, int min_depth, int max_depth, float term_prob) {
    tree *res = create_tree();
    res->dag = generate_program(t, method, max_depth, min_depth, term_prob, res);
    return res;
}


dag_node *generate_program(HashTable *table, int method, int max_depth, int min_depth, float term_prob, tree *stats) {
    stats->depth = INT_MAX;
    dag_node *res = generate_program_w(table, method, max_depth, max_depth - min_depth, term_prob, stats);
    stats->depth = max_depth - stats->depth;

    if (res == NULL) {
        printf(PREERR"Generated a NULL Program!\n");
    }
    return res;
}


dag_node *generate_program_w(HashTable *table, int method, int depth, int limit_dep, float term_prob, tree *stats) {
    uint32_t primitive;
    int arity = 0;
    c_node children = {NULL};

    // generate terminal node
    if((depth == 0) || (stats->n_prims >= MAX_TREE_NODES) || ((method == Grow) && (depth <= limit_dep) && (rand_float() < term_prob))) {
        stats->n_terms++;
        primitive = (rand() % TSET_LEN) + 1;
        
        if (primitive == Scalar) {
            arity = primitive_set[primitive].arity;

            float *scalar_ptr = alloc_scalar_ptrs(table, arity);
            children.f = scalar_ptr;
            for (int i = 0; i < arity; ++i) {
                *scalar_ptr++ = (rand_float() * DELTA_SCALAR) + MIN_SCALAR;
            }
        }

        stats->depth = min(depth, stats->depth);
    // generate function node
    } else {
        stats->n_prims++;
        primitive = (rand() % FSET_LEN) + TSET_END + 1;
        arity = primitive_set[primitive].arity;

        dag_node **child_ptr = alloc_child_ptrs(table, arity);
        children.n = child_ptr;
        for(int i = 0; i < arity; ++i) {
            *child_ptr++ = generate_program_w(table, method, depth - 1, limit_dep, term_prob, stats);
        }
    }
    return add_dag_node(table, primitive, children);
}


tree *str_to_tree(HashTable *table, const char* str) {
    tree *result = create_tree();
    result->dag = str_to_dag(table, &str, 0, result);
    return result;
}


// By design this function ignores spaces, a primitive name should not have a space
#define SKIP_CHAR(STRING, CHR) while(*STRING == CHR) { ++STRING; }
dag_node *str_to_dag(HashTable *table, const char **strp, int depth, tree *stats) {

    const char *str = *strp;
    SKIP_CHAR(str, ' '); // skip chars leading to the prim

    const char *start_prim = str;
    
    while(*str != '(' && *str != ')' && *str != ',' && *str != ' ') {
        ++str; // read the prim
    }

    size_t prim_size = min(str - start_prim, MAX_PRIM_NAME_SIZE);
    char chr_buffer[MAX_PRIM_NAME_SIZE + 1];
    chr_buffer[prim_size] = '\0';
    memcpy(chr_buffer, start_prim, prim_size);

    Prim *p = get_prim_index(chr_buffer);
    uint32_t primitive = p->id;
    int arity = p->arity;
    c_node children = {NULL};

    while(*str != '(' && *str != ')' && *str != ',') {
        ++str;  // advance to delimiter
    }

    if(arity > 0) {
        if(primitive == Scalar) {

            float *scalar_ptr = alloc_scalar_ptrs(table, arity);
            children.f = scalar_ptr;

            int index = 0;
            while(*str != ')') {
                ++str; // ++ for ( or ,

                const char *num_start = str;
                while(*str != ')' && *str != ',') {
                    ++str;  // advance to delimiter
                }
                size_t num_size = str - num_start;

                chr_buffer[num_size] = '\0';
                memcpy(chr_buffer, num_start, num_size);
                *scalar_ptr++ = atof(chr_buffer);

                ++index;
            }

            float value_to_rep = *(scalar_ptr - 1);
            for(int i = index; i < arity; ++i) { // extend
                *scalar_ptr++ = value_to_rep;
            }

        } else {

            dag_node **child_ptr = alloc_child_ptrs(table, arity);
            children.n = child_ptr;
            for(int i = 0; i < arity; ++i) {
                ++str; // pass delimiter
                *child_ptr++ = str_to_dag(table, &str, depth + 1, stats);
            }

        }
        ++str; // ++ for )
        SKIP_CHAR(str, ' ');

    }

    *strp = str;
    dag_node *res_ptr = add_dag_node(table, primitive, children);
    IMPLEMENT_STATS(res_ptr, stats, depth, max, IS_TERM(primitive, arity));
    return res_ptr;
}


// primitive can be set to -1 if you want to allow all primitives of that arity
int get_prim_same_arity(uint32_t primitive, int arity) {
    if (primitive == Scalar) {
        return Scalar;
    }

#define PRIM_TYPE int

    PRIM_TYPE *candidates;
    if (PSET_LEN < STACK_ALLOC) {
        candidates = (PRIM_TYPE *)alloca(PSET_LEN * sizeof(*candidates));
    } else {
        candidates = (PRIM_TYPE *)malloc(PSET_LEN * sizeof(*candidates));
    }

    int index = 0;
    int term = IS_TERM(primitive, arity);
    PRIM_TYPE start, end;
    if (term) {
        start = PSET_START;
        end = TSET_END;
    } else {
        start = TSET_END;
        end = PSET_END;
    }

    for(PRIM_TYPE i = start; i <= end; ++i) {
        if (primitive_set[i].arity == arity && i != primitive) {
            candidates[index++] = i;
        }
    }

    int chosen_one = candidates[(rand() % index)];
    if(PSET_LEN >= STACK_ALLOC) {
        free(candidates);
    }

    return chosen_one;
    #undef PRIM_TYPE
}


Prim *get_prim_index(const char *prim_str) {
    Prim* p;
    for(p = (Prim *)(primitive_set + 1); p < primitive_set + PSET_LEN && strcmp(prim_str, p->name); ++p);
    return p;
}


char *get_dag_expr(dag_node *t) {
    char *pname;
    if (t->primitive != 0xFFFF) {
        asprintf(&pname, "%s", primitive_set[t->primitive].name);
    } else {
        asprintf(&pname, "%s", PREERR"This is a placeholder node, should not appear in a tree!");
    }

    if (t->arity) {
        asprintf(&pname, "%s(", pname);
        if(t->primitive != Scalar) {
            //PRINT_CHILD("%s", print_dag_tree(GET_CHILD(t, i)), i);
            for(int i = 0; i < t->arity; ++i) {
                if(t->children.n[i] != NULL) {
                    asprintf(&pname, "%s%s", pname, get_dag_expr(t->children.n[i]));
                    if(i < t->arity - 1) asprintf(&pname, "%s, ", pname);
                }
            }
        } else {
            //PRINT_CHILD("%.3f", GET_SCALAR(t, i), 0);
            for(int i = 0; i < t->arity; ++i) {
                if(t->children.f != NULL) {
                    asprintf(&pname, "%s%.3f", pname, GET_SCALAR(t->children.f, i));
                    if(i < t->arity - 1) asprintf(&pname, "%s, ", pname);
                }
            }    
        }
        asprintf(&pname, "%s)", pname);
    }
    return pname;
}


void print_tree(tree *t, int print_stats, int do_fancy) {
    if (t != NULL) {

        if (t->dag != NULL) {
            //printf("Hey!\n");
            printf("\t%s\n", get_dag_expr(t->dag));
        } else {
            printf(PREWARN"This tree has no dag.\n");
        }
        if (do_fancy) {
            printf("In fancy print:\n");
            fancy_dag_print(t->dag);
        }

        if (print_stats) {
            printf("Individual stats:\n");
            int terms = t->n_terms;
            int prims = t->n_prims;
            printf("(Terminals, Primitives, total_nodes, depth, fitness): %d, %d, %d, %d, %.3f\n", 
            terms, prims, terms + prims, t->depth, t->fitness);
        }
    } else {
        printf("Tree is Null.\n");
    }
}


dag_node *get_dag_node_cnt(dag_node *node, uint32_t *count, int mode) {
    dag_node *res_dag = node;
    *count -= SUBTRACT_CNT(node, mode);

    if (!*count) {
        return node;
    }

    if (node->primitive != Scalar) {
        for (int i = 0; i < node->arity && *count; ++i) {
            res_dag = get_dag_node_cnt(node->children.n[i], count, mode);
        }
    }

    return res_dag;
}


int get_dag_node_dep(dag_node *node, int dep) {
    int max_dep = dep;
    if (node->primitive != Scalar) {
        for(int i = 0; i < node->arity; ++i) {
            max_dep = max(max_dep, get_dag_node_dep(node->children.n[i], dep + 1));
        }
    }
    return max_dep;
}


dag_node *copy_dag_node(dag_node *node, HashTable *to_table) {
    uint32_t primitive = node->primitive;
    int arity = node->arity;
    c_node children = {NULL};

    COPY_NODE(arity, primitive, node, to_table, children, copy_dag_node(node->children.n[i], to_table));
    
    return add_dag_node_index(to_table, primitive, children, node->candid);    
}


tree *copy_tree(tree *from, HashTable *to_table) {
    tree *result = create_tree();
    result->dag = copy_dag_node(from->dag, to_table);
    result->depth = from->depth;
    result->n_terms = from->n_terms;
    result->n_prims = from->n_prims;
    return result;
}


tree **generate_population(HashTable *t, int method, int min_depth, int max_depth, int n, float term_prob) {
    tree **population = (tree **)malloc(sizeof(tree *) * n);
    
    if (min_depth < 0) {
        printf(PREERR"Generate Population: Min depth (%d) lower than 0, setting to 0;\n", min_depth);
        min_depth = 0;
    }
    if (max_depth < 0) {
        printf(PREERR"Generate Population: Max depth (%d) lower than 0, setting to 0;\n", max_depth);
        max_depth = 0;
    }
    if (max_depth < min_depth) {
        printf(PREERR"Generate Population: Min depth (%d) lower than max depth (%d), switching.\n", min_depth, max_depth);
        int temp = max_depth;
        max_depth = min_depth;
        min_depth = temp;
    }

    if(method != Ramped) {
        for(int i = 0; i < n; ++i) {
            population[i] = generate_tree(t, method, min_depth, max_depth, term_prob);
        }
    } else {
        int divisions = max_depth - (min_depth - 1);
        int parts = floor((float)n / (float)divisions);
        int last_part = n - (divisions - 1) * parts;
        int load_balance_index = (max_depth + 1) - (last_part - parts);

        int num_parts = parts;
        int mfull = floor((float)num_parts / 2.0f);

        int pop_index = 0;
        for(int i = min_depth; i <= max_depth; ++i) {
            if (i == load_balance_index) {
                num_parts += 1;
                mfull += 1;
            }

            int cur_method = Full;
            for(int j = 0; j < num_parts; ++j) {
                if(j >= mfull) {
                    cur_method = 1 - cur_method; // switch between GROW and FULL
                }
                //printf("Generating ind with method %d, dep %d\n", cur_method, i);
                population[pop_index] = generate_tree(t, cur_method, min_depth, i, term_prob);
                ++pop_index;
            }
        }
    }

    return population;
}