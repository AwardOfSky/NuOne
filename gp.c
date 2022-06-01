#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <limits.h>
#include <omp.h>
#include "gp.h"


const Prim primitive_set[] = {
    {"placeholder", 0, 0}, // placeholder for 0 ID
    {"scalar", 3, Scalar},
    {"x", 0, X},
    {"y", 0, Y},
    {"add", 2, Add},
    {"sub", 2, Sub},
    {"div", 2, Div},
    {"mul", 2, Mul},
    {"sin", 1, Cos},
    {"cos", 1, Sin},
    {"tan", 1, Tan},
    {"if", 3, If}
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
            GENERATE_SCALAR(table, arity, children);
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


Prim *get_prim_index(const char *prim_str) {
    Prim* p;
    for(p = (Prim *)(primitive_set + 1); p < primitive_set + PSET_LEN && strcmp(prim_str, p->name); ++p);
    return p;
}



char *get_dag_expr(const dag_node *t) {
    // This is not a stack alloc but STCK_ALLOC is a nice init size
    // for the purposes of building a string to store the individual's genotype
    uint32_t res_size = STACK_ALLOC;
    uint32_t pos = 0;
    char *res = (char *)malloc(sizeof(char) * res_size);
    get_dag_expr_w(t, &res, &pos, &res_size);
    return res;
}


void safe_strcat(char **strp, char *to_add, uint32_t *pos, uint32_t *sptr) {
    uint32_t size_to_add = strlen(to_add);
    uint32_t size_after_add = *pos + size_to_add + 1;
    if (size_after_add > *sptr) {
        uint32_t next_size = max(size_after_add, *sptr << 1);
        realloc_char(strp, next_size, sptr);
    }
    strncpy(*strp + *pos, to_add, size_to_add + 1);
    *pos += size_to_add;
}


void get_dag_expr_w(const dag_node *t, char **strp, uint32_t *pos, uint32_t *sptr) {
    if (t->primitive != 0xFFFF) {
        safe_strcat(strp, primitive_set[t->primitive].name, pos, sptr);
    } else {
        safe_strcat(strp, PREERR"This is a placeholder node, should not appear in a tree!", pos, sptr);
    }

    if (t->arity) {
        safe_strcat(strp, "(", pos, sptr);

        if(t->primitive != Scalar) {
            for(int i = 0; i < t->arity; ++i) {
                if(t->children.n[i] != NULL) {                    
                    get_dag_expr_w(t->children.n[i], strp, pos, sptr);
                    if(i < t->arity - 1) safe_strcat(strp, ", ", pos, sptr);
                }
            }
        } else {

            for(int i = 0; i < t->arity; ++i) {
                if(t->children.f != NULL) {
                    
                    char scalar_to_add[STACK_ALLOC];
                    sprintf(scalar_to_add, "%.3f", GET_SCALAR(t->children.f, i));
                    safe_strcat(strp, scalar_to_add, pos, sptr);

                    if(i < t->arity - 1) safe_strcat(strp, ", ", pos, sptr);
                }
            }    
        }
        safe_strcat(strp, ")", pos, sptr);
    }

}


void print_tree(tree *t, int print_stats, int do_fancy) {
    if (t != NULL) {

        if (t->dag != NULL) {
            //printf("Hey!\n");
            if (PDEBUG) {
                char *temp_str = get_dag_expr(t->dag);
                printfd("\t%s\n", temp_str);
                free(temp_str);
            }

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


void print_population(tree **population, int n, int generation, int do_fancy) {
    if(generation < 0) {
        printf("\nIntial population:\n\n");
    } else {
        printf("\nPopulation of generation %d:\n\n", generation);
    }
    for(int i = 0; i < n; ++i) {
        printf("Individual #%d:", i);
        print_tree(population[i], 0, do_fancy);
    }
}


dag_node *get_dag_node_cnt(dag_node *node, uint32_t *count, const int mode) {
    dag_node *res_dag = node;

    int term = IS_TERM(node->primitive, node->arity);
    *count -= SUBTRACT_CNT(term, mode);
    //printf("(Prim, count, var) : (%s %d %d)\n", primitive_set[node->primitive].name, *count, SUBTRACT_CNT(node, mode));

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


// "exclude" - exclude mode
// (0 -> include all primitives with that arity)
// (1 -> exclude all primitives with that arity)
// ignore_specific to ignore a specific primitive, ignore_specific = -1 if you do not wish to ignore anything
// 
// search_set - which sets to search on
// (0 -> terminal set, only search for terminals)
// (1 -> function set, only search for functions)
// (>=2 -> search for both functions and terminals)
prim_index_type *get_list_same_arity(int arity, int *n, int exclude, int ignore_specific, int search_set) {
    prim_index_type *candidates = (prim_index_type *)malloc(PSET_LEN * sizeof(*candidates));

    int index = 0;
    prim_index_type start, end;
    if (search_set == 0) {
        start = PSET_START;
        end = TSET_END;
    } else if (search_set == 1) {
        start = TSET_END + 1;
        end = PSET_END;
    } else {
        start = PSET_START;
        end = PSET_END;
    }

    for(prim_index_type i = start; i <= end; ++i) {
        int p = primitive_set[i].arity;
        if (((exclude == 0 && p == arity) ||
             (exclude == 1 && p != arity)) && ignore_specific != i) {
            //printf("adding candidate[%d]: %d\n", index, i);
            candidates[index++] = i;
        }
    }

    *n = index;
    return candidates;
}


prim_index_type get_prim_same_arity(int arity, int exclude, int ignore_specific, int search_set) {
    int n;
    prim_index_type *candidates = get_list_same_arity(arity, &n, exclude, ignore_specific, search_set);
    int chosen_one = (n > 0) ? candidates[rand() % n] : -1;
    free(candidates);
    return chosen_one;
}


int cmp_tree(const void *a, const void *b) {
    fit_t aval = (*(tree **)a)->fitness;
    fit_t bval = (*(tree **)b)->fitness;
    int res = (bval > aval) ? -1 : ((bval < aval) ? 1 : 0);
    //printf("Comparing %.3f > %.3f for result %d\n", aval, bval, res);
    return res;
}


tree **sort_k_min_trees(tree **array, int n, int size) {
    tree **k_min = (tree **)malloc(size * sizeof(*k_min));
    qsort(array, n, sizeof(tree *), cmp_tree);
    memcpy(k_min, array, size * sizeof(tree *));
    return k_min;
}


tree **get_k_min_trees(tree **arr, int n, int size) {
    if (size) {
        //min(min(0.02*n + 260, 0.06*n + 40.0), 60000)
        if (size < 0) {
            return sel_k_min_trees(arr, n, size);
        } else {
            return sort_k_min_trees(arr, n, size);
        }
    } else {
        printf(PREERR"Got empty list of trees to order!\n");
        return NULL;
    }
}


tree **sel_k_min_trees(tree **arr, int n, int size) {
    fit_t max_of_mins = FLT_MAX;
    tree **k = (tree **)malloc((size + 1) * sizeof(tree *));

    k[0] = arr[0];
    for(int i = 1; i < n; ++i) {
        int k_els = min(i, size);
        fit_t cur_val = arr[i]->fitness;
        
        if (cur_val < max_of_mins) {
            
            int j;
            if (cur_val > k[0]->fitness && cur_val < k[k_els - 1]->fitness) {
                if (k_els > BINSEARCH_TRESHOLD) {
                    int inc = k_els >> 1;
                    j = inc;
                    while (k[j]->fitness < cur_val || k[j - 1]->fitness > cur_val) {
                        inc = (inc > 1) ? (inc >> 1) : 1;
                        j += (cur_val < k[j]->fitness) ? -inc : inc;
                    }
                } else {
                    for(j = 0; j < k_els && k[j]->fitness < cur_val; ++j) {}
                }
            } else {
                if (cur_val <= k[0]->fitness) {
                    j = 0;
                } else {
                    j = k_els;
                }
            }

            memmove(k + j + 1, k + j, (size - j) * sizeof(*k));
            k[j] = arr[i];

            if (k_els == size) {
                max_of_mins = k[size - 1]->fitness; 
            }
        }
    }

    return k;
}


void print_domain(fit_t *data, uint32_t n, const uint32_t *format) {
    printf("\n%.3f ", data[0]);
    for(int i = 1; i < n; ++i) {
        int index = 0;
        int divi = format[0];
        while( (i % divi) == 0) {
            divi *= format[++index];
        }
        for(int j = 0; j < (1 << (index - 1)); ++j) {
            putchar('\n');
        }
        printf("%.3f ", data[i]);
    }
    printf("\n");
}


#define _CN_RESULT_VAR result
#define _CN_CHILDS_VAR childs
#define GP_RES _CN_RESULT_VAR

//    _Pragma("omp parallel");

#define COMPUTE_NODE_TEMPLATE(CODE) do {                                                \
    _CN_RESULT_VAR = (fit_t *)malloc(sizeof(*_CN_RESULT_VAR) * run.fitness_cases);      \
    for(int i = 0; i < run.fitness_cases; ++i) {                                        \
        _CN_RESULT_VAR[i] = CODE;                                                       \
    }                                                                                   \
} while(0)


// Available GP functions (aside from the user defined)
#define compute_add(C1, C2) (C1[i] + C2[i])
#define compute_sub(C1, C2) (C1[i] - C2[i])
#define compute_div(C1, C2) ((C1[i] != 0.0 && C2[i] != 0.0) ? (C1[i] / C2[i]) : 0.0)
#define compute_mul(C1, C2) (C1[i] * C2[i])
#define compute_sin(C1) (sin(C1[i]))
#define compute_cos(C1) (cos(C1[i]))
#define compute_tan(C1) (tan(C1[i]))
#define compute_if(C1, C2, C3) ((C3[i] < 0) ? C1[i] : C2[i])


#define WAS_ALLOCATED(NODE) (((NODE)->primitive < VAR_START || (NODE)->primitive > VAR_END) && (!ENABLE_CACHE || (NODE)->candid <= -1))


#define CALL_COMPUTE_NODE_TEMPLATE(N_CHILDS, CHILD_ACC, COMPUTE_FUNC) do {  \
    fit_t *_CN_CHILDS_VAR[N_CHILDS];                                        \
    for(int i = 0; i < N_CHILDS; ++i) {                                     \
        _CN_CHILDS_VAR[i] = compute_node_cache(CHILD_ACC);                  \
    }                                                                       \
    if (PDEBUG) {                                                           \
        char *temp_str = get_dag_expr(t);                                   \
        printfd("allocating space for result of node %s\n", temp_str);      \
        free(temp_str);                                                     \
    }                                                                       \
    COMPUTE_NODE_TEMPLATE(COMPUTE_FUNC);                                    \
    for(int i = 0; i < N_CHILDS; ++i) {                                     \
        if (WAS_ALLOCATED(CHILD_ACC)) {                                     \
            free(_CN_CHILDS_VAR[i]);                                        \
        }                                                                   \
    }                                                                       \
}while(0)


fit_t *compute_node_cache(const dag_node *t) {
    const int is_cached = ENABLE_CACHE && t->candid > -1;
    
    //printfd("Entering fucntion! %s, mem: %d\n", primitive_set[t->primitive].name, allocated_mem);
    if (is_cached && cache.valid[t->candid]) {
        printfd("Retrieving from cache(%p), at pos(%d, %p), %d bytes\n", cache.data, t->candid, &cache.data[run.fitness_cases * t->candid], run.fitness_cases);
        cache.accesses++;
        return &cache.data[run.fitness_cases * t->candid];
    } else {

        fit_t *_CN_RESULT_VAR;
        printfd("Doing switch! %s\n", primitive_set[t->primitive].name);
        switch (t->primitive) {
            case X: default:
                _CN_RESULT_VAR = &cache.vars[0];
                break;
            case Y:
                _CN_RESULT_VAR = &cache.vars[run.fitness_cases];
                break;
            case Scalar:
                COMPUTE_NODE_TEMPLATE(t->children.f[0]);
                break;
            case Add: 
                CALL_COMPUTE_NODE_TEMPLATE(2, t->children.n[i], 
                    compute_add(_CN_CHILDS_VAR[0], _CN_CHILDS_VAR[1]));
                break;
            case Sub:
                CALL_COMPUTE_NODE_TEMPLATE(2, t->children.n[i], 
                    compute_sub(_CN_CHILDS_VAR[0], _CN_CHILDS_VAR[1]));
                //CALL_COMPUTE_NODE_TEMPLATE(2, t->children.n[i], compute_sub(childs[0], childs[1]));
                break;
            case Div:
                CALL_COMPUTE_NODE_TEMPLATE(2, t->children.n[i], 
                    compute_div(_CN_CHILDS_VAR[0], _CN_CHILDS_VAR[1]));
                //CALL_COMPUTE_NODE_TEMPLATE(2, t->children.n[i], compute_div(childs[0], childs[1]));
                break;
            case Mul:
                CALL_COMPUTE_NODE_TEMPLATE(2, t->children.n[i], 
                    compute_mul(_CN_CHILDS_VAR[0], _CN_CHILDS_VAR[1]));
                //CALL_COMPUTE_NODE_TEMPLATE(2, t->children.n[i], compute_mul(childs[0], childs[1]));
                break;
            case Cos:
                CALL_COMPUTE_NODE_TEMPLATE(1, t->children.n[i], 
                    compute_cos(_CN_CHILDS_VAR[0]));
                //CALL_COMPUTE_NODE_TEMPLATE(1, t->children.n[i], compute_cos(childs[0]));
                break;
            case Sin:
                CALL_COMPUTE_NODE_TEMPLATE(1, t->children.n[i], 
                    compute_sin(_CN_CHILDS_VAR[0]));
                //CALL_COMPUTE_NODE_TEMPLATE(1, t->children.n[i], compute_sin(childs[0]));
                break;
            case Tan:
                CALL_COMPUTE_NODE_TEMPLATE(1, t->children.n[i], 
                    compute_tan(_CN_CHILDS_VAR[0]));
                //CALL_COMPUTE_NODE_TEMPLATE(1, t->children.n[i], compute_tan(childs[0]));
                break;
            case If:
                CALL_COMPUTE_NODE_TEMPLATE(3, t->children.n[i], 
                    compute_if(_CN_CHILDS_VAR[0], _CN_CHILDS_VAR[1], _CN_CHILDS_VAR[2]));
                //CALL_COMPUTE_NODE_TEMPLATE(3, t->children.n[i], compute_if(childs[0], childs[1], childs[2]));
                break;
        }
        if (is_cached && !cache.valid[t->candid]) {
            //fit_t *ptr_to_cache = &cache.data[t->candid * run.fitness_cases];
            memcpy(&cache.data[t->candid * run.fitness_cases], _CN_RESULT_VAR, sizeof(fit_t) * run.fitness_cases);
            cache.valid[t->candid] = 1;

            if (PDEBUG) {
                char *temp_str = get_dag_expr(t);
                printfd("caching (and freeing) node %s (pos %d)\n", temp_str, t->candid);
                free(temp_str);
            }        
            
            free(_CN_RESULT_VAR);
            cache.validated++;
            return &cache.data[t->candid * run.fitness_cases];
        }
        return _CN_RESULT_VAR;
    }
}


fit_t compute_node(const dag_node *t) {
    int is_cached = ENABLE_CACHE && t->candid > -1;
    
    if (is_cached && cache.valid[t->candid]) {
        if (run.index == 0) {
            cache.accesses++;
            /*
            printf("Found a cached node: \n");
            print_dag_node(t);
            fancy_dag_print(t);
            printf("The cache state is :\n");
            print_cache(1);
            */
        }
        return cache.data[t->candid * run.fitness_cases + run.index];
    } else {

        fit_t result;
        switch(t->primitive) {
            case X: default:
                result = run.cur_vars[0];
                break;
            case Y:
                result = run.cur_vars[1];
                break;
            case Scalar:
                result = t->children.f[0];
                break;
            case Add: 
                result = compute_node(t->children.n[0]) + compute_node(t->children.n[1]);
                break;
            case Sub: 
                result = compute_node(t->children.n[0]) - compute_node(t->children.n[1]);
                break;
            case Div: {
                fit_t c1 = compute_node(t->children.n[0]);
                fit_t c2 = compute_node(t->children.n[1]);
                result = (c1 != 0.0 && c2 != 0.0) ? (c1 / c2) : 0.0;
                break;
            }
            case Mul:
                result = compute_node(t->children.n[0]) * compute_node(t->children.n[1]);
                break;
            case Cos:
                result = cos(compute_node(t->children.n[0]));
                break;
            case Sin:
                result = sin(compute_node(t->children.n[0]));
                break;
            case Tan:
                result = tan(compute_node(t->children.n[0]));
                break;
            case If: {
                fit_t c1 = compute_node(t->children.n[0]);
                fit_t c2 = compute_node(t->children.n[1]);
                fit_t c3 = compute_node(t->children.n[2]);
                result = (c3 < 0) ? c1 : c2;      
                }
                break;

        }
        if (is_cached && !cache.valid[t->candid]) {
            cache.data[t->candid * run.fitness_cases + run.index] = result;
            if (run.index == run.fitness_cases - 1) {
                cache.valid[t->candid] = 1;
                cache.validated++;
            }
        }
        return result;
    }
}


fit_t compute_node_debug(const dag_node *t) {
    fit_t result = compute_node(t);
    if (isinf(result) || isnan(result) || !isfinite(result)) {
        //printf("Bad operation (%s) between values %.3f and %.3f\n", primitive_set[t->primitive].name, c1, c2);
        printf("Bad operation (%s) between values: result: %.3f", primitive_set[t->primitive].name, result);
    }
    return result;
}


int calc_pop_fit(tree **population) {
    float best_fit = FLT_MAX;
    int best_ind = -1;
    for(int i = 0; i < run.pop_size; ++i) {
        //print_tree(population[i]);
        printfd("Calculating tree i: %d\n", i);
        float cur_fit = calc_tree_fit(population[i]);
        //float cur_fit = no_tree_fit(population[i]);
        if(cur_fit < best_fit) {
            best_fit = cur_fit;
            best_ind = i;
        }
    }
    if (best_ind <= -1) {
        printf(PREERR"Bad best index found while calculating fitness: %d\n", best_ind);
    }
    return best_ind;
}


fit_t no_tree_fit(tree *t) {
    t->fitness = 0;
    return 0;
}


fit_t calc_tree_fit(tree *t) {
    //print_tree(t, 0, 0);

    if (TENSOR_CACHE) {
        fit_t *domain = compute_node_cache(t->dag);
        //print_domain(domain, fitness_cases, resolution); //debug
        t->fitness = compute_fitness(domain);
        if (WAS_ALLOCATED(t->dag)) {
            free(domain);
        }
    } else {

        fit_t *domain = compute_domain(t->dag);
        //print_domain(domain, fitness_cases, resolution); //debug
        t->fitness = compute_fitness(domain);
        free(domain);
    }

    return t->fitness;
}


fit_t *compute_domain(dag_node *t) {
    int i;
    fit_t *domain = (fit_t *)malloc(sizeof(*domain) * run.fitness_cases);
    reset_cur_vars();

    for(i = 0; i < run.fitness_cases; ++i) {
        domain[i] = compute_node_debug(t);
        //domain[i] = 0;

        STEP_DOMAIN(run);
        //printf("Updating index, %d\n", run.index);
    }
    return domain;
}


fit_t compute_fitness(fit_t *values) {
    fit_t sq_sum = 0;
    for(int i = 0; i < run.fitness_cases; ++i) {
        sq_sum += (run.target[i] - values[i]) * (run.target[i] - values[i]);
    }
    return sqrt(min(sq_sum, FLT_MAX) / (fit_t)(run.fitness_cases));
}