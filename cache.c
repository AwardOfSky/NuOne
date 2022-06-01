#include <stdio.h>
#include <stdlib.h>
#include "gp.h"


void handle_candidates(dag_node *dag, HashTable *t, int ind) { 
    int node_imp = IMPORTANCE(dag->frequency, dag->n_offspring);

    if (node_imp >= MIN_CACHE_REQ) { // meets min cache req
        cache_node *get_dag_in_cand = NULL;
        if (dag->candid == -1) { // not in the cadidate list, let's add

            if (t->cand_index >= t->cand_size) {
                if (t->cand_index >= MAX_CACHE_CANDIDATES) {
                    //printf(PREWARN"Candidate list full.\n");
                    return;
                }

                //t->candidates = realloc_cache_node(t->candidates , t->cand_size * FACTOR_TO_GROW_CANDLIST, (uint32_t *)(&(t->cand_size)));
                realloc_cache_node(&(t->candidates), t->cand_size * FACTOR_TO_GROW_CANDLIST, (uint32_t *)(&(t->cand_size)));
            }
            //if (t->cand_index >= 4) printf("new cand node! %d %d\n", t->cand_index, ind);

            dag->candid = -t->cand_index++ - 2;
            get_dag_in_cand = &(t->candidates[-dag->candid - 2]);
            get_dag_in_cand->dag = dag; // get ref to this node

            // -1 means that this node is new (not in the old hashtable)
            // if node is not new and its cache index is valid
            // (was actually calculated in last gen, then add it, no need to recalculate)
            if (ind > -1 && cache.valid[ind]) {
                get_dag_in_cand->index = ind;
            } else {
                get_dag_in_cand->index = -1;
            }
        } else if (dag->candid < -1) { // already in the candidate list
            //if (t->cand_index >= 4) printf("already in the cand_list! %d %d\n", t->cand_index, ind);
            get_dag_in_cand = &(t->candidates[-dag->candid - 2]);
        } else { // positive candid means index in cache already
            printf(PREWARN"This Dag node is already in the cache (index %d). Maybe the cache was already build?\n", dag->candid);
        }

        // update importance even if already in cache to avoid having to updatelater
        // assumes # of repeats adds that are already candidates are lower
        // than the # of candidates, which might not be true at all
        // ( VICandidates will be added lots of times...)
        get_dag_in_cand->importance = node_imp;
    }
}


void init_cache(size_t n) {
    cache.vars = (fit_t *)malloc(DIMS * run.fitness_cases * sizeof(*cache.vars));
    if (ENABLE_CACHE) {
        cache.size = n;
        cache.data = (fit_t *)malloc(n * run.fitness_cases * sizeof(*(cache.data)));
        cache.valid = (char *)calloc(n, sizeof(*cache.valid));
        cache.occupied = 0;
    }
}


void free_cache() {
    free(cache.vars);
    if (ENABLE_CACHE) { 
        free(cache.data);
        free(cache.valid);
    }
}


void realloc_cache(uint32_t new_size) {
    realloc_char(&(cache.valid), new_size, (uint32_t *)(&(cache.size)));
    realloc_fit_t(&(cache.data), new_size * run.fitness_cases, NULL);
}


void print_cache(uint32_t *format, int print_invalid) {
    printf("\nCache state:\n");
    printf("Cache size: %d\n", cache.size);
    printf("Occupied: %d (%.3f%%)\n", cache.occupied, (((double)cache.occupied / (double)cache.size)) * 100.0);
    printf("Cached accessed %d times.\n", cache.accesses);
    printf("Cache nodes validated: %d\n", cache.validated);
    
    printf("Cache Table\t: (validity,\tdata)\n");
    for (int i = 0; i < cache.size; ++i) {
        printf("\nIndex %d\t\t: %d\n", i, cache.valid[i]);
        if (cache.valid[i] || print_invalid) {
            printf("Data:");
            print_domain(&cache.data[i * run.fitness_cases], run.fitness_cases, format);
        }
    }

    printf("\n");
}


int cmp_candidates(const void *a, const void *b) {
    int aval = ((cache_node *)a)->importance;
    int bval = ((cache_node *)b)->importance;
    int res = bval - aval;
    return res;
}

// int *const => const ptr to int
void subtract_cand_imp(HashTable *t, dag_node *dag, int dep_lim, int *imp, int *freq) {
    for (int i = 0; i < dag->arity; ++i) {
        dag_node *child = dag->children.n[i];
        if (child->candid < -1) {
            int am = t->candidates[-(child->candid) - 2].dag->n_offspring * (*freq);
            *imp -= am;
        } else if ((dep_lim != 1) && IS_FUNC(child->primitive, child->arity)) {
            subtract_cand_imp(t, child, dep_lim - 1, imp, freq);
        }
    }
}


// dep_lim = [1, +inf]. 0, -1, ... to disable depth limit
int get_cand_imp_to_sub(HashTable *t, dag_node *dag, int32_t dep_lim) {
    int res = 0;
    for(int i = 0; i < dag->arity; ++i) {
        dag_node *child = dag->children.n[i];
        if(child->candid < -1) {
            res += t->candidates[-child->candid - 2].dag->n_offspring;
        } else if (dep_lim != -1 && IS_FUNC(child->primitive, child->arity)) {
            res += get_cand_imp_to_sub(t, dag, dep_lim - 1);
        }
    }
    return res;
}


#define DEB_UPDATE 0

void build_cache(HashTable *t, int gen) {

    cache.validated = 0;
    cache.accesses = 0;

    printfd("Cache size: %d\n", cache.size); 
    // 0 - check if we need to grow the cache
    int shrink = 1;
    if ((t->cand_index >= cache.size) && !((gen + 1) % GROW_CACHE_N_GEN)) {
        
        uint32_t tentative_size = max(next_power_2(t->cand_index), cache.size * FACTOR_TO_GROW_CACHE);
        if (MAX_CACHE_SIZE != -1) {
            tentative_size = min(tentative_size, MAX_CACHE_SIZE);
        }
        if (MAX_CACHE_MEMORY != -1) {
            tentative_size = min(tentative_size, MAX_CACHE_MEMORY / (run.fitness_cases * sizeof(fit_t)));
        }
        if (tentative_size > cache.size) {
            shrink = 0;
            realloc_cache(tentative_size);
        } else {
            printf("\n"PREWARN"Cache is at 100%% capacity, cannot grow.\n");
            printf("Current cache (size, memory): (%d, %ld), out of (%d, %d)\n", 
                cache.size, cache.size * run.fitness_cases * sizeof(fit_t),
                MAX_CACHE_SIZE, MAX_CACHE_MEMORY);
            printf("Tentative cache (size, memory): (%d, %ld)\n\n", tentative_size, tentative_size * run.fitness_cases * sizeof(fit_t));
        }

        printfd("Grow: Cache size, cand index: %d %d\n", cache.size, t->cand_index);
    }


    // 0.5 - update importance list
    int i = 0;
    for(i = 0; i < t->cand_index; ++i) {
        dag_node *dag = t->candidates[i].dag;
        //printf("Updating importance of dag of index %d: %s\n", i, get_dag_expr(dag));
        if (IS_FUNC(dag->primitive, dag->arity)) {
            //subtract_cand_imp(t, dag, -1, &(t->candidates[i].importance), &(dag->frequency));
            t->candidates[i].importance -= get_cand_imp_to_sub(t, dag, -1) * dag->frequency;
        }
    }


    // 1 - sort candidates
    if (t->cand_index < 128) {
        qsort(t->candidates, t->cand_index, sizeof(*t->candidates), cmp_candidates);
    } else {
        RADIX_SORT_STRUCT(t->candidates, t->cand_index, 0, cache_node, int, .importance);
    }
    //for(int i = 0; i < t->cand_index; ++i) printf("%d\n", t->candidates[i].importance);

    // 2 - invalidate cache
    memset(cache.valid, 0, cache.size);
    int *voc;
    int voci;

    //check if we should (and can) shrink the cache
    if (1 && shrink && (t->cand_index * FACTOR_TO_SHRINK_CACHE < cache.size) && ((cache.size >> 1) > MIN_CACHE_SIZE) && !((gen + 1) % SHRINK_CACHE_N_GEN)) {
        printfd("Shrinking cache!\n");
        voc = (int *)malloc(sizeof(int) * t->cand_index); // arr of indices of voc (valid old candidates)
        voci = 0;
    } else {
        shrink = 0;
    }

    // 3 - set old candidates up until n_candidates to valid  
    int n_candidates = min(t->cand_index, cache.size);
    if (DEB_UPDATE) printf("\nCache size: %d, n candidates: %d\n", cache.size, n_candidates);
    cache.occupied = n_candidates;
    for (i = 0; i < n_candidates; ++i) {
        int index = t->candidates[i].index;
        if (index != -1) { // not -1 means its old dag
            if (DEB_UPDATE) printf("Found old dag, inserting at cache index: %d\n", index);
            cache.valid[index] = 1;
            t->candidates[i].dag->candid = index;
            //t->candidates[i].index = index;
            if (shrink) {
                voc[voci++] = i;
            }
        }
    }

    // 4 - assign empty spaces in cache for new dags
    int index_to_assign = 0;
    for(i = 0; i < n_candidates; ++i) {
        if (t->candidates[i].index == -1) { // new dag node
            while (cache.valid[index_to_assign]) {
                ++index_to_assign;
            }
            if (DEB_UPDATE) printf("Adding new dag, inserting at cache index: %d\n", index_to_assign);
            t->candidates[i].dag->candid = index_to_assign;
            t->candidates[i].index = index_to_assign++;
        }
    }


    // 5 - memcpy remaining vocs that are in higher indeces in order to shrink the cache
    if (shrink) {
        int cur_ind = index_to_assign;
        if (DEB_UPDATE) printf("Index to aasign %d, n candidates: %d\n", index_to_assign, n_candidates);
        for (i = 0; i < voci; ++i) {
            int voc_cache_index = t->candidates[voc[i]].index;
            if (DEB_UPDATE) printf("voc_cache index: %d, voc[%d]: %d\n", voc_cache_index, i, voc[i]);
            if (voc_cache_index >= (cache.size >> 1)) { // NOTE: this is not index to assing, right? It should be the size to shrink (next_power_2(t->cand_index))

                while (cache.valid[cur_ind]) {
                    ++cur_ind;
                }

                if (cur_ind >= n_candidates) {
                    printf(PREERR"Cached too many candidates: ");
                    printf("while inserting candidate at candidate list index %d, cache index %d, cur_index: %d\n", voc[i], voc_cache_index, cur_ind);
                }
                cache.valid[voc_cache_index] = 0;
                cache.valid[cur_ind] = 1;
                if (DEB_UPDATE) printf("Invalidating cache index %d and validating index %d\n", voc_cache_index, cur_ind);
                t->candidates[voc[i]].dag->candid = cur_ind;
                t->candidates[voc[i]].index = cur_ind;
                memcpy(&cache.data[run.fitness_cases * cur_ind], &cache.data[run.fitness_cases * voc_cache_index], run.fitness_cases * sizeof(*cache.data));
                ++cur_ind;
            }
        }
        free(voc);
        
        // 5.5 - actually shrink the cache
        realloc_cache(cache.size >> 1);
        printfd("Shrink: Cache size, cand index: %d %d\n", cache.size, t->cand_index);
    }
}


MAKE_REALLOC(cache_node);
MAKE_REALLOC(fit_t);

// fit_t *realloc_fit_t(fit_t *data, const uint32_t new_size, uint32_t *sptr) {
//     fit_t *tp = (fit_t*)realloc(data, sizeof(*tp) * new_size);
//     if (tp != NULL || new_size == 0) { 
//         data = tp; 
//         if (sptr != NULL) { 
//             *sptr = new_size;
//         }
//     } else { 
//         printf(PREERR"OOM: Failed to realloc array of type ""fit_t"", on memory address %p. ", data);
//         printf("Requested %lu bytes.\n", sizeof(*data) * new_size);
//     }
//     return data;
// };
