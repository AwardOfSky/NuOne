#include <stdio.h>
#include <stdlib.h>
#include "gp.h"

void handle_candidates1(dag_node *dag, HashTable *t, int ind) { 
    //printf("Entering cand handling\n");
    int node_imp = IMPORTANCE(dag->frequency, dag->n_offspring);

    if (node_imp >= MIN_CACHE_REQ) { // meets min cache req
        cache_node *get_dag_in_cand = NULL;
        //printf("dag candid: %d, ind: %d, cand_index: %u\n", dag->candid, ind, t->cand_index);
        if (dag->candid == -1) { // not in the cadidate list, let's add

            if (t->cand_index >= t->cand_size) {
                //realloc_candidate_list(t, t->cand_size * FACTOR_TO_GROW_CANDLIST);
                realloc_cache_node(t->candidates, t->cand_size * FACTOR_TO_GROW_CANDLIST, (size_t *)(&(t->cand_size)));
            }

            dag->candid = -t->cand_index++ - 2;
            get_dag_in_cand = &(t->candidates[-dag->candid - 2]);
            get_dag_in_cand->dag = dag; // get ref to this node

            // -1 means that this node is new (not in the old hashtable)
            // if node is not new and its cache index is valid
            // (was actually calculated in last gen, then add it, no need to recalculate)
            if (ind > -1 && cache->valid[ind]) {
                get_dag_in_cand->index = ind;
            } else {
                get_dag_in_cand->index = -1;
            }
        } else if (dag->candid < -1) { // already in the candidate list
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
    //printf("Leaving cand handling\n");
}



void handle_candidates(dag_node *dag, HashTable *t, int ind) { 
    int node_imp = IMPORTANCE(dag->frequency, dag->n_offspring);

    // meets minimum requirement to be a candidate for cache
    if (node_imp >= MIN_CACHE_REQ) {
        cache_node *get_dag_in_cand = NULL;

        // not in the cadidate list, let's add
        if (dag->candid == -1) {
            if (t->cand_index >= t->cand_size) {
                t->candidates = realloc_cache_node(t->candidates, t->cand_size * FACTOR_TO_GROW_CANDLIST, (size_t *)(&(t->cand_size)));
            }

            dag->candid = t->cand_index++;
            get_dag_in_cand = &(t->candidates[dag->candid]);
            get_dag_in_cand->dag = dag; // get ref to this node

            if (ind > -1 && cache->valid[ind]) {
                get_dag_in_cand->index = ind;
            } else {
                get_dag_in_cand->index = -1;
            }
        
        } else { // already in the candidate list
            get_dag_in_cand = &(t->candidates[dag->candid]);
        }

        get_dag_in_cand->importance = node_imp;
    }
}


Cache *create_cache(size_t fit_cases, size_t n) {
    if (!ENABLE_CACHE) {
        return NULL;
    }

    Cache *res = (Cache *)malloc(sizeof(*res));
    res->size = n;
    res->data = (fit_t *)malloc(n * fit_cases * sizeof(*res->data));
    res->vars = (fit_t *)malloc(DIMS * fit_cases * sizeof(*res->vars));
    res->valid = (char *)calloc(n, sizeof(*res->valid));
    res->occupied = 0;
    return res;
}


// TODO: if cache is global then maybe we shoould remove the argument to function
void free_cache(Cache *c) {
    if (ENABLE_CACHE) { 
        free(c->data);
        free(c->valid);
        free(c->vars);
        free(c);
    }
}


void realloc_cache(uint32_t new_size, size_t fit_cases) {
    cache->valid = realloc_char(cache->valid, new_size, (size_t *)(&(cache->size)));
    cache->data = realloc_fit_t(cache->data, new_size * fit_cases, NULL);
}


void print_cache(uint32_t *format, int print_invalid) {
    uint64_t fitness_cases = 1;
    for (int i = 0; i < DIMS; ++i) {
        fitness_cases *= format[i];
    }
    printf("\nCache state:\n");
    printf("Cache size: %d\n", cache->size);
    printf("Occupied: %d (%.3f%%)\n", cache->occupied, (((double)cache->occupied / (double)cache->size)) * 100.0);
    printf("Cached accessed %d times.\n", cache->accesses);
    printf("Cache nodes validated: %d\n", cache->validated);
    
    printf("Cache Table\t: (validity,\tdata)\n");
    for (int i = 0; i < cache->size; ++i) {
        printf("\nIndex %d\t\t: %d\n", i, cache->valid[i]);
        if (cache->valid[i] || print_invalid) {
            printf("Data:");
            print_domain(&cache->data[i * fitness_cases], fitness_cases, format);
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


void build_cache(HashTable *t, int gen, int fitness_cases) {

    cache->validated = 0;
    cache->accesses = 0;

    // 0 - check if we need to grow the cache
    int shrink = 1;
    printf("before Cache index, size: %d, %d\n", t->cand_index, cache->size);
    if ((t->cand_index * FACTOR_TO_GROW_CACHE > cache->size) && !((gen + 1) % UPDATE_CACHE_N_GEN)) {
        uint32_t tentative_size = cache->size << 1; // it can be this for now
        //printfd("@@@@@@@@@@@Growing cache! %d\n", tentative_size);

        if (tentative_size <= MAX_CACHE_SIZE) {
            shrink = 0;
            realloc_cache(tentative_size, fitness_cases);
        } else {
            printf(PREWARN"Cache is at 100%% capacity, cannot grow.\n");
        }
    }
    printf("after Cache index, size: %d, %d\n", t->cand_index, cache->size);

    //if (t->cand_index > RADIX_SORT_THRESH) {
    //    radix_strut(t->candidates, t->cand_index);
    //} else {
        //radix_strut(t->candidates, t->cand_index);
    //}
    // 1 - sort candidates
    qsort(t->candidates, t->cand_index, sizeof(*t->candidates), cmp_candidates);

    // 2 - invalidate cache
    memset(cache->valid, 0, cache->size);
    
    int *voc;
    int voci;

    //check if we should (and can) shrink the cache
    if (shrink && (t->cand_index * FACTOR_TO_SHRINK_CACHE < cache->size) && ((cache->size >> 1) > MIN_CACHE_SIZE) && !((gen + 1) % UPDATE_CACHE_N_GEN)) {
        printfd("Shrinking cache!\n");
        voc = (int *)malloc(sizeof(int) * t->cand_index); // arr of indices of voc (valid old candidates)
        voci = 0;
    } else {
        shrink = 0;
    }

    // 3 - set old candidates up until n_candidates to valid  
    int i;
    int n_candidates = min(t->cand_index, cache->size);
    //printfd("N candidates: %d\n", n_candidates);
    cache->occupied = n_candidates;
    for (i = 0; i < n_candidates; ++i) {
        int index = t->candidates[i].index;
        if (index != -1) { // not -1 means its old dag
            //printf("We should be passing here for (i, index): (%d, %d)\n", i, index);
            cache->valid[index] = 1;
            t->candidates[i].dag->candid = index;
            t->candidates[i].index = index;
            if (shrink) {
                voc[voci++] = i;
            }
        }
    }
    //printf("Ok here1!\n");
    //printf("Cand index: %d\n", t->cand_size);

    // 4 - assing empty spaces in cache for new dags
    int index_to_assign = 0;
    for(i = 0; i < n_candidates; ++i) {
        if (t->candidates[i].index == -1) { // new dag node
            //printfd("We should be here! %d %d\n", i, index_to_assign);
            while (cache->valid[index_to_assign]) {
                ++index_to_assign;
            }
            t->candidates[i].dag->candid = index_to_assign;
            t->candidates[i].index = index_to_assign++;
        }
    }

    // 5 - memcpy remaining vocs that are in higher indeces in order to shrink the cache
    if (shrink) {
        int cur_ind = index_to_assign;
        for (i = 0; i < voci; ++i) {
            int voc_cache_index = t->candidates[voc[i]].index;
            if (voc_cache_index > index_to_assign) { // NOTE: this is not index to assing, right? It should be the size to shrink (next_power_2(t->cand_index))
                cache->valid[voc_cache_index] = 0;
                cache->valid[cur_ind] = 1;
                t->candidates[voc[i]].dag->candid = cur_ind;
                t->candidates[voc[i]].index = cur_ind;
                memcpy(&cache->data[fitness_cases * cur_ind], &cache->data[fitness_cases * voc_cache_index], fitness_cases * sizeof(*cache->data));
                ++cur_ind;
            }
        }
        free(voc);

        if (cur_ind != n_candidates) {
            printf(PREERR"Incompatible indices while shrinking cache - current index: %d,\tn_candidates: %d\n", cur_ind, n_candidates);
        }
        
        // 5.5 - actually shrink the cache
        realloc_cache(cache->size >> 1, fitness_cases);
    }
    printf("after Cache index, size: %d, %d\n", t->cand_index, cache->size);
}


MAKE_REALLOC(cache_node);