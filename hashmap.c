#include <stdio.h>
#include <stdlib.h>
#include "gp.h"


HashTable *create_hashtable(uint32_t n) {
    HashTable *htable = (HashTable *)calloc(1, sizeof(HashTable));
    htable->size = n;
    htable->table = (b_node *)calloc(n + (OA_NODES - 1), sizeof(b_node));
    for (int i = VAR_START; i <= VAR_END; ++i) {
        htable->table[i - VAR_START].node.primitive = i;
        htable->table[i - VAR_START].node.candid = -1;
    }
    
    htable->candidates = (cache_node *)calloc(INIT_CAND_SIZE, sizeof(cache_node));
    
    htable->cand_size = INIT_CAND_SIZE;
    htable->n_nodes = (VAR_END - VAR_START) + 1;
    htable->child_ptrs = (dag_node ***)malloc(N_CHILDREN_ARR * sizeof(dag_node **));
    htable->scalar_ptrs = (float **)malloc(N_CHILDREN_ARR * sizeof(float *));
    htable->child_ptrs[0] = (dag_node **)malloc(CHILD_INIT_SIZE * sizeof(dag_node *));
    htable->scalar_ptrs[0] = (float *)malloc(SCALAR_INIT_SIZE * sizeof(float));
    htable->cur_child_arr_size = CHILD_INIT_SIZE;
    htable->cur_scalar_arr_size = SCALAR_INIT_SIZE;

    // +OA_NODES safeguards
    for(int i = 0; i < OA_NODES - 1; ++i) {
        htable->table[n + i].node.primitive = 0xFFFF; // not empty, 0xFFFF should to the safeguards 
        htable->table[n + i].node.arity = 0;
        htable->table[n + i].node.children.f = htable->scalar_ptrs[0]; // make something not NULL to be accepted to memcmp 
    }
    return htable;
}


void print_table_stats(b_node *table, int n, int do_extra_stats) {
    int sum_slots = 0;
    int max_slots = INT_MIN;
    int min_slots = INT_MAX;
    float *data_slots = (float *)malloc(n * sizeof(float));
    int colisions = 0;
    int unused_slots = 0;
    int buckets_with_n[10] = {0};

    for(int i = 0; i < n; ++i) {
        b_node *cur = &table[i];
        int slots = 0;
        if (cur->node.primitive == 0) {
            unused_slots++;
            slots--;
        }

        while(cur != NULL) {
            ++slots;   
            cur = cur->next;
        }

        if (slots > 1) {
            colisions += slots - 1;
        }

        if(do_extra_stats && slots < 10) {
            buckets_with_n[slots]++;
        }

        data_slots[i] = (float)slots;
        sum_slots += slots;
        min_slots = min(min_slots, slots);
        max_slots = max(max_slots, slots);
    }
    float mean_slots = (float)(sum_slots) / n;
    float std_slots = compute_std_float(data_slots, n, sum_slots);

    int alloced_mem = (sizeof(HashTable) + (n * sizeof(b_node))) + (sizeof(b_node) * sum_slots);
    // print statistics
    printf("\n======================================\n");
    printf("Statistics of slots in each bucket: \n");
    printf("Size of table: %d\n", n);
    printf("%% of colisions: %.3f %%\n", ((double)colisions / (double)sum_slots) * 100.0);    
    printf("Colisions: %d\n", colisions);
    printf("Nodes: %d\n", sum_slots);
    printf("Mean: %.3f\n", mean_slots);
    printf("Std: %.3f\n", std_slots);
    printf("Min: %d\n", min_slots);
    printf("Max: %d\n", max_slots);
    printf("Unused buckets/slots: %d (%d bytes)\n", unused_slots, unused_slots << 5);
    printf("Memory allocated: %d bytes.\n", alloced_mem);
    printf("%% of unused memory: %.3f %%\n", (float)(unused_slots << 5) / (float)(alloced_mem) * 100.0);
    if (do_extra_stats) {
        printf("Bucket stats:\n\n");
        for(int i = 0; i < 10; ++i) {
            printf("%d bucks with %d els\n", buckets_with_n[i], i);
        }
    }
    printf("======================================\n\n");

    free(data_slots);
}


void free_scalar_ptrs(HashTable *t) {
    FREE_PTRS(t, scalar_ptrs, scalar_arr);
}


void free_child_ptrs(HashTable *t) {
    FREE_PTRS(t, child_ptrs, child_arr);
}


void free_hashtable(HashTable *t) {
    delete_dag_table(t->table, t->size + (OA_NODES - 1));
    free_scalar_ptrs(t);
    free_child_ptrs(t);
    free(t->candidates);
    free(t);
}


void print_dag_table(b_node *table, int n) {
    printf("\nHashtable: \n");
    for(int i = 0; i < n + (OA_NODES - 1); ++i) {
        b_node *cur = &table[i];
        printf("Bucket %d:\n", i);
        while(cur != NULL && cur->node.primitive != 0) {
            print_dag_node(&cur->node);
            cur = cur->next;
        }
    }
}


void delete_dag_table(b_node *table, uint32_t tab_size) {
    for(int i = 0; i < tab_size; ++i) {
        b_node *cur = table[i].next;

        while(cur != NULL) {
            b_node *temp = cur->next;
            free(cur);
            cur = temp;
        }
    }
    //memset(table, 0, sizeof(b_node) * tab_size);
    free(table);
}


uint32_t hash_dag_node(uint32_t primitive, int arity, c_node children) {
    uint32_t hash_code = primitive * 0x9e3779b9;
    
    if (primitive != Scalar) {
        for(int i = 0; i < arity; ++i) {
            hash_code ^= (GET_CHILD(children.n, i)) + (hash_code >> 6);
        }
    } else {
        for(int i = 0; i < arity; ++i) {
            hash_code ^= (GET_SCALAR_INT(children.f, i)) + (hash_code >> 6);
        }
    }
    
    return hash_code;
}


dag_node **alloc_child_ptrs(HashTable *t, int slots) {
    ALLOC_PTRS(t, child_ptrs, child_ind, child_arr, cur_child_arr_size, dag_node*, "Child");
}


float *alloc_scalar_ptrs(HashTable *t, int slots) {
    ALLOC_PTRS(t, scalar_ptrs, scalar_ind, scalar_arr, cur_scalar_arr_size, float, "Scalar");
}


void print_candidate_list(HashTable *t, int print_cache_cands) {
    printf("\nCandidate list:\n");
    int to_print = print_cache_cands ? t->cand_index : t->cand_size;
    printf("Size of Candidate List: %d\n", t->cand_size);
    printf("Candidate List index: %d\n", t->cand_index);
    printf("\nIndex:\t(DAG ptr,\tCache Index,\tImportance)\n\n");
    for(int i = 0; i < to_print; ++i) {
        printf("%d:\t(%p,\t%d,\t%d)\n", i, t->candidates[i].dag, t->candidates[i].index, t->candidates[i].importance);
    }
    printf("\n");
}


dag_node *add_dag_node(HashTable *t, uint32_t primitive, c_node children) {
    return add_dag_node_index(t, primitive, children, -1);
}


dag_node *add_dag_node_index(HashTable *t, uint32_t primitive, c_node children, int ind) {
    dag_node *dag;

    // if not a variable of the terminal set
    if (primitive > VAR_END || primitive < VAR_START) {
        int arity = primitive_set[primitive].arity;
        int index = hash_dag_node(primitive, arity, children) % t->size;

        b_node *cur_node = t->table + index;
        dag = &cur_node->node;
        int i, is_empty;

        // search OA_NODES in liner probing fashion
        for (i = 0; i < OA_NODES && !((is_empty = (cur_node->node.primitive == 0)) || DAG_EQ(dag, arity, primitive, children)); ++i) {
            dag = &(++cur_node)->node;
        }

        // if neither the node or an empty slot was found, search by chaining
        if (i == OA_NODES) {
            b_node *last_node = cur_node - OA_NODES;
            cur_node = last_node->next;
            while(cur_node != NULL) {
                dag = &cur_node->node; // remove dag_node *
                if (DAG_EQ(dag, arity, primitive, children)) {
                    break;
                }
                last_node = cur_node;
                cur_node = cur_node->next;
            }
            if (cur_node == NULL) {
                last_node->next = (b_node *)malloc(sizeof(b_node));
                last_node->next->next = NULL;
                dag = &last_node->next->node;
            }
        }

        // node not found, add it to the hashtable
        if (is_empty || cur_node == NULL) {
            t->n_nodes++;

            int nodes = IS_FUNC(primitive, arity);
            if (nodes) { // if it is a primitive
                for (int i = 0; i < arity; ++i) {
                    nodes += children.n[i]->n_offspring;
                }
            }

            dag->children = children;
            dag->primitive = primitive;
            dag->arity = arity;
            dag->frequency = 1;
            dag->n_offspring = nodes;
            dag->candid = -1;
            
        } else { // node found, free "children" arg passed (we dont need anymore), handle cache 
            dag->frequency++;
            if (primitive != Scalar) { // free up pointer for childs
                t->child_ind = max(0, t->child_ind - arity);
            } else {
                t->scalar_ind = max(0, t->scalar_ind - arity);
            }
            
            if (ENABLE_CACHE) {
                handle_candidates(dag, t, ind);
            }
        }

    } else { // otimization for terminal set variables (x, y, z...)
        dag = &t->table[primitive - VAR_START].node;
        dag->frequency++;
    }

    return dag;
}