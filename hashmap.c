#include <stdio.h>
#include <stdlib.h>
#include "gp.h"


HashTable *create_hashtable(uint32_t n) {
    HashTable *htable = (HashTable *)calloc(1, sizeof(HashTable));
    htable->size = n;
    //htable->table = (b_node *)calloc(n, sizeof(b_node));
    htable->table = (b_node *)calloc(n + (OA_NODES - 1), sizeof(b_node));
    for (int i = VAR_START; i <= VAR_END; ++i) {
        htable->table[i - VAR_START].node.primitive = i;
        htable->table[i - VAR_START].node.candid = -1;
    }
    
    //htable->candidates = (cache_node *)calloc(INIT_CAND_SIZE, sizeof(cache_node));
    htable->candidates = NULL;
    
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
    for(int i = 0; i < n + (OA_NODES - 1); ++i) {
        b_node *cur = &table[i];
        printf("Bucket :%d %d\n", i, cur->node.primitive);
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
                //handle_candidates(dag, t, ind);
            }
        }

    } else { // otimization for terminal set variables (x, y, z...)
        dag = &t->table[primitive - VAR_START].node;
        dag->frequency++;
    }

    return dag;
}