//gcc -o engine.exe -Wall -Ofast engine.c gp.c hashmap.c utils.c genetics.c
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "gp.h"

#define STRESS(N, FUNC_CALL) do {                                   \
    clock_t start = clock();                                        \
    for(int i = 0; i < N; ++i) {                                    \
        FUNC_CALL;                                                  \
        printf("%u %d\n\n", seed, i);                               \
    }                                                               \
    double time_spent = (double)(clock() - start) / CLOCKS_PER_SEC; \
    printf("Time spent: %.3fs\n", time_spent);                      \
} while(0)


void test_tree_creation(uint32_t table_size, int stress);
void test_tree_generation(uint32_t table_size, int method, int mind, int maxd);
void test_subtree_crossover(uint32_t table_size, int stress);
void test_subtree_mutation(uint32_t table_size, int stress);


int main() {

    unsigned int seed = 1644863329; // reproducibility
    //unsigned int seed = time(NULL); // reproducibility

    srand(seed);
    printf("Random seed: %u\n", seed);

    
    //test_tree_generation(16, Full, 0, 0);
    
    //test_subtree_crossover(16, 0);
    //STRESS(10000, test_subtree_crossover(1600, 1));
    //STRESS(10000, test_tree_creation(16, 1));

    //test_subtree_mutation(16, 0);
    STRESS(10000, test_subtree_mutation(16, 1));

    printf("Exiting!\n");

    return 0;
}


void test_subtree_crossover(uint32_t table_size, int stress) {
    HashTable *t1 = create_hashtable(table_size);
    HashTable *t2 = create_hashtable(table_size);

    tree *parent1;
    tree *parent2;
    int d = 5;
    if (stress) {
        int m1 = rand() % 3;
        int m2 = rand() % 3;
        int md1 = rand() % (d + 1);
        int md2 = rand() % (d + 1);
        float tp1 = rand_float();
        float tp2 = rand_float();
        parent1 = generate_tree(t1, m1, 0, md1, tp1);
        parent2 = generate_tree(t1, m2, 0, md2, tp2);
        printf("(met, mdep, term_prob)- p1: %d %d %.6f | p2: %d %d %.6f \n", m1, md1, tp1, m2, md2, tp2);
    } else {
        int m1 = rand() % 3;
        int m2 = rand() % 3;
        int md1 = rand() % 6;
        int md2 = rand() % 6;
        float tp1 = rand_float();
        float tp2 = rand_float();
        printf("(met, mdep, term_prob)- p1: %d %d %.6f | p2: %d %d %.6f \n", m1, md1, tp1, m2, md2, tp2);
        m1 = 1;
        m2 = 2;
        md1 = 4;
        md2 = 0;
        tp1 = 0.453;
        tp2 = 0.322;
        parent1 = generate_tree(t1, m1, 0, md1, tp1);
        parent2 = generate_tree(t1, m2, 0, md2, tp2);

        //const char *str_tree1 = "div(div(add(add(x, y), add(x, x)), add(sub(x, y), div(scalar, scalar(-3.451, 0.021, 2.303)))), mul(sub(div(y, scalar(-1.390, -3.560, 1.223)), add(x, y)), add(div(scalar(2.418, -1.152, 4.430), y), div(scalar(1.098, -3.826, 4.727), y))))";
        //const char *str_tree1 = "div(scalar, scalar(-3.451, 0.021, 2.303))";
        //const char *str_tree2 = "x";
        //parent2 = str_to_tree(t2, str_tree2);
        //parent1 = str_to_tree(t1, str_tree1);
    }
    //tree *parent1 = generate_tree(t1, Full, mind, maxd, 0.33);
    //tree *parent2 = generate_tree(t1, Grow, mind, maxd, 0.33);

    if (!stress) printf("\nParent 1: \n");
    print_tree(parent1, 1 - stress, 1 - stress);
    if (!stress) printf("\nParent 2: \n");
    print_tree(parent2, 1 - stress, 1 - stress);

    tree *result;
    if (stress) {
        result = subtree_crossover(parent1, parent2, t2, d);
    } else {
        //tree *result = subtree_crossover_d(parent1, parent2, t2, 7, 3, All_Prims, All_Prims, 10);
        //tree *result = subtree_crossover_d(parent1, parent2, t2, 3, 1, All_Prims, All_Prims, 3);
        result = subtree_crossover_d(parent1, parent2, t2, 3, 1, 1, 1, d);
    }

    if (!stress) printf("\nResult: \n");
    print_tree(result, 1 - stress, 1 - stress);

    free_hashtable(t1);
    free_hashtable(t2);
    printf("Subtree crossover test ended!\n");
}


// tested with 0 depths
void test_tree_generation(uint32_t table_size, int method, int mind, int maxd) {
    HashTable *t1 = create_hashtable(table_size);
    print_dag_table(t1->table, table_size);

    tree *result = generate_tree(t1, method, mind, maxd, 0.33);
    printf("Expression:\n");
    printf("%s\n", get_dag_expr(result->dag));

    printf("\nCurrent Hashtable:\n");
    print_dag_table(t1->table, table_size);

    free_hashtable(t1);
    printf("Tree generation test ended!\n");
}


void test_tree_creation(uint32_t table_size, int stress) {
    HashTable *t1 = create_hashtable(table_size);
    //print_dag_table(t1->table, table_size);

    tree *result;
    if (!stress) {
        const char *genotype = "mul(add(scalar(1.0), scalar(1.0)), add(scalar(1.0), scalar(1.0)))";
        printf("Expression:\n%s\n", genotype);
        result = str_to_tree(t1, genotype);
    } else {
        int d = 8;
        int met = rand() % 3;
        int dep = rand() % (d + 1);
        int mdep = d;
        int term_prob = rand_float();
        result = generate_tree(t1, met, dep, mdep, term_prob);
    }
    
    if (!stress) printf("Confirm expression:\n");
    printf("%s\n", get_dag_expr(result->dag));
    if (!stress) printf("Now in fancy:\n");
    if (!stress) fancy_dag_print(result->dag);
    
    if (!stress) printf("Current Hashtable:\n");
    if (!stress) print_dag_table(t1->table, table_size);

    free_hashtable(t1);
    printf("Tree creation test ended!\n");
}


void test_subtree_mutation(uint32_t table_size, int stress) {
    HashTable *t1 = create_hashtable(table_size);
    HashTable *t2 = create_hashtable(table_size);

    tree *parent, *result;
    if (!stress) {
        const char *genotype = "mul(add(scalar(1.0), scalar(1.0)), add(scalar(1.0), scalar(1.0)))";
        parent = str_to_tree(t1, genotype);
        int c = 3;
        int m = All_Prims;
        result = subtree_mutation_d(parent, t2, c, m, 5);
    } else {
        int d = 6;
        int met = rand() % 3;
        int dep = rand() % (d + 1);
        int mdep = d;
        int term_prob = rand_float();
        parent = generate_tree(t1, met, dep, mdep, term_prob);
        result = subtree_mutation(parent, t2, d);
    }

    if (!stress) printf("\nParent: \n");
    print_tree(parent, 1 - stress, 1 - stress);
    if (!stress) printf("\nResult: \n");
    print_tree(result, 1 - stress, 1 - stress);

    if (!stress && 1) {
        printf("\nOffspring Hashtable:\n");
        print_dag_table(t2->table, table_size);
    }

    free_hashtable(t1);
    free_hashtable(t2);
    printf("Subtree mutation test ended!\n");
}

