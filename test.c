//gcc -o test.exe -Wall -Ofast test.c gp.c hashmap.c utils.c genetics.c engine.c cache.c
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


#define CLOSE_(PRT_TAB) do {                    \
    if (!stress) printf("\nParent: \n");        \
    print_tree(parent, 1, 1 - stress);          \
    if (!stress) printf("\nResult: \n");        \
    print_tree(result, 1, 1 - stress);          \
    if (!stress && PRT_TAB) {                   \
        printf("\nOffspring Hashtable:\n");     \
        print_dag_table(t2->table, table_size); \
    }                                           \
    free_hashtable(t1);                         \
    free_hashtable(t2);                         \
} while(0)


void test_tree_creation(uint32_t table_size, int stress);
void test_tree_generation(uint32_t table_size, int method, int mind, int maxd);
void test_subtree_crossover(uint32_t table_size, int stress);
void test_subtree_mutation(uint32_t table_size, int stress);
void test_delete_mutation(uint32_t table_size, int stress);
void test_insert_mutation(uint32_t table_size, int stress);
void test_list_same_arity();
void test_sel_min(int n, int size);
void test_point_mutation(uint32_t table_size, int stress);
void test_evolution_v1();
void test_evolution_v2();
void test_realloc();
void test_cache(size_t table_size);


int main() {

    //unsigned int seed = 1644863329; // reproducibility
    unsigned int seed = 1645327584; // reproducibility
    //unsigned int seed = time(NULL); // reproducibility

    srand(seed);
    printf("Random seed: %u\n", seed);

    
    //test_tree_generation(16, Full, 0, 0);
    
    //test_subtree_crossover(16, 0);
    //STRESS(10000, test_subtree_crossover(1600, 1));
    //STRESS(10, test_tree_creation(16, 1));

    //test_subtree_mutation(16, 0);
    //STRESS(10000, test_subtree_mutation(16, 1));
    
    //test_delete_mutation(16, 0);
    //test_insert_mutation(16, 0);
    //test_list_same_arity();

    //test_point_mutation(16, 0);
    //STRESS(10, test_point_mutation(16, 1));
    

    //STRESS(10000, test_tree_creation(16, 1)); // passed
    //STRESS(10000, test_subtree_crossover(16, 1)); //passed
    //STRESS(10000, test_subtree_mutation(16, 1)); //passed
    //STRESS(10000, test_delete_mutation(16, 1)); //passed
    //STRESS(10000, test_insert_mutation(16, 1)); //passed
    //STRESS(10000, test_point_mutation(16, 1)); //passed

    //test_point_mutation(16, 0);

    // test_sel_min(10000, 10);
    
    //test_evolution_v1();
    test_evolution_v2();
    //test_realloc();
    
    //test_cache(16);

    printf("Exiting!\n");

    return 0;
}


void test_cache(size_t table_size) {
    uint32_t cache_s = 16;
    // tested shrink with cache size 16 (factor grow, shrink 2, 4) and grow with 3 (factor grow, shrink 2, 4)

    init_engine(1);
    for (int i = 0; i < DIMS; ++i) {
        run.resolution[i] = 4;
        run.MIN_DOMAIN[i] = -5;
        run.MAX_DOMAIN[i] = 5;
    }
    setup(cache_s);
    //print_domain(run.target, run.fitness_cases, run.resolution);
    //print_params(run);

    HashTable *t1 = create_hashtable(table_size);
    HashTable *t2 = create_hashtable(table_size);
    //const char *genotype = "mul(sin(x), sin(x))";
    //const char *genotype = "add(if(sin(x), sin(x), tan(y)), mul(add(y, y), add(y, y)))";
    //const char *genotype1 = "sub(mul(sin(x), add(y, y)), tan(y))";

    const char *genotype = "add(if(sin(x), sin(x), tan(y)), if(sin(x), tan(y), add(y, y)))";
    const char *genotype1 = "sub(mul(sin(x), add(y, y)), mul(add(y, y), add(y, y)))";

    tree *parent = str_to_tree(t1, genotype);
    tree *parent1 = str_to_tree(t1, genotype1);

    print_dag_table(t1->table, t1->size);
    print_candidate_list(t1, 1);


    printf("Starting gen 1!\n");
    printf("Building cache...\n");
    build_cache(t1, 0);
    printf("Done!\n\n");
    printf("Updated node list: ");
    print_candidate_list(t1, 1);
    fit_t fitness;
    printf("First tree\n");
    print_tree(parent, 0, 0);
    fitness = calc_tree_fit(parent);
    
    printf("Second tree\n");
    print_tree(parent1, 0, 0);
    fitness = calc_tree_fit(parent1);
    
    print_cache(run.resolution, 1);

    printf("Starting gen 2!\n");
    //const char *gen2genotype = "add(if(sin(x), tan(y), add(y, y)), if(sin(x), tan(y), add(y, y)))";
    //tree *parent2 = str_to_tree(t2, gen2genotype);

    tree *parent2 = subtree_crossover_d(parent, parent, t2, 2, 9, All_Prims, All_Prims, 10);
    printf("\n");
    print_tree(parent2, 1, 0);

    printf("Building cache gen 2...\n");
    build_cache(t2, 1);
    printf("Done!\n\n");
    printf("Updated node list: ");
    print_candidate_list(t2, 1);
    fitness = calc_tree_fit(parent2);
    printf("Done!\n");
    print_cache(run.resolution, 1);


    printf("The fitness was (do not optimize): %.3f\n", fitness);
    printf("\nFreeing cache:...");
    printf("Done!\nFreeing other stuff...");
    free(t1);
    free(t2);
    free(parent);
    cleanup();
    printf("Done!\nCleanup done, let's exit!\n");
}


void test_realloc() {
    int size = 5; 
    float *arr = malloc(size * sizeof(float));
    printf("Size of array: %d\n", size);
    arr = realloc_float(arr, 100000000000, (size_t *)(&size));
    //arr = realloc(arr, 100000000000)
    //arr = realloc_float(arr, 0, (size_t *)(&size));
    //arr = realloc(arr, 0)
    printf("Size of array: %d\n", size);
    free(arr);
}


void test_evolution_v1() {
    init_engine(1);

    // specific vars
    run.generations = 5;
    run.pop_size = 20;
    run.tournament_size = 3;
    run.gen_method = Ramped;
    run.allowed_depth.min = 2;
    run.allowed_depth.max = 10;
    run.init_depth.min = 2;
    run.init_depth.max = 6;
    run.debug = 1;
    run.elitism = 0.2;
    for (int i = 0; i < DIMS; ++i) {
        run.resolution[i] = 256;
        run.MIN_DOMAIN[i] = -5;
        run.MAX_DOMAIN[i] = 5;
    }

    setup(10);

    //print_domain(run.target, run.fitness_cases, run.resolution);
    print_params();
    
    evolve();
    cleanup();
}


void test_evolution_v2() {
    init_engine(1);

    // specific vars
    run.generations = 1000;
    run.pop_size = 100;
    run.tournament_size = 3;
    run.gen_method = Ramped;
    run.allowed_depth.min = 2;
    run.allowed_depth.max = 25;
    run.init_depth.min = 2;
    run.init_depth.max = 6;
    run.debug = 1;
    run.elitism = 0.2;
    for (int i = 0; i < DIMS; ++i) {
        run.resolution[i] = 256;
        run.MIN_DOMAIN[i] = -5;
        run.MAX_DOMAIN[i] = 5;
    }

    setup(10);

    //print_domain(run.target, run.fitness_cases, run.resolution);
    print_params();
    
    evolve();
    cleanup();
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
        parent = generate_tree(t1, met, 0, dep, term_prob);
        result = subtree_mutation(parent, t2, mdep);
    }

    CLOSE_(1);
    printf("Subtree mutation test ended!\n");
}


void test_delete_mutation(uint32_t table_size, int stress) {
    HashTable *t1 = create_hashtable(table_size);
    HashTable *t2 = create_hashtable(table_size);

    tree *parent, *result;
    if(!stress) {
        //const char *genotype = "mul(add(scalar(1.0), scalar(1.0)), add(scalar(1.0), scalar(1.0)))";
        //const char *genotype = "add(mul(x, y), y)";
        //const char *genotype = "mul(x, y)";
        const char *genotype = "y";

        parent = str_to_tree(t1, genotype);
        int c = 2;
        result = delete_mutation_d(parent, t2, c);

    } else {
        int d = 6;
        int met = rand() % 3;
        int dep = rand() % (d + 1);
        int term_prob = rand_float();
        parent = generate_tree(t1, met, 0, dep, term_prob);
        result = delete_mutation(parent, t2);
    }

    CLOSE_(0);
    printf("Delete mutation test ended!\n");
}


void test_insert_mutation(uint32_t table_size, int stress) {
    HashTable *t1 = create_hashtable(table_size);
    HashTable *t2 = create_hashtable(table_size);

    tree *parent, *result;
    if(!stress) {
        //const char *genotype = "mul(add(scalar(1.0), scalar(1.0)), add(scalar(1.0), scalar(1.0)))";
        const char *genotype = "add(mul(y, y), y)";
        //const char *genotype = "mul(x, y)";
        //const char *genotype = "y";
        parent = str_to_tree(t1, genotype);
        int c = 2;
        result = insert_mutation_d(parent, t2, c, All_Prims, 2);
    } else {
        int d = 6;
        int met = rand() % 3;
        int dep = rand() % 6;
        int mdep = d;
        int term_prob = rand_float();
        parent = generate_tree(t1, met, 0, dep, term_prob);
        result = insert_mutation(parent, t2, mdep);
    }

    CLOSE_(0);
    printf("Insert mutation test ended!\n");
}


// PRINT_ARR(ARR, N, FORMAT)
// const Prim primitive_set[] = {
//     {"placeholder", 0, 0}, // placeholder for 0 ID
//     {"scalar", 3, Scalar},
//     {"x", 0, X},
//     {"y", 0, Y},
//     {"add", 2, Add},
//     {"sub", 2, Sub},
//     {"div", 2, Div},
//     {"mul", 2, Mul}
// };
// "exclude" - exclude mode
// (0 -> include all primitives with that arity)
// (1 -> exclude all primitives with that arity)
// ignore_specific to ignore a specific primitive, ignore_specific = -1 if you do not wish to ignore anything
// 
// search_set - which sets to search on
// (0 -> terminal set, only search for terminals)
// (1 -> function set, only search for functions)
// (>=2 -> search for both functions and terminals)
//enum PSet_IDs {Scalar = 1, X = 2, Y = 3, Add = 4, Sub = 5, Div = 6, Mul = 7, PSET_START = 1, TSET_END = 3, PSET_END = 7};
void test_list_same_arity() {
    int n;
    prim_index_type *cands = get_list_same_arity(0, &n, 0, 10, 0);
    printf("N is : %d\n", n);
    PRINT_ARR(cands, n, "%d");
    printf("\n");
    free(cands);

    cands = get_list_same_arity(0, &n, 0, 2, 0);
    printf("N is : %d\n", n);
    PRINT_ARR(cands, n, "%d");
    printf("\n");
    free(cands);

    cands = get_list_same_arity(2, &n, 0, 2, 0);
    printf("N is : %d\n", n);
    PRINT_ARR(cands, n, "%d");
    printf("\n");
    free(cands);

    cands = get_list_same_arity(2, &n, 0, 2, 1); // also for 2
    printf("N is : %d\n", n);
    PRINT_ARR(cands, n, "%d");
    printf("\n");
    free(cands);

    cands = get_list_same_arity(0, &n, 1, 10, 0);
    printf("N is : %d\n", n);
    PRINT_ARR(cands, n, "%d");
    printf("\n");
    free(cands);

    cands = get_list_same_arity(0, &n, 1, 1, 0);
    printf("N is : %d\n", n);
    PRINT_ARR(cands, n, "%d");
    printf("\n");
    free(cands);

    cands = get_list_same_arity(2, &n, 1, 2, 0);
    printf("N is : %d\n", n);
    PRINT_ARR(cands, n, "%d");
    printf("\n");
    free(cands);

    cands = get_list_same_arity(2, &n, 1, 2, 1); // also for 2
    printf("N is : %d\n", n);
    PRINT_ARR(cands, n, "%d");
    printf("\n");
    free(cands);
    
    printf("End of prim test!\n");
}


void test_sel_min(int n, int size) {
    int *array = (int *)malloc(sizeof(*array) * n);
    for(int i = 0; i < n; ++i) {
        array[i] = rand() % RAND_MAX;
    }
    int *result = sel_k_min(array, n, size);
    PRINT_ARR(array, min(n, 10), "%d");
    PRINT_ARR(result, min(size, 10), "%d");
    free(array);
    free(result);
}


void test_point_mutation(uint32_t table_size, int stress) {
    HashTable *t1 = create_hashtable(table_size);
    HashTable *t2 = create_hashtable(table_size);

    tree *parent, *result;
    if(!stress) {
        //const char *genotype = "mul(add(scalar(1.0), scalar(1.0)), add(scalar(1.0), scalar(1.0)))";
        const char *genotype = "add(mul(y, y), y)";
        //const char *genotype = "mul(x, y)";
        //const char *genotype = "y";
        parent = str_to_tree(t1, genotype);
        int c = 3;
        result = point_mutation_d(parent, t2, c, All_Prims);
    } else {
        int d = 6;
        int met = rand() % 3;
        int dep = rand() % (d + 1);
        int term_prob = rand_float();
        parent = generate_tree(t1, met, 0, dep, term_prob);
        //printf("Here we are!\n");
        result = point_mutation(parent, t2);
        //printf("1Here we are!\n");
    }

    CLOSE_(0);
    printf("Point mutation test ended!\n");
}