#include <stdio.h>
#include <stdlib.h>
#include "gp.h"


// fitness func
float pagie_poly(float x, float y) {
    return ((1 / (1 + (1 / (x * x * x * x)))) + (1 / (1 + (1 / (y * y * y * y)))));
}


void set_default_params(Engine *run) {
    run->generations = 50;
    run->tournament_size = 3;
    run->term_prob = 0.33;
    run->mut_rate = 0.1;
    run->cross_rate = 0.9;
    run->pop_size = 50;
    run->gen_method = Ramped;
    run->allowed_depth.min = 1;
    run->allowed_depth.max = 15;
    run->init_depth.min = 1;
    run->init_depth.max = 15;
    run->debug = 0;
    run->max_retries = 3;
    run->initial_hashmap_size = 131071;
    run->elitism = 0.02;
    for (int i = 0; i < DIMS; ++i) {
        run->resolution[i] = 16;
        run->MIN_DOMAIN[i] = -1;
        run->MAX_DOMAIN[i] = 1;
    }
}


Engine *create_params(int set_default) {
    Engine *run = (Engine *)malloc(sizeof(*run));
    if (set_default) {
        set_default_params(run);
    }
    return run;
}


void setup(Engine *run, int cache_size) {

    // setup fitness cases, step
    run->fitness_cases = 1;
    for(uint32_t i = 0; i < DIMS; ++i) {
        run->STEP[i] = DOMAIN_DELTA(run, i) / (float)(run->resolution[i] - 1);
        run->fitness_cases *= run->resolution[i];
    }
    run->target = (float *)malloc(sizeof(*(run->target)) * run->fitness_cases);
    //run->domain = (float *)malloc(sizeof(*(run->domain)) * run->fitness_cases);

    // init cache
    // run->cache = create_cache(cache_size);


    // init vars for target loop
    float cur_vars[DIMS];
    for (uint32_t i = 0; i < DIMS; ++i) {
        cur_vars[i] = run->MIN_DOMAIN[i];
    }
    
    // setup target
    for(uint32_t i = 0; i < run->fitness_cases; ++i) {
        run->target[i] = pagie_poly(cur_vars[0], cur_vars[1]);
        
        // setup vars in cache
        //for (int j = 0; j < DIMS; ++j) {
        //   run->cache->vars[j * run->fitness_cases + i] = run->cur_vars_g[j];
        //}
        //cache->vars[0] = cur_vars_g[0];

        STEP_DOMAIN(run, cur_vars);
    }
}


void print_params(Engine *run) {
    printf("\n========== Evolution parameters: ==========\n");
    printf("Generations:\t%u\n", run->generations);
    printf("Tour size:\t%u\n", run->tournament_size);
    printf("Terminal prob:\t%.3f\n", run->term_prob);
    printf("Mutation prob:\t%.3f\n", run->mut_rate);
    printf("Crossover prob:\t%.3f\n", run->cross_rate);
    printf("Pop size:\t%u\n", run->pop_size);
    printf("Pop gen method:\t");
    int gen_method = run->gen_method;
    if (gen_method == Ramped) {
        printf("Ramped\n");
    } else if (gen_method == Grow) {
        printf("Grow\n");
    } else {
        printf("Full\n");
    }
    printf("Min allow dep:\t%u\n", run->allowed_depth.min);
    printf("Max allow dep:\t%u\n",run->allowed_depth.max);
    printf("Min init dep:\t%u\n", run->init_depth.min);
    printf("Max init dep:\t%u\n", run->init_depth.max);
    printf("Debug level:\t%d\n", run->debug);
    printf("Max retries:\t%d\n", run->max_retries);
    printf("Init map size:\t%u\n", run->initial_hashmap_size);
    printf("Elitism %%:\t%.3f\n", run->elitism);
    printf("Problem dimensionality: %d\n", DIMS);
    printf("Domain Range: from (");
    for (uint32_t i = 0; i < DIMS; ++i) {
        printf("%.3f", run->MIN_DOMAIN[i]);
        if (i != DIMS - 1) printf(", ");
    }
    printf(") to (");
    for (uint32_t i = 0; i < DIMS; ++i) {
        printf("%.3f", run->MAX_DOMAIN[i]);
        if (i != DIMS - 1) printf(", ");
    }
    printf(")\n");
    printf("Domain resolution: [");
    for (uint32_t i = 0; i < DIMS; ++i) {
        printf("%d", run->resolution[i]);
        if (i != DIMS - 1) printf(", ");
    }
    printf("]\n");
    printf("Domain step: [");
    for (uint32_t i = 0; i < DIMS; ++i) {
        printf("%.3f", run->STEP[i]);
        if (i != DIMS - 1) printf(", ");
    }
    printf("]\n");
    printf("Total number of fitness cases: %d\n", run->fitness_cases);
    printf("===========================================\n");
}


void free_engine(Engine *run) {
    // free(run->cache);
    free(run);
}


void cleanup(Engine *run) {
    printf("Cleaning up...\n");
    free_engine(run);
    printf("Done!\n");
}


void print_gen_statistics(int gen, tree **population, int best_index) {
    printf("[Gen %d]\tBest fitness: %f on index %d\n", gen, population[best_index]->fitness, best_index);
}



int evolve(Engine *run) {

    // general run vars
    int best_index;
    int debug = run->debug;
    int elitism_n = run->elitism * run->pop_size;
    //if(debug) printf("Elitism n: %d\n", elitism_n);

    //create and evaluate initial pop
    HashTable *curt = create_hashtable(run->initial_hashmap_size);
    HashTable *newt;

    tree **population = generate_population(curt, run->gen_method,
                                            run->init_depth.min, run->init_depth.min,
                                            run->pop_size, run->term_prob);
    if (debug) print_population(population, run->pop_size, 0, 0);
    tree **new_pop = (tree **)malloc(sizeof(tree *) * run->pop_size);
    tree **best_trees = (tree **)malloc(sizeof(tree *) * elitism_n);

    if (ENABLE_CACHE) {
        if (debug) {
            printf("Building cache!\n");
            //print_candidate_list(curt, 1);
        }
        //build_cache(curt, 0);
        if (debug) {
            printf("Built cache!\n");
            //print_candidate_list(curt, 1);
            //print_cache(0);
        }
    }


    best_index = calc_pop_fit(run, population);
    print_gen_statistics(0, population, best_index);

    //run for generations (only stop criteria for now)
    for(int gen = 1; gen < run->generations; ++gen) {

        // generate new hashtable
        uint32_t new_size = next_power_2(curt->n_nodes) - 1;
        newt = create_hashtable(new_size);
        if (debug) printf("Generated new hashtable, size %d.\n", new_size);

        // and insert best in next pop
        best_trees = get_k_min_trees(population, run->pop_size, elitism_n);
        if (debug) {
            printf("Best trees: ");
            print_population(best_trees, elitism_n, gen, 0);
        }

        for(int i = 0; i < elitism_n; ++i) {
            new_pop[i] = copy_tree(best_trees[i], newt);
        }
        if (debug) printf("Inserting elitist individuals in new generation\n");

        // recombine remaining individuals
        for(int ind = elitism_n; ind < run->pop_size; ++ind) {
            if (debug) printf("\nComputing individual %d\n", ind);
            int depth = -1;
            int retries = 0;

            tree *parent1;
            tree *indiv = NULL;
            while ((depth > run->allowed_depth.max || depth < run->allowed_depth.min) && retries < 3) {
                
                parent1 = tournament(run, population);
                if (debug) printf("Tournament sel for parent 1\n"); 

                int rng = rand_float();
                if(rng < run->cross_rate) {
                    tree *parent2 = tournament(run, population);
                    if (debug) printf("Tournament sel for parent 2\n"); 
                    indiv = subtree_crossover(parent1, parent2, newt, run->allowed_depth.max);
                    if (debug) printf("Crossover\n"); 
                } else if (rng < run->mut_rate + run->cross_rate) {
                    indiv = mutation(parent1, newt, run->allowed_depth.max);
                    if (debug) printf("Mutation\n"); 
                } else {
                    indiv = copy_tree(parent1, newt);
                    if (debug) printf("Copy tree\n"); 
                }

                depth = indiv->depth;
                ++retries;
            }

            if (indiv != NULL) {
                new_pop[ind] = indiv;
            } else {
                printf(PREERR"Evolve: Failed to generate individual index %d of new pop.\n", ind);
            }

            if(debug) {
                printf("Retries for ind %d: %d\n", ind, retries);
                //print_candidate_list(newt, 1);
                //print_cache(0);
            }
        }


        // calculate fitness of new population
        if (ENABLE_CACHE) {
            //build_cache(newt, gen);
            if (debug) printf("Building cache!\n");
        }
        best_index = calc_pop_fit(run, new_pop);
        print_gen_statistics(gen, new_pop, best_index);

        // advance hastable and delete current table
        HashTable *tempt = curt;
        curt = newt;
        free_hashtable(tempt);

        // switch populations
        tree **tempp = population;
        population = new_pop;
        new_pop = tempp;
    }

    return 0;
}