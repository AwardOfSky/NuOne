#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include "gp.h"

// fitness func
domain_t pagie_poly(domain_t x, domain_t y) {
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
    run->max_retries = 10;
    run->initial_hashmap_size = INIT_HASHMAP_SIZE;
    run->elitism = 0.02;
    for (int i = 0; i < DIMS; ++i) {
        run->resolution[i] = 16;
        run->MIN_DOMAIN[i] = -1;
        run->MAX_DOMAIN[i] = 1;
    }
    run->cur_gen = 0;
    run->stats = (Stats *)malloc(sizeof(*(run->stats)));
    init_stats(run->stats);
}


void init_stats(Stats *stats) {
    stats->hist_len = INIT_MAX_GENS;
    for (int i = 0; i < STAT_ARRS; ++i) {
        stats->dep_hist[i] = (double *)malloc(INIT_MAX_GENS * sizeof(*(stats->dep_hist[i])));
        stats->fit_hist[i] = (double *)malloc(INIT_MAX_GENS * sizeof(*(stats->fit_hist[i])));
        stats->node_hist[i] = (double *)malloc(INIT_MAX_GENS * sizeof(*(stats->node_hist[i])));
    }
}


void free_stats(Stats *stats) {
    for (int i = 0; i < STAT_ARRS; ++i) {
        free(stats->dep_hist[i]);
        free(stats->fit_hist[i]);
        free(stats->node_hist[i]);
    }
}


Engine *create_params(int set_default) {
    Engine *run = (Engine *)malloc(sizeof(*run));
    if (set_default) {
        set_default_params(run);
    }
    return run;
}


void reset_cur_vars(Engine *run) {
    for(int i = 0; i < DIMS; ++i) {
        run->cur_vars[i] = run->MIN_DOMAIN[i];
    }
    run->index = 0;
    //printf("Reseting index, %d\n", run->index);
}


void setup(Engine *run, int cache_size) {

    // setup fitness cases, step
    run->fitness_cases = 1;
    for(uint32_t i = 0; i < DIMS; ++i) {
        run->STEP[i] = DOMAIN_DELTA(run, i) / (domain_t)(run->resolution[i] - 1);
        run->fitness_cases *= run->resolution[i];
    }
    run->target = (domain_t *)malloc(sizeof(*(run->target)) * run->fitness_cases);
    //run->domain = (domain_t *)malloc(sizeof(*(run->domain)) * run->fitness_cases);

    // init cache
    cache = create_cache(run->fitness_cases, cache_size);


    // init vars for target loop
    reset_cur_vars(run);
    
    // setup target
    for(uint32_t i = 0; i < run->fitness_cases; ++i) {
        run->target[i] = pagie_poly(run->cur_vars[0], run->cur_vars[1]);
        
        // setup vars in cache
        if (ENABLE_CACHE) {
            for (int j = 0; j < DIMS; ++j) {
                cache->vars[j * run->fitness_cases + i] = run->cur_vars[j];
            }
        }

        STEP_DOMAIN(run);
    }
}


void print_params(Engine *run) {
    printf("\n========== Engine parameters: ==========\n");
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
    free_cache(cache);
    free_stats(run->stats);
    free(run->stats);
    free(run);
}


void cleanup(Engine *run) {
    printf("Cleaning up...\n");
    free_engine(run);
    printf("Done!\n");
}


uint64_t calculate_stats(Engine *run, tree **population) {
    fit_t sum_fit = 0.0;
    Stats *s = run->stats;
    int sum_dep = 0;
    int sum_nodes = 0;
    int sum_prims = 0;
    int n = run->pop_size;
    int gen = run->cur_gen;

    for(int i = 0; i < n; ++i) {
        sum_fit += population[i]->fitness;
        sum_dep += population[i]->depth;
        sum_prims += population[i]->n_prims;
        sum_nodes += population[i]->n_prims + population[i]->n_terms;
    }
    double mean_fit = (double)sum_fit / (double)(n);
    double mean_dep = (double)sum_dep / (double)(n);
    double mean_nodes = (double)sum_nodes / (double)(n);
    double std_fit = 0.0;
    double std_dep = 0.0;
    double std_nodes = 0.0;

    for(int i = 0; i < n; ++i) {
        std_fit = pow(population[i]->fitness - mean_fit, 2);
        std_dep = pow(population[i]->depth - mean_dep, 2);
        std_nodes = pow((population[i]->n_terms + population[i]->n_prims) - mean_nodes, 2);
    }


    s->fit_hist[Average][gen] = mean_fit;
    s->fit_hist[StandardDev][gen] = (double)sqrt(std_fit / (double)(n));
    s->fit_hist[Best][gen] = population[run->best_ind_gen]->fitness;

    s->dep_hist[Average][gen] = mean_dep;
    s->dep_hist[StandardDev][gen] = (double)sqrt(std_dep / (double)(n));
    s->dep_hist[Best][gen] = population[run->best_ind_gen]->depth;

    s->node_hist[Average][gen] = mean_nodes;
    s->node_hist[StandardDev][gen] = (double)sqrt(std_nodes / (double)(n));
    s->node_hist[Best][gen] = population[run->best_ind_gen]->n_terms + population[run->best_ind_gen]->n_prims;
    return run->fitness_cases * sum_prims;
}


void print_gen_statistics(Engine *run, tree **population, double duration) {
    duration = max(FLT_MIN, duration);
    uint64_t gpops = (uint64_t)(calculate_stats(run, population) / duration);
    Stats *s = run->stats;
    uint32_t g = run->cur_gen;
    if (!g) {
        printf("\n\n");
        printf("[       |              FITNESS             |               DEPTH              |                     NODES                        | TIMINGS(s) ]\n");
        printf("[  gen  |    avg    ,    std    ,   best   |    avg    ,    std    ,   best   |    avg    ,    std    ,   best   ,   GPops(/s)   |    eval    ]\n");
    }
    printf("[ %-6d| %-10.6f, %-10.6f, %-9.6f| %-10.6f, %-10.6f, %-9d| %-10.6f, %-10.6f, %-9d, %-14I64d| %-10.3f ]\n",
        run->cur_gen, s->fit_hist[Average][g], s->fit_hist[StandardDev][g], s->fit_hist[Best][g],
        s->dep_hist[Average][g], s->dep_hist[StandardDev][g], (int)(s->dep_hist[Best][g]),
        s->node_hist[Average][g], s->node_hist[StandardDev][g], (int)(s->node_hist[Best][g]), gpops, duration);
}


int evolve(Engine *run) {

    // general run vars
    clock_t engine_start = clock();
    int debug = run->debug;
    int elitism_n = run->elitism * run->pop_size;
    //if(debug) printf("Elitism n: %d\n", elitism_n);

    //create and evaluate initial pop
    HashTable *curt = create_hashtable(run->initial_hashmap_size);
    HashTable *newt;

    double gen_duration;
    tree **population = generate_population(curt, run->gen_method,
                                            run->init_depth.min, run->init_depth.max,
                                            run->pop_size, run->term_prob);

    if (debug) print_population(population, run->pop_size, 0, 0);
    
    tree **new_pop = (tree **)malloc(sizeof(tree *) * run->pop_size);
    tree **best_trees = (tree **)malloc(sizeof(tree *) * elitism_n);

    if (ENABLE_CACHE) {
        if (debug) {
            printf("Building cache...");
            //print_candidate_list(curt, 1);
        }
        //build_cache(curt, 0);
        if (debug) {
            printf("Done!\n");
            //print_candidate_list(curt, 1);
            //print_cache(0);
        }
    }

    clock_t gen_start = clock();
    run->best_ind_gen = calc_pop_fit(run, population);
    gen_duration = ((double) (clock() - gen_start)) / CLOCKS_PER_SEC;
    print_gen_statistics(run, population, gen_duration);

    //run for generations (only stop criteria for now)
    for(run->cur_gen = 1; run->cur_gen < run->generations; run->cur_gen++) {

        // generate new hashtable
        uint32_t new_size = next_power_2(curt->n_nodes) - 1;
        newt = create_hashtable(new_size);
        if (PDEBUG) printf("Generated new hashtable, size %d.\n", new_size);

        // and insert best in next pop
        best_trees = get_k_min_trees(population, run->pop_size, elitism_n);
        if (PDEBUG) {
            printf("Best trees: ");
            print_population(best_trees, elitism_n, run->cur_gen, 0);
        }

        for(int i = 0; i < elitism_n; ++i) {
            if (PDEBUG) printf("Inserting elitist individuals in new generation\n");
            new_pop[i] = copy_tree(best_trees[i], newt);
        }

        // recombine remaining individuals
        for(int ind = elitism_n; ind < run->pop_size; ++ind) {
            if (PDEBUG) printf("\nComputing individual %d\n", ind);
            int depth = -1;
            int retries = 0;

            tree *parent1;
            tree *indiv = NULL;
            while ((depth > run->allowed_depth.max || depth < run->allowed_depth.min) && retries < run->max_retries) {
                
                parent1 = tournament(run, population);

                printfd("\nParent 1: ");
                if (PDEBUG) print_tree(parent1, 0, 0);

                int rng = rand_float();
                if(rng < run->cross_rate) {
                    tree *parent2 = tournament(run, population);
                    printfd("Parent 2: ");
                    if (PDEBUG) print_tree(parent2, 0, 0);
                    indiv = subtree_crossover(parent1, parent2, newt, run->allowed_depth.max);
                    printfd("Crossover\n"); 
                } else if (rng < run->mut_rate + run->cross_rate) {
                    indiv = mutation(parent1, newt, run->allowed_depth.max);
                    printfd("Mutation\n"); 
                } else {
                    indiv = copy_tree(parent1, newt);
                    printfd("Copy tree\n"); 
                }


                printfd("Resulting indiv: ");
                if (PDEBUG) print_tree(indiv, 1, 0);

                depth = indiv->depth;
                ++retries;
            }

            if (indiv != NULL) {
                new_pop[ind] = indiv;
            } else {
                printf(PREERR"Evolve: Failed to generate individual index %d of new pop.\n", ind);
            }

            if (PDEBUG) {
                printf("Retries for ind %d: %d\n", ind, retries);
                //print_candidate_list(newt, 1);
                //print_cache(0);
            }
        }


        // calculate fitness of new population
        if (ENABLE_CACHE) {
            //build_cache(newt, gen);
            if (PDEBUG) printf("Building cache!\n");
        }
        gen_start = clock();
        run->best_ind_gen = calc_pop_fit(run, new_pop);
        gen_duration = ((double) (clock() - gen_start)) / CLOCKS_PER_SEC;
        print_gen_statistics(run, new_pop, gen_duration);

        // advance hastable and delete current table
        HashTable *tempt = curt;
        curt = newt;
        free_hashtable(tempt);

        // switch populations
        tree **tempp = population;
        population = new_pop;
        new_pop = tempp;
    }

    double engine_dur = ((double) (clock() - engine_start)) / CLOCKS_PER_SEC;
    printf("\nBest individual (fitness %.3f): %s\n", population[run->best_ind_gen]->fitness, get_dag_expr(population[run->best_ind_gen]->dag));
    printf("\nTotal engine time: %.3fs\n\n", engine_dur);


    return 0;
}