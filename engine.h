#ifndef ENGINE_H
#define ENGINE_H

#define DIMS 2
#define INIT_MAX_GENS 1024
#define STAT_ARRS 3


// 32 bytes
typedef struct Stats {
    float *fit_hist[STAT_ARRS];
    float *dep_hist[STAT_ARRS];
    float *node_hist[STAT_ARRS];
    size_t hist_len;
} Stats;

//TODO: add vars like gen and timers to the engine struct
typedef struct Engine {
    uint32_t generations;
    uint32_t tournament_size;
    float term_prob;
    float mut_rate;
    float cross_rate;
    float elitism;
    uint32_t pop_size;
    int gen_method;
    struct allowed_depth {
        uint32_t max;
        uint32_t min; 
    } allowed_depth;
    struct init_depth {
        uint32_t max;
        uint32_t min;
    } init_depth;
    int debug;
    int max_retries;
    uint32_t initial_hashmap_size;
    int req_to_cache;
    uint32_t resolution[DIMS];
    float MIN_DOMAIN[DIMS];
    float MAX_DOMAIN[DIMS];
    float STEP[DIMS];
    float *target;
    uint32_t fitness_cases;
    Stats *stats;
    uint32_t best_ind_gen;
    uint32_t cur_gen;
} Engine;



float pagie_poly(float x, float y);
void set_default_params(Engine *params);
void init_stats(Stats *stats);
void free_stats(Stats *stats);
Engine *create_params(int set_default);
void setup(Engine *run, int cache_size);
void print_params(Engine *run);
void free_engine(Engine *run);
void cleanup(Engine *run);
uint64_t calculate_stats(Engine *run, tree **population);
void print_gen_statistics(Engine *run, tree **population, double duration);
int evolve(Engine *run);


#endif