#ifndef UTILS_H
#define UTILS_H

#include <stdarg.h>
#include <stdint.h>
#include <math.h>
#include "rsort.h"

#define PREERR "[Error]:\t"
#define PREWARN "[Warning]:\t"

#define STACK_ALLOC 1024
#define BINSEARCH_TRESHOLD 32
#define NUM_THREADS 8
// bits in a byte
#define CHAR_BITS 8

#define max(a,b)                \
({ __typeof__ (a) _a = (a);     \
    __typeof__ (b) _b = (b);    \
    _a > _b ? _a : _b; })


#define min(a,b)                \
({ __typeof__ (a) _a = (a);     \
    __typeof__ (b) _b = (b);    \
    _a < _b ? _a : _b; })


#define abs(a)                  \
({ __typeof__ (a) _a = (a);     \
    _a > 0 ? _a : -_a; })


#define clip(a, n, b)           \
({ __typeof__ (a) _a = (a);     \
    __typeof__ (b) _b = (b);    \
    __typeof__ (n) _n = (n);    \
    _b > _n ? ((_a < _n) ? _n : _a) : _b; })
    

#define PRINT_ARR(ARR, N, FORMAT) do {  \
    for(size_t _i = 0; _i < N; ++_i) {  \
        printf(FORMAT" ", ARR[_i]);     \
    }                                   \
    printf("\n");                       \
} while(0)


// Computes the standard deviation of a list of type T
#define DECLARE_COMPUTE_STD(T)                          \
    double compute_std_##T(T *data, const size_t n, T sum)
#define MAKE_COMPUTE_STD(T)                             \
    DECLARE_COMPUTE_STD(T) {                            \
        double mean = (double)sum / (double)n;          \
        double SD = 0.0;                                \
        for (size_t i = 0; i < n; ++i) {                \
            SD += pow(data[i] - mean, 2);               \
        }                                               \
        return sqrt(SD / (double)n);                    \
    }


// Computes the mean of a list of type T, is last argument is not NULL, store sum in "sum"
#define DECLARE_COMPUTE_MEAN(T)                         \
    double compute_mean_##T(T *data, const size_t n, T *sum)
#define MAKE_COMPUTE_MEAN(T)                            \
    DECLARE_COMPUTE_MEAN(T)  {                          \
        T s = (T)0;                                     \
        for(size_t i = 0; i < n; ++i) {                 \
            s += data[i];                               \
        }                                               \
        if (sum != NULL) {                              \
            *sum = s;                                   \
        }                                               \
        return (double)(s) / (double)(n);               \
    }


// if (!strcmp("fit_t", #T)) {
//     printf("We are on the float thing\n");
// }

// reallocing arrays of arbitrary type. Note: realloc returns NULL if size is 0 but no error, hence the check

// #define DECLARE_REALLOC(T)                                                                              
//     T *realloc_##T(T *data, const uint32_t new_size, uint32_t *sptr)
// #define MAKE_REALLOC(T)                                                                                 
//     DECLARE_REALLOC(T) {                                                                                
//         T *tp = (T*)realloc(data, sizeof(*tp) * new_size);                                              
//         if (tp != NULL || new_size == 0) {                                                              
//             data = tp;                                                                                  
//             if (sptr != NULL) {                                                                         
//                 *sptr = new_size;                                                                       
//             }                                                                                           
//         } else {                                                                                        
//             printf(PREERR"OOM: Failed to realloc array of type "#T", on memory address %p. ", data);    
//             printf("Requested %lu bytes.\n", sizeof(*data) * new_size);                               
//         }                                                                                               
//         return data;                                                                                    
//     }                            


#define DECLARE_REALLOC(T)                                                                              \
    void realloc_##T(T **data, const uint32_t new_size, uint32_t *sptr)
#define MAKE_REALLOC(T)                                                                                 \
    DECLARE_REALLOC(T) {                                                                                \
        T *tp = (T*)realloc(*data, sizeof(*tp) * new_size);                                              \
        if (tp != NULL || new_size == 0) {                                                              \
            *data = tp;                                                                                  \
            if (sptr != NULL) {                                                                         \
                *sptr = new_size;                                                                       \
            }                                                                                           \
        } else {                                                                                        \
            printf(PREERR"OOM: Failed to realloc array of type ""char"", on memory address %p. ", *data);\
            printf("Requested %lu bytes.\n", sizeof(**data) * new_size);                            \
        }                                                                                               \
    }      \



//DECLARE_COMPUTE_STD(int);
DECLARE_COMPUTE_STD(float);
//DECLARE_COMPUTE_MEAN(int);
//DECLARE_COMPUTE_MEAN(float);
DECLARE_REALLOC(float);
DECLARE_REALLOC(char);

float rand_float();
// int vasprintf(char **strp, const char *format, va_list ap);
// int asprintf(char **strp, const char *format, ...);
uint32_t next_power_2(uint32_t n);
int *range_random_sample(int n, int size);
int *sel_k_min(int *arr, int n, int size);


#endif
