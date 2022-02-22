#ifndef UTILS_H
#define UTILS_H

#include <stdarg.h>
#include <stdint.h>

#define PREERR "[Error]:\t"
#define PREWARN "[Warning]:\t"

#define STACK_ALLOC 1024
#define BINSEARCH_TRESHOLD 32

#define max(a,b)                \
({ __typeof__ (a) _a = (a);     \
    __typeof__ (b) _b = (b);    \
    _a > _b ? _a : _b; })


#define min(a,b)                \
({ __typeof__ (a) _a = (a);     \
    __typeof__ (b) _b = (b);    \
    _a < _b ? _a : _b; })


#define clip(a, n, b)           \
({ __typeof__ (a) _a = (a);     \
    __typeof__ (b) _b = (b);    \
    __typeof__ (n) _n = (n);    \
    _b > _n ? ((_a < _n) ? _n : _a) : _b; })
    

#define PRINT_ARR(ARR, N, FORMAT) do {  \
    for(size_t i = 0; i < N; ++i) {     \
        printf(FORMAT" ", ARR[i]);      \
    }                                   \
    printf("\n");                       \
} while(0)


#define RADIX_STRUCT(STRUCT_TYPE, ARR, MEMBER) do {                 \
    STRUCT_TYPE *h = (STRUCT_TYPE *)malloc(n * sizeof(*h));         \
    STRUCT_TYPE *pointers[0x100];                                   \
    STRUCT_TYPE *s;                                                 \
    int swap = 0;                                                   \
    for(int mbit = 0; mbit < sizeof(int) << 3; mbit += 8) {         \
        int bucket[0x100] = {0};                                    \
        int i;                                                      \
        for (i = 0; i < n; ++i) {                                   \
            ++bucket[(ARR[i].MEMBER >> mbit) & 0xFF];               \
        }                                                           \
        s = h;                                                      \
        int next = 0;                                               \
        for (i = 0; i < 0x100; ++i) {                               \
            if(bucket[i] == n) {                                    \
                next = 1;                                           \
                break;                                              \
            }                                                       \
        }                                                           \
        if (next) {                                                 \
            continue;                                               \
        }                                                           \
        for (i = 0xFF; i >= 0; s += bucket[i--]) {                  \
            pointers[i] = s;                                        \
        }                                                           \
        for (i = 0; i < n; ++i) {                                   \
            *pointers[(ARR[i].MEMBER >> mbit) & 0xFF]++ = ARR[i];   \
        }                                                           \
        STRUCT_TYPE *temp = h;                                      \
        h = ARR;                                                    \
        ARR = temp;                                                 \
        swap = 1 - swap;                                            \
    }                                                               \
                                                                    \
    if(swap) {                                                      \
        memcpy(h, ARR, sizeof(*h) * n);                             \
        free(ARR);                                                  \
    } else {                                                        \
        free(h);                                                    \
    }                                                               \
} while(0)


float rand_float();
int vasprintf(char **strp, const char *format, va_list ap);
int asprintf(char **strp, const char *format, ...);
uint32_t next_power_2(uint32_t n);
int *range_random_sample(int n, int size);
int *sel_k_min(int *arr, int n, int size);

#endif