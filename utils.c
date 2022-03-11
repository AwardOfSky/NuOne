#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "utils.h"


float rand_float() {
    return (float)rand() / (float)RAND_MAX;
}


int vasprintf(char **strp, const char *format, va_list ap) {
    int len = _vscprintf(format, ap);
    if (len == -1) return -1;
    char *str = (char*)malloc((size_t) len + 1);
    if (!str) return -1;
    int retval = vsnprintf(str, len + 1, format, ap);
    if (retval == -1) {
        free(str);
        return -1;
    }
    *strp = str;
    return retval;
}


int asprintf(char **strp, const char *format, ...) {
    va_list ap;
    va_start(ap, format);
    int retval = vasprintf(strp, format, ap);
    va_end(ap);
    return retval;
}


uint32_t next_power_2(uint32_t n) {
    --n;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    return ++n;
}


int *range_random_sample(int n, int size) {
    if(size < n && size > 0) {
        int sampling_size;
        int *result = (int *)malloc(sizeof(int) * size);
        char *helper;
        
        if(n < STACK_ALLOC) {
            helper = (char *)alloca(n * sizeof(char));
            memset(helper, 0, n);
        } else {
            helper = (char *)calloc(n, sizeof(char));
        }

        // when the next voodoo magic is true, exclude selected items instead of including
        if(size <= n * 0.555) {
            sampling_size = size;

            for(int i = 0; i < sampling_size; ++i) {
                int index = rand() % n;
                while(helper[index] == 1) {
                    index = rand() % n;
                }
                helper[index] = 1;
                result[i] = index;
            }
        } else {
            sampling_size = n - size;

            for(int i = 0; i < sampling_size; ++i) {
                int index = rand() % n;
                while(helper[index] == 1) {
                    index = rand() % n;
                }
                helper[index] = 1;
            }

            int running = 0;
            for(int i = 0; i < size; ++i) {
                while(helper[running] == 1) {
                    ++running;
                }
                result[i] = running++;
            }            
        }

        if(n >= STACK_ALLOC) {
            free(helper);
        }

        return result; 
    } else { // this error should not be possible
        printf(PREERR"Wrong number of elements passed to sampling function!\n");
        return NULL;
    }
}


int *sel_k_min(int *arr, int n, int size) {
    int max_of_mins = INT32_MAX;
    int *k = (int *)malloc((size + 1) * sizeof(*k));

    k[0] = arr[0];
    for(int i = 1; i < n; ++i) {
        int k_els = min(i, size);
        //printf("\n");
        //PRINT_ARR(k, k_els, "%d");
        int cur_val = arr[i];
        
        if (cur_val < max_of_mins) {
            
            int j;
            if (cur_val > k[0] && cur_val < k[k_els - 1]) {
                //printf("Adding middle!\n");
                if (k_els > BINSEARCH_TRESHOLD) {
                    int inc = k_els >> 1;
                    j = inc;
                    while (k[j] < cur_val || k[j - 1] > cur_val) {
                        inc = (inc > 1) ? (inc >> 1) : 1;
                        j += (cur_val < k[j]) ? -inc : inc;
                    }
                } else {
                    for(j = 0; j < k_els && k[j] < cur_val; ++j) {}
                }
            } else {

                if (cur_val <= k[0]) {
                    //printf("Adding beggining!\n");
                    j = 0;
                } else {
                    //printf("Adding end!\n");
                    j = k_els;
                }
            }

            //printf("(j, k_els, cur_val, (size - 1 - j)), (%d %d %d %d)\n", j, k_els, cur_val, (size - 1 - j));
            memmove(k + j + 1, k + j, (size - j) * sizeof(*k));
            k[j] = arr[i];

            if (k_els == size) {
                max_of_mins = k[size - 1]; 
            }
        }
    }

    return k;
}

//MAKE_COMPUTE_STD(int); // compute_std_int
MAKE_COMPUTE_STD(float); // compute_std_float
//MAKE_COMPUTE_MEAN(int); // compute_mean_int
//MAKE_COMPUTE_MEAN(float); // compute_mean_float
MAKE_REALLOC(float);
MAKE_REALLOC(char);
