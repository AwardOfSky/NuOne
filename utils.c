#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"


float rand_float() {
    return (float)rand() / (float)RAND_MAX;
}


int vasprintf(char **strp, const char *format, va_list ap) {
    int len = _vscprintf(format, ap);
    if (len == -1)
        return -1;
    char *str = (char*)malloc((size_t) len + 1);
    if (!str)
        return -1;
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