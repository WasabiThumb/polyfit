#include "polyfit.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define POINT_COUNT 8
#define TARGET_FN(x) cos((4.0 / M_PI) * (x))

int run(pf_len_t order, double *xs, double *ys) {
    printf("- Order: %u\n", order);

    double *coefficients = malloc(sizeof(double) * order);
    if (coefficients == NULL) {
        fprintf(stderr, "failed to allocate coefficients buffer\n");
        return 1;
    }

    pf_err_t err = pf_fit_easy(POINT_COUNT, xs, ys, order, coefficients);
    if (err != PF_OK) {
        fprintf(stderr, "fit failed: %s\n", pf_strerror(err));
        free(coefficients);
        return 1;
    }

    size_t str_len = pf_strpoly(order, coefficients, NULL, 0);
    char *str = malloc(str_len);
    if (str == NULL) {
        free(coefficients);
        return 1;
    }
    pf_strpoly(order, coefficients, str, str_len);
    printf("  f(x) = %s\n", str);
    free(str);

    double margin = 0;
    for (pf_len_t x=0; x < POINT_COUNT; x++) {
        margin += fabs(ys[x] - pf_eval(order, coefficients, (double) x));
    }
    printf("  error = %f\n\n", margin);

    free(coefficients);
    return 0;
}

int main() {
    double xs[POINT_COUNT];
    double ys[POINT_COUNT];

    double dx;
    for (pf_len_t i=0; i < POINT_COUNT; i++) {
        xs[i] = dx = (double) i;
        ys[i] = TARGET_FN(dx);
    }

    for (pf_len_t i=1; i <= POINT_COUNT; i++) {
        if (run(i, xs, ys) == 1) return 1;
    }
    return 0;
}

