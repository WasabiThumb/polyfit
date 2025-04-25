/*
 * polyfit.c
 * Polynomial regression library
 *
 * Copyright 2025 Wasabi Codes
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "polyfit.h"
#include "matrix.h"

#include <stdio.h> // snprintf
#include <math.h>  // pow


const char *pf_strerror(pf_err_t code) {
    switch (code) {
        case PF_OK:
            return "OK";
        case PF_EALLOC:
            return "Allocation failure";
        case PF_EPARAM:
            return "Bad parameter";
        case PF_ESOLVE:
            return "Unable to solve";
        default:
            return "Unknown error";
    }
}

//

struct pf_points_by_array {
    __attribute__((unused)) pf_points_t super;
    const double *xs;
    const double *ys;
};

/** @internal */
static double pf_points_array_x(const pf_points_t *self, pf_len_t idx) {
    return ((struct pf_points_by_array *) self)->xs[idx];
}

/** @internal */
static double pf_points_array_y(const pf_points_t *self, pf_len_t idx) {
    return ((struct pf_points_by_array *) self)->ys[idx];
}

//

pf_err_t pf_fit(const pf_points_t *points, pf_len_t order, double *coefficients) {
    const pf_len_t degree = order - 1;
    double tx, ty;
    pf_err_t err;

    // Check parameters
    if (points == NULL || coefficients == NULL) {
        return PF_EPARAM;
    }
    pf_len_t point_count = points->count;
    if (point_count < order) {
        return PF_EPARAM;
    }

    // Create matrix A
    matrix_t *mat_a = matrix_raw(point_count, order);
    if (mat_a == NULL) {
        return PF_EALLOC;
    }
    for (pf_len_t r=0; r < point_count; r++) {
        for (pf_len_t c=0; c < order; c++) {
            tx = points->x(points, r);
            matrix_set_val(mat_a, r, c, pow(tx, (double) (degree - c)));
        }
    }

    // Create matrix B
    matrix_t *mat_b = matrix_create(point_count, 1);
    if (mat_b == NULL) {
        err = PF_EALLOC;
        goto ex_1;
    }
    for (pf_len_t r=0; r < point_count; r++) {
        ty = points->y(points, r);
        matrix_set_val(mat_b, r, 0, ty);
    }

    // Create matrix AT (transpose A)
    matrix_t *mat_at = matrix_transpose(mat_a);
    if (mat_at == NULL) {
        err = PF_EALLOC;
        goto ex_2;
    }

    // Create matrix ATA (multiply AT & A)
    matrix_t *mat_ata = matrix_product(mat_at, mat_a);
    if (mat_ata == NULL) {
        err = PF_EALLOC;
        goto ex_3;
    }

    // Create matrix ATB (multiply AT & B)
    matrix_t *mat_atb = matrix_product(mat_at, mat_b);
    if (mat_atb == NULL) {
        err = PF_EALLOC;
        goto ex_4;
    }

    // Solve the system of linear equations
    // (AT)Ax = (AT)b
    for (pf_len_t c=0; c < mat_ata->cols; c++) {
        pf_len_t pr = c; // pivot row
        double prv = matrix_get_val(mat_ata, pr, c);
        if (prv == 0.0) {
            err = PF_ESOLVE;
            goto ex_5;
        }
        for (pf_len_t r=0; r < mat_ata->rows; r++) {
            if (r == pr) continue;
            double trv = matrix_get_val(mat_ata, r, c);
            double factor = trv / prv;
            for (pf_len_t c2=0; c2 < mat_ata->cols; c2++) {
                *matrix_get(mat_ata, r, c2) -= matrix_get_val(mat_ata, pr, c2) * factor;
            }
            *matrix_get(mat_atb, r, 0) -= matrix_get_val(mat_atb, pr, 0) * factor;
        }
    }
    for (pf_len_t c=0; c < mat_ata->cols; c++) {
        pf_len_t pr = c; // pivot row
        double prv = matrix_get_val(mat_ata, pr, c);
        *matrix_get(mat_ata, pr, c) /= prv;
        *matrix_get(mat_atb, pr, 0) /= prv;
    }

    for (pf_len_t i=0; i < order; i++) {
        coefficients[i] = matrix_get_val(mat_atb, i, 0);
    }

    err = PF_OK;
    ex_5: matrix_destroy(mat_atb);
    ex_4: matrix_destroy(mat_ata);
    ex_3: matrix_destroy(mat_at);
    ex_2: matrix_destroy(mat_b);
    ex_1: matrix_destroy(mat_a);
    return err;
}

pf_err_t pf_fit_easy(pf_len_t count, const double *xs, const double *ys, pf_len_t order, double *coefficients) {
    struct pf_points_by_array points = {
            .super = {
                    .count = count,
                    .x     = pf_points_array_x,
                    .y     = pf_points_array_y
            },
            .xs = xs,
            .ys = ys
    };
    return pf_fit(&points.super, order, coefficients);
}

//

size_t pf_strpoly(pf_len_t order, const double *coefficients, char *buf, size_t len) {
    if (order == 0) return 0;

    size_t head = 0;
    int ok = 1;
    pf_len_t exp = order;
    double co;
    int wr;

    for (pf_len_t i=0; i < order; i++) {
        exp--;
        co = coefficients[i];

        if (i != 0) {
            char op;
            if (co < 0.0f) {
                op = '-';
                co = -co;
            } else {
                op = '+';
            }
            wr = snprintf(buf + head, len, " %c ", op);
            head += wr;
            if (wr < len) {
                len -= wr;
            } else {
                len = 0;
                ok = 0;
            }
        }

        if (exp == 0) {
            wr = snprintf(buf + head, len, "%f", co);
        } else if (exp == 1) {
            wr = snprintf(buf + head, len, "%fx", co);
        } else {
            wr = snprintf(buf + head, len, "%fx^%d", co, exp);
        }
        head += wr;
        if (wr < len) {
            len -= wr;
        } else {
            len = 0;
            ok = 0;
        }
    }

    return ok ? head : head + 1;
}

//

double pf_eval(pf_len_t order, const double *coefficients, double x) {
    if (order == 0) return 0;
    double n = coefficients[order - 1];
    if (order == 1) return n;
    double q = 1;
    pf_len_t i = order - 2;
    while (1) {
        q *= x;
        n += q * coefficients[i];
        if (i == 0) break;
        i--;
    }
    return n;
}
