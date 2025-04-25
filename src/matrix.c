/*
 * matrix.c
 * Matrix type for use in polyfit calculation
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

#include "matrix.h"

#include <stdlib.h> // malloc(), free()
#include <string.h> // memset()
#ifndef NDEBUG
#include <stdio.h>  // printf()
#endif


#define MATRIX_ALLOC(nrows, ncols)                          \
        const unsigned int nr = nrows;                      \
        const unsigned int nc = ncols;                      \
        size_t data_len = sizeof(double) * nr * nc;         \
        void *buf = malloc(sizeof(matrix_t) + data_len);    \
        if (buf == NULL) return NULL;                       \
        matrix_t *matrix = (matrix_t *) buf;                \
        double *data = (double *) (buf + sizeof(matrix_t)); \
        matrix->rows = nr;                                  \
        matrix->cols = nc;                                  \
        matrix->data = data


matrix_t *matrix_raw(unsigned int rows, unsigned int cols) {
    MATRIX_ALLOC(rows, cols);
    return matrix;
}

matrix_t *matrix_create(unsigned int rows, unsigned int cols) {
    MATRIX_ALLOC(rows, cols);
    memset(data, 0, data_len);
    return matrix;
}

matrix_t *matrix_transpose(const matrix_t *src) {
    MATRIX_ALLOC(src->cols, src->rows);
    for (int r=0; r < nr; r++) {
        for (int c=0; c < nc; c++) {
            matrix_set(matrix, r, c, matrix_get(src, c, r));
        }
    }
    return matrix;
}

matrix_t *matrix_product(const matrix_t *m1, const matrix_t *m2) {
#ifndef NDEBUG
    if (m1 == NULL || m2 == NULL || m1->cols != m2->rows) {
        fprintf(stderr, "illegal state in matrix_product\n");
        return NULL;
    }
#endif
    MATRIX_ALLOC(m1->rows, m2->cols);
    memset(data, 0, data_len);
    for (int i=0; i < nr; i++) {
        for (int j=0; j < nc; j++) {
            for (int k=0; k < m1->cols; k++) {
                *matrix_get(matrix, i, j) += matrix_get_val(m1, i, k) * matrix_get_val(m2, k, j);
            }
        }
    }
    return matrix;
}

void matrix_destroy(matrix_t *matrix) {
    free(matrix);
}

#ifndef NDEBUG
void matrix_dbg(const matrix_t *matrix) {
    for (int r=0; r < matrix->rows; r++) {
        for (int c=0; c < matrix->cols; c++) {
            printf("   %f", matrix_get_val(matrix, r, c));
        }
        printf("\n");
    }
}
#endif
