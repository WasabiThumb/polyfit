/*
 * matrix.h
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

#pragma once

/**
 * A matrix of double-precision floats,
 * stored in row-major order
 */
typedef struct matrix {
    unsigned int rows;
    unsigned int cols;
    double      *data;
} matrix_t;

/**
 * Allocates a matrix and leaves its contents uninitialized
 * @param rows Row count, not explicitly validated
 * @param cols Column count, not explicitly validated
 * @return The matrix, or NULL if failed to allocate
 */
matrix_t *matrix_raw(unsigned int rows, unsigned int cols);

/**
 * Allocates a matrix and initializes its contents to 0
 * @param rows Row count, not explicitly validated
 * @param cols Column count, not explicitly validated
 * @return The matrix, or NULL if failed to allocate
 */
matrix_t *matrix_create(unsigned int rows, unsigned int cols);

/**
 * Allocates a new matrix which represents a transposed form of the given matrix
 * @param src The matrix to transpose
 * @return A new matrix, or NULL if failed to allocate
 */
matrix_t *matrix_transpose(const matrix_t *src);

/**
 * Allocates a new matrix which represents the product of the given matrices
 * @param m1 The first matrix
 * @param m2 The second matrix
 * @return A new matrix, or NULL if failed to allocate
 */
matrix_t *matrix_product(const matrix_t *m1, const matrix_t *m2);

/**
 * Destroys a matrix previously created by matrix_create
 */
void matrix_destroy(matrix_t *matrix);

/**
 * Returns a reference to a value from a matrix
 * @param row The row index of the value
 * @param col The column index of the value
 * @return A pointer to the value within the matrix data
 */
static inline double *matrix_get(const matrix_t *matrix, unsigned int row, unsigned int col) {
    return &matrix->data[(row * matrix->cols) + col];
}

/**
 * Returns a value from a matrix
 * @param row The row index of the value
 * @param col The column index of the value
 * @return The value at the given row and column index within the matrix
 */
static inline double matrix_get_val(const matrix_t *matrix, unsigned int row, unsigned int col) {
    return *matrix_get(matrix, row, col);
}

/**
 * Sets a value in a matrix
 * @param row The row index of the value
 * @param col The column index of the value
 * @param val A pointer to the value to set
 */
static inline void matrix_set(matrix_t *matrix, unsigned int row, unsigned int col, const double *val) {
    *matrix_get(matrix, row, col) = *val;
}

/**
 * Sets a value in a matrix
 * @param row The row index of the value
 * @param col The column index of the value
 * @param val The value to set
 */
static inline void matrix_set_val(matrix_t *matrix, unsigned int row, unsigned int col, double val) {
    *matrix_get(matrix, row, col) = val;
}

#ifdef NDEBUG
/**
 * Prints information about the given matrix to stdout
 * only if debugging is enabled.
 */
static inline void matrix_dbg(const matrix_t *matrix) { }
#else
/**
 * Prints information about the given matrix to stdout
 * only if debugging is enabled.
 */
void matrix_dbg(const matrix_t *matrix);
#endif
