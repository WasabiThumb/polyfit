/*
 * polyfit.h
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

#include <stdlib.h>
#include <stdint.h>

#ifndef POLYFIT_H
#define POLYFIT_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Status code for polyfit
 * @see PF_OK
 * @see PF_EALLOC
 * @see PF_EPARAM
 * @see PF_ESOLVE
 */
typedef uint_fast8_t pf_err_t;

/** Success */
#define PF_OK     ((pf_err_t) 0)

/** Allocation failure */
#define PF_EALLOC ((pf_err_t) 1)

/** Bad parameter */
#define PF_EPARAM ((pf_err_t) 2)

/** Unable to solve */
#define PF_ESOLVE ((pf_err_t) 3)

/**
 * Returns a human-readable string describing a
 * polyfit status code.
 */
const char *pf_strerror(pf_err_t code);

//

/**
 * Type for polyfit indices & array lengths
 */
typedef unsigned int pf_len_t;

/**
 * Abstraction for a sequence of 2D points.
 * The "x" component is the independent variable.
 * The "y" component is the dependent variable.
 */
typedef struct pf_points {
    /** Number of points contained in this sequence */
    pf_len_t count;

    /**
     * Accessor for the X coordinates of points in this sequence
     * @param self Reference to self
     * @param idx Point index
     */
    double (*x)(const struct pf_points *self, pf_len_t idx);

    /**
     * Accessor for the Y coordinates of points in this sequence
     * @param self Reference to self
     * @param idx Point index
     */
    double (*y)(const struct pf_points *self, pf_len_t idx);
} pf_points_t;

//

/**
 * Fits an N-order polynomial to a set of input points
 * @param points Input points. Must have a length greater than or equal to "order".
 * @param order Order of the polynomial to compute
 * @param coefficients Buffer to receive the polynomial coefficients
 * @return PF_OK (0) on success
 * @see pf_fit_easy
 */
pf_err_t pf_fit(const pf_points_t *points, pf_len_t order, double *coefficients);

/**
 * Fits an N-order polynomial to a set of input points
 * @param count Number of input points. May not be less than "order".
 * @param xs X coordinates of input points,
 *           length equal to sizeof(int) * count
 * @param ys Y coordinates of input points,
 *           length equal to sizeof(int) * count
 * @param order Order of the polynomial to compute
 * @param coefficients Buffer to receive the polynomial coefficients, from highest to lowest order,
 *                     length equal to sizeof(double) * order
 * @return PF_OK (0) on success
 * @see pf_fit
 */
pf_err_t pf_fit_easy(pf_len_t count, const double *xs, const double *ys, pf_len_t order, double *coefficients);

//

/**
 * Writes into "buf" a string representing a polynomial.
 * @param order The order of the polynomial
 * @param coefficients The coefficients of the polynomial, from highest to lowest order,
 *                     length equal to sizeof(double) * order
 * @param buf The buffer to receive the string
 * @param len The length of the buffer. If 0, "buf" may be NULL.
 * @return If successful, returns the length in chars of the string written to "buf",
 *         excluding the terminating null character.
 *         If "buf" is too small to contain the string, returns the length in chars
 *         of a buffer that is sufficiently large to contain the string including the
 *         terminating null character.
 */
size_t pf_strpoly(pf_len_t order, const double *coefficients, char *buf, size_t len);

//

/**
 * Computes the value of a polynomial expression for a given x
 * @param order The order of the polynomial
 * @param coefficients The coefficients of the polynomial, from highest to lowest order,
 *                     length equal to sizeof(double) * order
 * @param x The value of x
 * @return The value of the polynomial expression, f(x)
 */
double pf_eval(pf_len_t order, const double *coefficients, double x);

#ifdef __cplusplus
}
#endif

#endif //POLYFIT_H
