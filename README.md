# polyfit
C11/C++ library for least-squares polynomial regression.
Aims to have a pleasant, portable API
with plenty of fluff.

## Comparison
|                                                     | [WasabiThumb/polyfit](https://github.com/WasabiThumb/polyfit) | [natedomin/polyfit](https://github.com/natedomin/polyfit) | [henryfo/polyfit](https://github.com/henryfo/polyfit) |
| --------------------------------------------------: | :-----------------------------------------------------------: | :-------------------------------------------------------: | :---------------------------------------------------: |
| **Test coverage**                                   | ✅                                                             | ✅                                                         | ✅                                                    |
| **Unbounded**                                       | ✅                                                             | ❌<sup>1</sup>                                             | ✅                                                    |
| [**Stringification**](#stringification)             | ✅                                                             | ❌                                                         | ✅                                                    |
| **Memory safety**                                   | ✅                                                             | ✅                                                         | ❌<sup>2</sup>                                        |
| [**Error handling**](#error-handling)               | ✅                                                             | ✅<sup>3</sup>                                             | ✅<sup>3</sup> <sup>4</sup>                           |
| [**Point set abstraction**](#point-set-abstraction) | ✅                                                             | ❌                                                         | ❌                                                    |
| [**Evaluation**](#evaluation)                       | ✅                                                             | ❌                                                         | ❌<sup>5</sup>                                        |
| **Non-allocating**                                  | ❌                                                             | ✅                                                         | ❌                                                    |

1. maxOrder can be increased, but bounded by the stack
2. Buffers are [not properly freed when errors arise](https://github.com/henryfo/polyfit/blob/a4e4ab87857279eeb4af3678ebcb3bfbe09c4569/src/polyfit.c#L152)
3. Lacks strerror
4. Elides multiple error states
5. See [open issue](https://github.com/natedomin/polyfit/issues/1)

## API
### Regression
##### pf_err_t [pf_fit](https://github.com/WasabiThumb/polyfit/blob/8f04f3024ce8cfa7216a225d8cc5e7ee141e1d4d/include/polyfit.h#L98)(const pf_points_t *points, unsigned int order, double *coefficients)
> Fits an N-order polynomial to a set of input points (see [point set abstraction](#point-set-abstraction)).

##### pf_err_t [pf_fit_easy](https://github.com/WasabiThumb/polyfit/blob/8f04f3024ce8cfa7216a225d8cc5e7ee141e1d4d/include/polyfit.h#L113)(unsigned int count, const double *xs, const double *ys, unsigned int, double *coefficients)
> Fits an N-order polynomial to a set of input points.

### Error Handling
##### const char *[pf_strerror](https://github.com/WasabiThumb/polyfit/blob/8f04f3024ce8cfa7216a225d8cc5e7ee141e1d4d/include/polyfit.h#L55)(pf_err_t code)
> Returns a human-readable string describing a polyfit status code.

### Stringification
##### size_t [pf_strpoly](https://github.com/WasabiThumb/polyfit/blob/8f04f3024ce8cfa7216a225d8cc5e7ee141e1d4d/include/polyfit.h#L130)(pf_len_t order, const double *coefficients, char *buf, size_t len)
> Writes into "buf" a string representing a polynomial.

### Evaluation
##### double [pf_eval](https://github.com/WasabiThumb/polyfit/blob/8f04f3024ce8cfa7216a225d8cc5e7ee141e1d4d/include/polyfit.h#L142)(pf_len_t order, const double *coefficients, double x);
> Computes the value of a polynomial expression for a given x

## Point set abstraction
While optimal for many use cases, the ``const double *xs, const double *ys`` signature
may cause unecessary allocations and/or ``memcpy`` calls. Consider the case where you have
an array of vector structs:
```c
typedef struct {
    double x;
    double y;
} v2;

v2 arr[N_POINTS];
```
This array is fundamentally incompatible with the aforementioned signature. Instead of mangling this
data into the format demanded by the signature, we can instead use the ``pf_points_t`` abstraction:
```c
struct pf_points_by_v2s {
    pf_points_t super;
    const v2 *v2s;
};

static double pf_points_by_v2s_x(const pf_points_t *self, pf_len_t idx) {
    return ((struct pf_points_by_v2s *) self)->v2s[idx].x;
}

static double pf_points_by_v2s_y(const pf_points_t *self, pf_len_t idx) {
    return ((struct pf_points_by_v2s *) self)->v2s[idx].y;
}

//

struct pf_points_by_v2s points_impl;
points_impl.super.count = N_POINTS;
points_impl.super.x = pf_points_by_v2s_x;
points_impl.super.y = pf_points_by_v2s_y;
points_impl.v2s = arr;

pf_points_t *points = &points_impl.super;
// pass "points" to pf_fit
```
This involves 0 heap allocations, and is a *reasonable* tradeoff given that each point
is only read approximately once by the internal algorithm. This also means that the x/y accessor may be fairly
expensive while incurring no additional penalty.
