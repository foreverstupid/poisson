#ifndef HPC_POISSON_MATRIX_MODULE_H
#define HPC_POISSON_MATRIX_MODULE_H

#include <stdlib.h>
#include <math.h>
#include "definitions.h"

/*
 * Matrix of scalars.
 */
typedef struct Matrix
{
    scalar_t *data;     /* matrix values in a row-major order */
    int nx;             /* size along X-axis (count of columns) */
    int ny;             /* size along Y-axis (count of rows) */
} Matrix;



/*
 * Creates a new matrix of the given size.
 */
Matrix *new_matrix(int nx, int ny);

/*
 * Creates a new matrix that contains samples of the given scalar
 * function of two arguments.
 */
Matrix *get_matrix_from_func2(Func2 func, const Range2 *ranges);

/*
 * Returns the copy of the given matrix.
 */
Matrix *copy_matrix(const Matrix *other);

/*
 * Disposes all resources that are taken by the given matrix.
 */
void delete_matrix(Matrix *matrix);

/*
 * Help function for getting an element corresponding to xi and yj, i.e.
 * m[j, i].
 */
static inline scalar_t at(const Matrix *m, int i, int j)
{
    return m->data[j * m->nx + i];
}

/*
 * Help fucntion for setting an element corresponding to xi and yj, i.e.
 * m[j, i].
 */
static inline void set(Matrix *m, int i, int j, scalar_t val)
{
    m->data[j * m->nx + i] = val;
}

/*
 * Returns the dot product of two matricies of the same size.
 */
scalar_t dot_product(
    const Matrix *m1,
    const Matrix *m2,
    scalar_t h1,
    scalar_t h2);

/*
 * Returns the square of the norm of the given matrix.
 */
inline static scalar_t get_squared_norm(const Matrix *matrix, scalar_t h1, scalar_t h2)
{
    return dot_product(matrix, matrix, h1, h2);
}

/*
 * Multiplies the given matrix by the given scalar componentwise.
 */
void multiply(Matrix *matrix, scalar_t f);

/*
 * Subtracts one matrix from another, storing the result in the first
 * operand.
 * Note: the size of matricies should be the same.
 */
void sub(Matrix *m1, const Matrix *m2);

/*
 * Calculate linear combination of the given matricies, storing result
 * into the first operand: res = m1 + t*m2
 * Note: the size of matricies should be the same.
 */
void linear_combination(
    Matrix *res,
    const Matrix *m1,
    scalar_t t,
    const Matrix *m2);

/*
 * Returns C-norm of the matrixs difference.
 * Note: the size of matricies should be the same.
 */
scalar_t get_difference_cnorm(const Matrix *m1, const Matrix *m2);

#endif