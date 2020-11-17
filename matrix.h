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
    /*
     * Matrix values in the row-major order.
     */
    scalar_t *data;

    /*
     * The size of the matrix along X-axis (count of columns).
     */
    int nx;

    /*
     * The size of the matrix along Y-axis (count of rows).
     */
    int ny;
} Matrix;



/*
 * The mask defines which elements of the matrix should be used in the
 * operations.
 */
typedef struct MatrixMask
{
    /*
     * Defines the first using element along X-axis.
     */
    int x0;

    /*
     * Defines the last using element along X-axis.
     */
    int x1;

    /*
     * Defines the first using element along Y-axis.
     */
    int y0;

    /*
     * Defines the last using element along Y-axis.
     */
    int y1;
} MatrixMask;



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
    const MatrixMask *mask,
    scalar_t h1,
    scalar_t h2);

/*
 * Returns the square of the norm of the given matrix.
 */
inline static scalar_t get_squared_norm(
    const Matrix *matrix,
    const MatrixMask *mask,
    scalar_t h1,
    scalar_t h2)
{
    return dot_product(matrix, matrix, mask, h1, h2);
}

/*
 * Multiplies the given matrix by the given scalar componentwise.
 */
void multiply(Matrix *matrix, scalar_t f, const MatrixMask *mask);

/*
 * Subtracts one matrix from another, storing the result in the first
 * operand.
 * Note: the size of matricies should be the same.
 */
void sub(Matrix *m1, const Matrix *m2, const MatrixMask *mask);

/*
 * Calculate linear combination of the given matricies, storing result
 * into the first operand: res = m1 + t*m2
 * Note: the size of matricies should be the same.
 */
void linear_combination(
    Matrix *res,
    const Matrix *m1,
    scalar_t t,
    const Matrix *m2,
    const MatrixMask *mask);

/*
 * Returns C-norm of the matrix.
 */
scalar_t get_cnorm(const Matrix *m, const MatrixMask *mask);

/*
 * Gets the row of the matrix as an array.
 */
void get_row(const Matrix *m, int i, scalar_t *row);

/*
 * Sets the row of the matrix from the array.
 * Note: an array should contain an appropriate count of elements.
 */
void set_row(Matrix *m, int i, scalar_t *row);

/*
 * Gets the column of the matrix as an array.
 */
void get_column(const Matrix *m, int i, scalar_t *column);

/*
 * Sets the column of the matrix from the array.
 * Note: an array should contain an appropriate count of elements.
 */
void set_column(Matrix *m, int i, scalar_t *column);

#endif