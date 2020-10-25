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
    int ncols;          /* count of columns */
    int nrows;          /* count of rows */
} Matrix;



/*
 * Creates a new matrix of the given size.
 */
Matrix *new_matrix(int nrows, int ncols);

/*
 * Creates a new matrix that contains samples of the given scalar
 * function of two arguments.
 */
Matrix *get_matrix_from_func2(Func2 func, const Range2 *ranges);

/*
 * Returns the copy of the given matrix.
 */
Matrix *copy(const Matrix *other);

/*
 * Disposes all resources that are taken by the given matrix.
 */
void delete_matrix(Matrix *matrix);

/*
 * Returns the dot product of two matricies of the same size.
 */
scalar_t dot_product(const Matrix *m1, const Matrix *m2);

/*
 * Returns the norm of the given matrix.
 */
inline scalar_t get_norm(const Matrix *matrix)
{
    return sqrt(dot_product(matrix, matrix));
}

/*
 * Multiplies the given matrix by the given scalar componentwise.
 */
void multiply(Matrix *matrix, scalar_t f);

/*
 * Adds one matrix to another, storing the result in the first operand.
 * Note: the size of matricies should be the same.
 */
void add(Matrix *m1, const Matrix *m2);

/*
 * Subtracts one matrix from another, storing the result in the first
 * operand.
 * Note: the size of matricies should be the same.
 */
void sub(Matrix *m1, const Matrix *m2);

/*
 * Returns C-norm of the matrixs difference.
 * Note: the size of matricies should be the same.
 */
scalar_t get_difference_cnorm(const Matrix *m1, const Matrix *m2);

#endif