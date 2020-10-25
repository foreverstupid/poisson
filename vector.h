#ifndef HPC_POISSON_VECTOR_MODULE_H
#define HPC_POISSON_VECTOR_MODULE_H

#include <stdlib.h>
#include <math.h>
#include "definitions.h"

/*
 * Vector of scalars.
 */
typedef struct Vector
{
    scalar_t *data;
    int length;
} Vector;



/*
 * Creates a new vector of the given length.
 */
Vector *new_vector(int length);

/*
 * Creates a new vector that contains samples of the given scalar
 * function.
 */
Vector *get_vector_from_func(Func func, const Range *range);

/*
 * Returns the copy of the given vector.
 */
Vector *copy(const Vector *other);

/*
 * Disposes all resources that are taken by the given vector.
 */
void delete_vector(Vector *vector);

/*
 * Returns the dot product of two vectors of the same length.
 */
scalar_t dot_product(const Vector *v1, const Vector *v2);

/*
 * Returns the norm of the given vector.
 */
inline scalar_t get_norm(const Vector *vector)
{
    return sqrt(dot_product(vector, vector));
}

/*
 * Multiplies the given vector by the given scalar componentwise.
 */
void multiply(Vector *vector, scalar_t f);

/*
 * Adds one vector to another, storing the result in the first operand.
 * Note: the length of vectors should be the same.
 */
void add(Vector *v1, Vector *v2);

/*
 * Subtracts one vector from another, storing the result in the first
 * operand.
 * Note: the length of vectors should be the same.
 */
void sub(Vector *v1, Vector *v2);

/*
 * Returns C-norm of the vectors difference.
 * Note: the length of vectors should be the same.
 */
scalar_t get_difference_norm(const Vector *v1, const Vector *v2);

#endif