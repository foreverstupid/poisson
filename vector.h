#ifndef HPC_POISSON_VECTOR_MODULE_H
#define HPC_POISSON_VECTOR_MODULE_H

#include <stdlib.h>
#include <math.h>

/*
 * The type of the scalar.
 */
typedef double scalar_t;

/*
 * Vector of scalars.
 */
typedef struct Vector
{
    scalar_t *data;
    int length;
} Vector;

/*
 * Scalar function of a single argument.
 */
typedef scalar_t (*Func)(scalar_t);



/*
 * Creates a new vector of the given length.
 */
Vector *new_vector(int length);

/*
 * Creates a new vector that contains samples of the given scalar
 * function.
 * Note: samples include start and end points.
 */
Vector *get_vector_from_func(
    Func func,
    scalar_t start,
    scalar_t end,
    int length);

/*
 * Returns the copy of the given vector.
 */
Vector *copy(Vector *other);

/*
 * Disposes all resources that are taken by the given vector.
 */
void delete_vector(Vector *vector);

/*
 * Returns the dot product of two vectors of the same length.
 */
scalar_t dot_product(Vector *v1, Vector *v2);

/*
 * Returns the norm of the given vector.
 */
scalar_t get_norm(Vector *vector);

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

#endif