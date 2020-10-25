#ifndef HPC_POISSON_DEFINITIONS_MODULE_H
#define HPC_POISSON_DEFINITIONS_MODULE_H

/*
 * Almost zero value.
 */
#define EPS 1e-10

/*
 * The type of the scalar.
 */
typedef double scalar_t;

/*
 * Linear range of scalar value samples.
 * Note: start and end values are parts of the range.
 */
typedef struct Range
{
    scalar_t start;
    scalar_t end;
    int count;
} Range;

/*
 * Gets the step of the range.
 */
inline scalar_t get_step(const Range *range)
{
    return (range->end - range->start) / (range->count - 1);
}

/*
 * Scalar function of a single argument.
 */
typedef scalar_t (*Func)(scalar_t);

/*
 * Scalar function of two arguments.
 */
typedef scalar_t (*Func2)(scalar_t, scalar_t);

/*
 * The type of the boundary conditions.
 */
typedef enum BoundaryType
{
    /*
     * u = phi
     */
    first,

    /*
     * ku' = phi
     */
    second,

    /*
     * ku' + u = phi
     */
    third
} BoundaryType;

/*
 * Contains input data for the Poisson boundary value problem.
 */
typedef struct Problem
{
    /*
     * Potential function.
     */
    Func2 q;

    /*
     * Coefficient function of the Laplass operator.
     */
    Func2 k;

    /*
     * Right-part function.
     */
    Func2 f;

    /*
     * Boundary values function.
     */
    Func2 phi;

    /*
     * The type of the boundary conditions.
     */
    BoundaryType boundary_type;

    /*
     * The range of the X-values.
     */
    Range x_range;

    /*
     * The range of the Y-values.
     */
    Range y_range;
} Problem;

#endif