#ifndef HPC_POISSON_DEFINITIONS_MODULE_H
#define HPC_POISSON_DEFINITIONS_MODULE_H

/*
 * Almost zero value.
 */
#define EPS 1e-12

/*
 * Format of the scalar value.
 */
#define S_FORMAT "%lf"

/*
 * The type of the scalar.
 */
typedef double scalar_t;

/*
 * Scalar function of a single argument.
 */
typedef scalar_t (*Func)(scalar_t);

/*
 * Scalar function of two arguments.
 */
typedef scalar_t (*Func2)(scalar_t, scalar_t);



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
 * XY-ranges of scalar value samples.
 */
typedef struct Range2
{
    Range xrange;
    Range yrange;
} Range2;

/*
 * Gets the step of the range.
 */
static inline scalar_t get_step(const Range *range)
{
    return (range->end - range->start) / (range->count - 1);
}



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
 * Contains the information about the boundary condition.
 */
typedef struct BoundaryCondition
{
    Func phi;
    BoundaryType type;
} BoundaryCondition;

/*
 * Contains information about known functions that are parts of the
 * problem.
 */
typedef struct FunctionsInfo
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
} FunctionsInfo;

/*
 * Information about the equation area.
 */
typedef struct Area
{
    scalar_t x1;
    scalar_t x2;
    scalar_t y1;
    scalar_t y2;
} Area;

/*
 * Contains full info about the boundary conditions.
 * Note: boundary consitions are suppsoed to be consistent, i.e.
 *      phi_L(y1) = phi_B(x1)
 *      phi_L(y2) = phi_T(x1)
 *      phi_R(y1) = phi_B(x2)
 *      phi_R(y2) = phi_T(x2)
 */
typedef struct BoundaryInfo
{
    BoundaryCondition left;
    BoundaryCondition right;
    BoundaryCondition top;
    BoundaryCondition bottom;
} BoundaryInfo;

/*
 * Contains input data for the Poisson boundary value problem.
 */
typedef struct Problem
{
    /*
     * Known functions of the equation.
     */
    FunctionsInfo funcs;

    /*
     * Boundary conditions.
     */
    BoundaryInfo boundary;

    /*
     * The equation area.
     */
    Area area;
} Problem;



/*
 * Contains solving configuration info.
 */
typedef struct SolvingInfo
{
    int x_grid_count;
    int y_grid_count;

    /*
     * If C-norm of difference between iterations is less than this
     * value, the solving process will stop.
     */
    scalar_t eps;

    /*
     * The solving process will print iterations with this frequency
     * (i.e. if it is 10 then every 10-th iteration will be printed).
     * If it equals to -1 then none iteration will be printed.
     */
    int iteration_print_frequency;

    /*
     * The directory that will contain output information.
     */
    const char *output_dir;
} SolvingInfo;

#endif