#ifndef HPC_POISSON_SOLVE_MODULE_H
#define HPC_POISSON_SOLVE_MODULE_H

//#include <mpi.h>
#include "definitions.h"
#include "matrix.h"
#include "operator.h"

/*
 * Solves the given Poisson boundary values problem,
 * returning the found function samples.
 */
Matrix *solve(const Problem *problem, const SolvingInfo *config);

#endif