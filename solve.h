#ifndef HPC_POISSON_SOLVE_MODULE_H
#define HPC_POISSON_SOLVE_MODULE_H

#include "definitions.h"
#include "matrix.h"

/*
 * Solves the given Poisson boundary values problem,
 * returning the found function samples.
 * The <eps> parameter defines the exit
 * consition: when iterations are close to each other
 * in C-norm by this value, the iterative process will
 * be completed.
 */
Matrix *solve(const Problem *problem, scalar_t eps);

#endif