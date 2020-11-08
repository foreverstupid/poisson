#ifndef HPC_POISSON_SOLVE_MODULE_H
#define HPC_POISSON_SOLVE_MODULE_H

#include <mpi.h>
#include <stdio.h>
#include "definitions.h"
#include "matrix.h"
#include "operator.h"
#include "process.h"
#include "output.h"

/*
 * Solves the given Poisson boundary values problem,
 * storing the result according to the given log configuration.
 */
void solve(const Problem *problem, const SolvingConfig *config);

#endif