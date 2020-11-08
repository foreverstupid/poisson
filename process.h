#ifndef HPC_POISSON_PROCESS_MODULE_H
#define HPC_POISSON_PROCESS_MODULE_H

#include <mpi.h>
#include <stdio.h>
#include "definitions.h"
#include "operator.h"


/*
 * Configuration of the process, solves the subproblem.
 */
typedef struct ProcessInfo
{
    /*
     * The cart communicator.
     */
    MPI_Comm comm;

    /*
     * The rank of the current process in the cart communicator.
     */
    int rank;

    /*
     * The rank of the left process neighbour in the cart communicator.
     * MPI_PROC_NULL if the current process has no left neighbour.
     */
    int left_neighbour;

    /*
     * The rank of the right process neighbour in the cart communicator.
     * MPI_PROC_NULL if the current process has no right neighbour.
     */
    int right_neighbour;

    /*
     * The rank of the bottom process neighbour in the cart communicator.
     * MPI_PROC_NULL if the current process has no bottom neighbour.
     */
    int bottom_neighbour;

    /*
     * The rank of the top process neighbour in the cart communicator.
     * MPI_PROC_NULL if the current process has no top neighbour.
     */
    int top_neighbour;

    /*
     * Operator of the subproblem.
     */
    Operator *op;

    /*
     * Log configuration. It is defined only for the main process.
     */
    LogConfig log;

    /*
     * The mask for matrix operation that prevents changing values of
     * the boundary condition of the first order.
     */
    MatrixMask mask;

    /*
     * Difference between iterations that indicates that the process of
     * solving is over.
     */
    scalar_t eps;
} ProcessInfo;



/*
 * Initializes process information.
 */
ProcessInfo *get_process_info(
    const Problem *global_problem,
    const SolvingConfig *config);

#endif