#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "definitions.h"
#include "solve.h"
#include "output.h"
#include "problem.h"



/*
 * Fills the MPI information.
 */
void set_mpi_config(MpiConfig *mpi)
{
    int proc_count;
    int rank;
    int dims[] = { 0, 0 };
    int periods[] = { 0, 0 };

    MPI_Comm_size(MPI_COMM_WORLD, &proc_count);
    MPI_Dims_create(proc_count, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &(mpi->grid_comm));

    if (mpi->grid_comm == MPI_COMM_NULL)
    {
        return;
    }

    mpi->x_proc_count = dims[0];
    mpi->y_proc_count = dims[1];

    MPI_Comm_rank(mpi->grid_comm, &rank);
    MPI_Cart_coords(mpi->grid_comm, rank, 2, dims);

    mpi->x_proc_idx = dims[0];
    mpi->y_proc_idx = dims[1];
}



/*
 * Measures the program running time.
 */
void measure_running_time(double start, int is_main)
{
    double local_duration;
    double global_duration;

    local_duration = MPI_Wtime() - start;
    MPI_Reduce(
        &local_duration, &global_duration,
        1, MPI_DOUBLE, MPI_MAX,
        0, MPI_COMM_WORLD);

    if (is_main)
    {
        printf("Running time: %lfs\n", global_duration);
    }
}



/*
 * Prints the run information.
 */
void print_run_info(char **argv)
{
    printf(
        "Running configuration:\n"
        "    grid size:          %s x %s\n"
        "    required precision: %s\n"
        "    log frequency:      %s\n"
        "    data dir:           %s\n\n",
        argv[1], argv[2],
        argv[3],
        argv[4],
        argv[5]);
}



/*
 * Performs all needed initializations, returning the rank of the process
 * in the MPI_COMM_WORLD communicator. Returns -1 in the case of errors.
 */
int init(int argc, char **argv, SolvingConfig *config)
{
    OutputInitResult init_res;
    int rank;

    if (argc < 5)
    {
        fprintf(stderr, "Too few cmd args\n");
        return -1;
    }

    if (MPI_Init(&argc, &argv) != 0)
    {
        fprintf(stderr, "Cannot init MPI");
        return -1;
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
    {
        print_run_info(argv);
    }

    sscanf(argv[1], "%d", &(config->num.x_grid_count));
    sscanf(argv[2], "%d", &(config->num.y_grid_count));
    sscanf(argv[3], SF, &(config->num.eps));
    sscanf(argv[4], "%d", &(config->log.iteration_print_frequency));

    set_mpi_config(&(config->mpi));
    if (config->mpi.grid_comm == MPI_COMM_NULL)
    {
        /* not an error so return rank, but do nothing */
        printf("Process with rank %d is extra\n", rank);
        return rank;
    }

    printf("Output module initialization...\n");
    init_res = init_output(
        argv[5],
        config->mpi.x_proc_idx,
        config->mpi.y_proc_idx);

    config->log.write_matrix = write_matrix;
    config->log.log_message = log_info;

    if (init_res != success)
    {
        fprintf(stderr, "Can not init output module: %d\n", init_res);
        return -1;
    }

    return rank;
}



int main(int argc, char **argv)
{
    int rank;
    SolvingConfig config;

    rank = init(argc, argv, &config);
    if (rank == -1)
    {
        MPI_Finalize();
        return 1;
    }

    if (config.mpi.grid_comm == MPI_COMM_NULL)
    {
        MPI_Finalize();
        return 0;
    }

    double start = MPI_Wtime();
    solve(&main_problem, &config);
    measure_running_time(start, rank == 0);

    dispose_output();
    MPI_Finalize();

    return 0;
}
