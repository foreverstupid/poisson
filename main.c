#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>

#include "definitions.h"
#include "solve.h"
#include "output.h"

#define X1 0.0
#define X2 4.0
#define Y1 0.0
#define Y2 3.0

scalar_t k(scalar_t x, scalar_t y)
{
    return x + y + 4.0;
}

scalar_t q(scalar_t x, scalar_t y)
{
    return x + y;
}

scalar_t F(scalar_t x, scalar_t y)
{
    scalar_t sq = sqrt(4.0 + x * y);
    return
        q(x, y) * sq -
        0.5 * (x + y) / sq +
        0.25 * (x * x + y * y) * (4 + x + y) / (sq * sq * sq);
}

scalar_t left(scalar_t t)
{
    scalar_t sq = sqrt(4 + X1 * t);
    return -0.5 * t * k(X1, t) / sq + sq;
}

scalar_t right(scalar_t t)
{
    scalar_t sq = sqrt(4 + X2 * t);
    return 0.5 * t * k(X2, t) / sq + sq;
}

scalar_t bottom(scalar_t t)
{
    scalar_t sq = sqrt(4 + t * Y1);
    return -0.5 * t * k(t, Y1) / sq;
}

scalar_t top(scalar_t t)
{
    scalar_t sq = sqrt(4 + t * Y2);
    return 0.5 * t * k(t, Y2) / sq;
}



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

    mpi->x_proc_count = dims[0];
    mpi->y_proc_count = dims[1];

    MPI_Comm_rank(mpi->grid_comm, &rank);
    MPI_Cart_coords(mpi->grid_comm, rank, 2, dims);

    mpi->x_proc_idx = dims[0];
    mpi->y_proc_idx = dims[1];
}



int main(int argc, char **argv)
{ 
    Problem problem = {
        .funcs = {
            .f = &F,
            .k = &k,
            .q = &q
        },
        .boundary = {
            .left = {
                .phi = left,
                .type = third
            },
            .right = {
                .phi = right,
                .type = third
            },
            .top = {
                .phi = top,
                .type = second
            },
            .bottom = {
                .phi = bottom,
                .type = second
            }
        },
        .area = {
            .x1 = X1,
            .x2 = X2,
            .y1 = Y1,
            .y2 = Y2
        }
    };

    SolvingConfig config;
    InitResult init_res;

    sscanf(argv[1], "%d", &(config.num.x_grid_count));
    sscanf(argv[2], "%d", &(config.num.y_grid_count));
    sscanf(argv[3], SF, &(config.num.eps));
    sscanf(argv[4], "%d", &(config.log.iteration_print_frequency));
    config.log.write_matrix = write_matrix;
    config.log.log_message = log_info;

    printf("Prepare to init MPI...\n");
    if (MPI_Init(&argc, &argv) != 0)
    {
        fprintf(stderr, "Cannot init MPI");
        return 1;
    }

    set_mpi_config(&(config.mpi));
    printf("Output module initialization...\n");
    init_res = init_output(
        argv[5],
        config.mpi.x_proc_idx,
        config.mpi.y_proc_idx);

    if (init_res != success)
    {
        MPI_Finalize();
        fprintf(stderr, "Can not init output module: %d\n", init_res);
        return 1;
    }

    solve(&problem, &config);
    dispose_output();
    MPI_Finalize();

    return 0;
}
