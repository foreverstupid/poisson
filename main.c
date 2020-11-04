#include <stdio.h>
#include <math.h>
//#include <mpi.h>
#include "definitions.h"
#include "solve.h"
#include "output.h"

#define SOLUTION_FILE_NAME "solution.csv"

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

    SolvingInfo config;
    char solution_file[80];

    sscanf(argv[1], "%d", &(config.x_grid_count));
    sscanf(argv[2], "%d", &(config.y_grid_count));
    sscanf(argv[3], S_FORMAT, &(config.eps));
    sscanf(argv[4], "%d", &(config.iteration_print_frequency));
    config.output_dir = argv[5];

    //MPI_Init(&argc, &argv);
    Matrix *solution = solve(&problem, &config);
    //MPI_Finalize();
    sprintf(solution_file, "%s/" SOLUTION_FILE_NAME, argv[5]);
    write_as_csv(solution, solution_file);
    delete_matrix(solution);

    return 0;
}
