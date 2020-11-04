#include <stdio.h>
#include <math.h>
//#include <mpi.h>
#include "definitions.h"
#include "solve.h"
#include "output.h"

#define X1 0.0
#define X2 4.0
#define Y1 0.0
#define Y2 3.0

scalar_t k(scalar_t x, scalar_t y)
{
    return 1.0;
}

scalar_t q(scalar_t x, scalar_t y)
{
    return 0.0;
}

scalar_t F(scalar_t x, scalar_t y)
{
    return -sin(x);
}

scalar_t left(scalar_t t)
{
    return sin(X1);
}

scalar_t right(scalar_t t)
{
    return sin(X2);
}

scalar_t bottom(scalar_t t)
{
    return sin(t);
}

scalar_t top(scalar_t t)
{
    return sin(t);
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
                .type = first
            },
            .right = {
                .phi = right,
                .type = first
            },
            .top = {
                .phi = top,
                .type = first
            },
            .bottom = {
                .phi = bottom,
                .type = first
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

    sscanf(argv[1], "%d", &(config.x_grid_count));
    sscanf(argv[2], "%d", &(config.y_grid_count));
    sscanf(argv[3], S_FORMAT, &(config.eps));

    //MPI_Init(&argc, &argv);
    Matrix *solution = solve(&problem, &config);
    //MPI_Finalize();
    write_as_csv(solution, argv[4]);
    delete_matrix(solution);

    return 0;
}
