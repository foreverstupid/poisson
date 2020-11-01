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

scalar_t u(scalar_t x, scalar_t y)
{
    return sqrt(4 + x * y);
}

scalar_t k(scalar_t x, scalar_t y)
{
    return 4 + x + y;
}

scalar_t q(scalar_t x, scalar_t y)
{
    scalar_t t = x + y;
    return t > 0 ? t : 0.0;
}

scalar_t F(scalar_t x, scalar_t y)
{
    scalar_t uval = u(x, y);
    scalar_t u3 = uval * uval * uval;

    return
        q(x, y) * uval -
        0.5 * (x + y) / uval +
        0.25 * (x * x + y * y) * (4 + x + y) / u3;
}

scalar_t left(scalar_t t)
{
    return u(X1, t);
}

scalar_t right(scalar_t t)
{
    return u(X2, t);
}

scalar_t bottom(scalar_t t)
{
    return u(t, Y1);
}

scalar_t top(scalar_t t)
{
    return u(t, Y2);
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

    sscanf(argv[1], "%d", &(config.x_grid_count));
    sscanf(argv[2], "%d", &(config.y_grid_count));
    sscanf(argv[3], "%lf", &(config.eps));

    //MPI_Init(&argc, &argv);
    Matrix *solution = solve(&problem, &config);
    //MPI_Finalize();
    write_as_csv(solution, argv[4]);
    delete_matrix(solution);

    return 0;
}
