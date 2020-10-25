#include <stdio.h>
#include "definitions.h"
#include "solve.h"

scalar_t F(scalar_t x, scalar_t y)
{
    return sin(x * y);
}

scalar_t phi(scalar_t x, scalar_t y)
{
    return cos(2 * M_PI * (x + y));
}

scalar_t k(scalar_t x, scalar_t y)
{
    return x;
}

scalar_t q(scalar_t x, scalar_t y)
{
    return x + y;
}


int main(int argc, char **argv)
{
    Problem problem = {
        .funcs = {
            .boundary_type = first,
            .phi = &phi,
            .f = &F,
            .k = &k,
            .q = &q
        },
        .ranges = {
            .xrange = {
                .count = 100,
                .start = -2,
                .end = 2
            },
            .yrange = {
                .count = 100,
                .start = -2,
                .end = 2
            }
        }
    };

    Matrix *solution = solve(&problem, 1e-6);
    delete_matrix(solution);

    return 0;
}
