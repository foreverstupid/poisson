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
        .boundary_type = first,
        .phi = &phi,
        .f = &F,
        .k = &k,
        .q = &q,
        .x_range = {
            .count = 100,
            .start = -2,
            .end = 2
        },
        .y_range = {
            .count = 100,
            .start = -2,
            .end = 2
        }
    };

    Vector *solution = solve(&problem, 1e-6);
    delete_vector(solution);

    return 0;
}
