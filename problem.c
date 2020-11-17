#include "problem.h"

static scalar_t k(scalar_t x, scalar_t y)
{
    return x + y + 4.0;
}

static scalar_t q(scalar_t x, scalar_t y)
{
    return x + y;
}

static scalar_t F(scalar_t x, scalar_t y)
{
    scalar_t sq = sqrt(4.0 + x * y);
    return
        q(x, y) * sq -
        0.5 * (x + y) / sq +
        0.25 * (x * x + y * y) * (4 + x + y) / (sq * sq * sq);
}

static scalar_t left(scalar_t t)
{
    scalar_t sq = sqrt(4 + X1 * t);
    return -0.5 * t * k(X1, t) / sq + sq;
}

static scalar_t right(scalar_t t)
{
    scalar_t sq = sqrt(4 + X2 * t);
    return 0.5 * t * k(X2, t) / sq + sq;
}

static scalar_t bottom(scalar_t t)
{
    scalar_t sq = sqrt(4 + t * Y1);
    return -0.5 * t * k(t, Y1) / sq;
}

static scalar_t top(scalar_t t)
{
    scalar_t sq = sqrt(4 + t * Y2);
    return 0.5 * t * k(t, Y2) / sq;
}



/*
 * The global problem for solving.
 */
const Problem main_problem = {
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