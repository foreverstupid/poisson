#include "operator.h"

/*
 * Help function returns X part of the Laplas operator that is defined
 * by the equation in the point (xi, yj).
 */
inline static scalar_t laplas_xpart(
    int i,
    int j,
    const Operator *op,
    const Matrix *u)
{
    return 
        1.0 / op->h1 *
        (
            at(op->A, i + 1, j) / op->h1 *
            (
                at(u, i + 1, j) - at(u, i, j)
            ) -
            at(op->A, i, j) / op->h1 *
            (
                at(u, i, j) - at(u, i - 1, j)
            )
        );
}

/*
 * Help function returns Y part of the Laplas operator that is defined
 * by the equation in the point (xi, yj).
 */
inline static scalar_t laplas_ypart(
    int i,
    int j,
    const Operator *op,
    const Matrix *u)
{
    return 
        1.0 / op->h2 *
        (
            at(op->B, i, j + 1) / op->h2 *
            (
                at(u, i, j + 1) - at(u, i, j)
            ) -
            at(op->B, i, j) / op->h2 *
            (
                at(u, i, j) - at(u, i, j - 1)
            )
        );
}

/*
 * Applies operator action on the left boundary, storing result in the
 * first operand: v = Au.
 */
static void left_boundary_perform(
    Matrix *v,
    const Operator *op,
    const Matrix *u)
{
    int j;
    scalar_t tmp;
    scalar_t alpha_part;

    if (op->left_type == first)
    {
        for (j = 0; j < op->PhiL->ny; j++)
        {
            set(v, 0, j, op->PhiL->data[j]);
        }
    }
    else
    {
        alpha_part = op->left_type == second ? 0.0 : 2.0 / op->h1;
        for (j = 1; j < op->PhiL->ny - 1; j++)
        {
            tmp =
                -2.0 / op->h1 *
                (
                    at(op->A, 1, j) * (at(u, 1, j) - at(u, 0, j))
                ) +
                (at(op->Q, 0, j) + alpha_part) * at(u, 0, j) -
                laplas_ypart(0, j, op, u);

            set(v, 0, j, tmp);
        }
    }
}

/*
 * Applies operator action on the right boundary, storing result in the
 * first operand: v = Au.
 */
static void right_boundary_perform(
    Matrix *v,
    const Operator *op,
    const Matrix *u)
{
    int j;
    scalar_t tmp;
    scalar_t alpha_part;
    int M = v->nx - 1;

    if (op->right_type == first)
    {
        for (j = 0; j < op->PhiL->ny; j++)
        {
            set(v, M, j, op->PhiL->data[j]);
        }
    }
    else
    {
        alpha_part = op->left_type == second ? 0.0 : 2.0 / op->h1;
        for (j = 1; j < op->PhiL->ny - 1; j++)
        {
            tmp =
                2.0 / op->h1 *
                (
                    at(op->A, M, j) * (at(u, M, j) - at(u, M - 1, j))
                ) +
                (at(op->Q, M, j) + alpha_part) * at(u, M, j) -
                laplas_ypart(M, j, op, u);

            set(v, M, j, tmp);
        }
    }
}

/*
 * Applies operator action on the bottom boundary, storing result in the
 * first operand: v = Au.
 */
static void bottom_boundary_perform(
    Matrix *v,
    const Operator *op,
    const Matrix *u)
{
    int i;
    scalar_t tmp;
    scalar_t alpha_part;

    if (op->left_type == first)
    {
        for (i = 0; i < op->PhiL->nx; i++)
        {
            set(v, i, 0, op->PhiL->data[i]);
        }
    }
    else
    {
        alpha_part = op->left_type == second ? 0.0 : 2.0 / op->h2;
        for (i = 1; i < op->PhiL->nx - 1; i++)
        {
            tmp =
                -2.0 / op->h2 *
                (
                    at(op->B, i, 1) * (at(u, i, 1) - at(u, i, 0))
                ) +
                (at(op->Q, i, 0) + alpha_part) * at(u, i, 0) -
                laplas_xpart(i, 0, op, u);

            set(v, i, 0, tmp);
        }
    }
}

/*
 * Applies operator action on the top boundary, storing result in the
 * first operand: v = Au.
 */
static void top_boundary_perform(
    Matrix *v,
    const Operator *op,
    const Matrix *u)
{
    int i;
    scalar_t tmp;
    scalar_t alpha_part;
    int N = u->ny - 1;

    if (op->left_type == first)
    {
        for (i = 0; i < op->PhiL->nx; i++)
        {
            set(v, i, N, op->PhiL->data[i]);
        }
    }
    else
    {
        alpha_part = op->left_type == second ? 0.0 : 2.0 / op->h2;
        for (i = 1; i < op->PhiL->nx - 1; i++)
        {
            tmp =
                2.0 / op->h2 *
                (
                    at(op->B, i, N) * (at(u, i, N) - at(u, i, N - 1))
                ) +
                (at(op->Q, i, N) + alpha_part) * at(u, i, N) -
                laplas_xpart(i, N, op, u);

            set(v, i, N, tmp);
        }
    }
}

void apply(Matrix *v, const Operator *op, const Matrix *u)
{
    int i;
    int j;
    scalar_t tmp;

    /* calculate inner points */
    for (i = 1; i < u->nx - 1; i++)
    {
        for (j = 1; j < u->ny - 1; j++)
        {
            tmp =
                -laplas_xpart(i, j, op, u) - laplas_ypart(i, j, op, u) +
                at(op->Q, i, j) * at(u, i, j);
            
            set(v, i, j, tmp);
        }
    }

    left_boundary_perform(v, op, u);
    right_boundary_perform(v, op, u);
    bottom_boundary_perform(v, op, u);
    top_boundary_perform(v, op, u);
}