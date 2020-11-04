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
        (
            at(op->A, i + 1, j) * (at(u, i + 1, j) - at(u, i, j)) -
            at(op->A, i, j) * (at(u, i, j) - at(u, i - 1, j))
        ) /
        (op->h1 * op->h1);
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
        (
            at(op->B, i, j + 1) *
                (at(u, i, j + 1) - at(u, i, j)) -
            at(op->B, i, j) *
                (at(u, i, j) - at(u, i, j - 1))
        ) /
        (op->h2 * op->h2);
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
                -2.0 * at(op->A, 1, j) /
                    (op->h1 * op->h1) *
                    (at(u, 1, j) - at(u, 0, j)) +
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
        for (j = 0; j < op->PhiR->ny; j++)
        {
            set(v, M, j, op->PhiR->data[j]);
        }
    }
    else
    {
        alpha_part = op->right_type == second ? 0.0 : 2.0 / op->h1;
        for (j = 1; j < op->PhiR->ny - 1; j++)
        {
            tmp =
                2.0 * at(op->A, M, j) /
                    (op->h1 * op->h1) *
                    (at(u, M, j) - at(u, M - 1, j)) +
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

    if (op->bottom_type == first)
    {
        for (i = 0; i < op->PhiB->nx; i++)
        {
            set(v, i, 0, op->PhiB->data[i]);
        }
    }
    else
    {
        alpha_part = op->bottom_type == second ? 0.0 : 2.0 / op->h2;
        for (i = 1; i < op->PhiB->nx - 1; i++)
        {
            tmp =
                -2.0 * at(op->B, i, 1) /
                    (op->h2 * op->h2) *
                    (at(u, i, 1) - at(u, i, 0)) +
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

    if (op->top_type == first)
    {
        for (i = 0; i < op->PhiT->nx; i++)
        {
            set(v, i, N, op->PhiT->data[i]);
        }
    }
    else
    {
        alpha_part = op->top_type == second ? 0.0 : 2.0 / op->h2;
        for (i = 1; i < op->PhiT->nx - 1; i++)
        {
            tmp =
                2.0 * at(op->B, i, N) /
                    (op->h2 * op->h2) *
                    (at(u, i, N) - at(u, i, N - 1)) +
                (at(op->Q, i, N) + alpha_part) * at(u, i, N) -
                laplas_xpart(i, N, op, u);

            set(v, i, N, tmp);
        }
    }
}



/*
 * Applies operator action on the left-bottom corner of the area,
 * storing result in the first operand: v = Au.
 */
static void lb_corner_perform(
    Matrix *v,
    const Operator *op,
    const Matrix *u)
{
    scalar_t tmp;
    scalar_t a1_part;
    scalar_t a2_part;

    if (op->bottom_type != first && op->left_type != first)
    {
        a1_part = op->left_type == second ? 0.0 : 2.0 / op->h1;
        a2_part = op->bottom_type == second ? 0.0 : 2.0 / op->h2;

        tmp =
            -2.0 * at(op->A, 1, 0) /
                (op->h1 * op->h1) *
                (at(u, 1, 0) - at(u, 0, 0)) -
            2.0 * at(op->B, 0, 1) /
                (op->h2 * op->h2) *
                (at(u, 0, 1) - at(u, 0, 0)) +
            at(u, 0, 0) * (at(op->Q, 0, 0) + a1_part + a2_part);
        
        set(v, 0, 0, tmp);
    }
}



/*
 * Applies operator action on the right-bottom corner of the area,
 * storing result in the first operand: v = Au.
 */
static void rb_corner_perform(
    Matrix *v,
    const Operator *op,
    const Matrix *u)
{
    scalar_t tmp;
    scalar_t a1_part;
    scalar_t a2_part;
    int M = u->nx - 1;

    if (op->bottom_type != first && op->right_type != first)
    {
        a1_part = op->right_type == second ? 0.0 : 2.0 / op->h1;
        a2_part = op->bottom_type == second ? 0.0 : 2.0 / op->h2;

        tmp =
            2.0 * at(op->A, M, 0) /
                (op->h1 * op->h1) *
                (at(u, M, 0) - at(u, M - 1, 0)) -
            2.0 * at(op->B, M, 1) /
                (op->h2 * op->h2) *
                (at(u, M, 1) - at(u, M, 0)) +
            at(u, M, 0) * (at(op->Q, M, 0) + a1_part + a2_part);
        
        set(v, M, 0, tmp);
    }
}



/*
 * Applies operator action on the right-top corner of the area,
 * storing result in the first operand: v = Au.
 */
static void rt_corner_perform(
    Matrix *v,
    const Operator *op,
    const Matrix *u)
{
    scalar_t tmp;
    scalar_t a1_part;
    scalar_t a2_part;
    int M = u->nx - 1;
    int N = u->ny - 1;

    if (op->top_type != first && op->right_type != first)
    {
        a1_part = op->right_type == second ? 0.0 : 2.0 / op->h1;
        a2_part = op->top_type == second ? 0.0 : 2.0 / op->h2;

        tmp =
            2.0 * at(op->A, M, N) /
                (op->h1 * op->h1) *
                (at(u, M, N) - at(u, M - 1, N)) +
            2.0 * at(op->B, M, N) /
                (op->h2 * op->h2) *
                (at(u, M, N) - at(u, M, N - 1)) +
            at(u, M, N) * (at(op->Q, M, N) + a1_part + a2_part);
        
        set(v, M, N, tmp);
    }
}



/*
 * Applies operator action on the left-top corner of the area,
 * storing result in the first operand: v = Au.
 */
static void lt_corner_perform(
    Matrix *v,
    const Operator *op,
    const Matrix *u)
{
    scalar_t tmp;
    scalar_t a1_part;
    scalar_t a2_part;
    int N = u->ny - 1;

    if (op->top_type != first && op->left_type != first)
    {
        a1_part = op->left_type == second ? 0.0 : 2.0 / op->h1;
        a2_part = op->top_type == second ? 0.0 : 2.0 / op->h2;

        tmp =
            -2.0 * at(op->A, 1, N) /
                (op->h1 * op->h1) *
                (at(u, 1, N) - at(u, 0, N)) +
            2.0 * at(op->B, 0, N) /
                (op->h2 * op->h2) *
                (at(u, 0, N) - at(u, 0, N - 1)) +
            at(u, 0, N) * (at(op->Q, 0, N) + a1_part + a2_part);
        
        set(v, 0, N, tmp);
    }
}



void apply(Matrix *v, const Operator *op, const Matrix *u)
{
    int i;
    int j;
    scalar_t tmp;

    /* calculate inner points */
    for (j = 1; j < u->ny - 1; j++)
    {
        for (i = 1; i < u->nx - 1; i++)
        {
            tmp =
                -laplas_xpart(i, j, op, u) -
                laplas_ypart(i, j, op, u) +
                at(op->Q, i, j) * at(u, i, j);
            
            set(v, i, j, tmp);
        }
    }

    /* calculate boundaries */
    left_boundary_perform(v, op, u);
    right_boundary_perform(v, op, u);
    bottom_boundary_perform(v, op, u);
    top_boundary_perform(v, op, u);

    /* calculate corners */
    lb_corner_perform(v, op, u);
    rb_corner_perform(v, op, u);
    rt_corner_perform(v, op, u);
    lt_corner_perform(v, op, u);
}