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
            at(op->A, i + 1, j) *
                (at(u, i + 1, j) - at(u, i, j)) -
            at(op->A, i, j) *
                (at(u, i, j) - at(u, i - 1, j))
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

    if (op->left != first)
    {
        alpha_part = op->left == second ? 0.0 : 2.0 / op->h1;
        for (j = 1; j < v->ny - 1; j++)
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

    if (op->right != first)
    {
        alpha_part = op->right == second ? 0.0 : 2.0 / op->h1;
        for (j = 1; j < v->ny - 1; j++)
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

    if (op->bottom != first)
    {
        alpha_part = op->bottom == second ? 0.0 : 2.0 / op->h2;
        for (i = 1; i < v->nx - 1; i++)
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

    if (op->top != first)
    {
        alpha_part = op->top == second ? 0.0 : 2.0 / op->h2;
        for (i = 1; i < v->nx - 1; i++)
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

    if (op->bottom != first && op->left != first)
    {
        a1_part = op->left == second ? 0.0 : 2.0 / op->h1;
        a2_part = op->bottom == second ? 0.0 : 2.0 / op->h2;

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

    if (op->bottom != first && op->right != first)
    {
        a1_part = op->right == second ? 0.0 : 2.0 / op->h1;
        a2_part = op->bottom == second ? 0.0 : 2.0 / op->h2;

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

    if (op->top != first && op->right != first)
    {
        a1_part = op->right == second ? 0.0 : 2.0 / op->h1;
        a2_part = op->top == second ? 0.0 : 2.0 / op->h2;

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

    if (op->top != first && op->left != first)
    {
        a1_part = op->left == second ? 0.0 : 2.0 / op->h1;
        a2_part = op->top == second ? 0.0 : 2.0 / op->h2;

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



void apply_first_order_boundary(Matrix *v, const Operator *op)
{
    int i;
    int M = v->nx - 1;
    int N = v->ny - 1;

    if (op->left == first)
    {
        for (i = 0; i <= N; i++)
        {
            set(v, 0, i, at(op->F, 0, i));
        }
    }

    if (op->right == first)
    {
        for (i = 0; i <= N; i++)
        {
            set(v, M, i, at(op->F, M, i));
        }
    }

    if (op->bottom == first)
    {
        for (i = 0; i <= M; i++)
        {
            set(v, i, 0, at(op->F, i, 0));
        }
    }

    if (op->top == first)
    {
        for (i = 0; i <= M; i++)
        {
            set(v, i, N, at(op->F, i, N));
        }
    }
}



/*
 * Initializes boundary conditions.
 */
static void init_boundary(Operator *op, const Problem *problem)
{
    int i;
    int M = op->F->nx - 1;
    int N = op->F->ny - 1;
    scalar_t t;
    scalar_t tmp;

    op->left = problem->boundary.left.type;
    t = problem->area.y1;
    for (i = 0; i <= N; i++)
    {
        tmp = problem->boundary.left.phi(t);
        t += op->h2;

        if (op->left != first)
        {
            tmp = at(op->F, 0, i) + 2.0 / op->h1 * tmp;
        }

        set(op->F, 0, i, tmp);
    }

    op->right = problem->boundary.right.type;
    t = problem->area.y1;
    for (i = 0; i <= N; i++)
    {
        tmp = problem->boundary.right.phi(t);
        t += op->h2;

        if (op->right != first)
        {
            tmp = at(op->F, M, i) + 2.0 / op->h1 * tmp;
        }

        set(op->F, M, i, tmp);
    }

    op->bottom = problem->boundary.bottom.type;
    t = problem->area.x1;
    for (i = 0; i <= M; i++)
    {
        tmp = problem->boundary.bottom.phi(t);
        t += op->h1;

        if (op->bottom != first)
        {
            tmp = at(op->F, i, 0) + 2.0 / op->h2 * tmp;
        }

        set(op->F, i, 0, tmp);
    }

    op->top = problem->boundary.top.type;
    t = problem->area.x1;
    for (i = 0; i <= M; i++)
    {
        tmp = problem->boundary.top.phi(t);
        t += op->h1;

        if (op->top != first)
        {
            tmp = at(op->F, i, N) + 2.0 / op->h2 * tmp;
        }

        set(op->F, i, N, tmp);
    }
}



Operator *new_operator(
    const Problem *problem,
    const NumericalConfig *config)
{
    Range xrange = {
        .count = config->x_grid_count,
        .start = problem->area.x1,
        .end = problem->area.x2
    };

    Range yrange = {
        .count = config->y_grid_count,
        .start = problem->area.y1,
        .end = problem->area.y2
    };

    Range2 ranges = { .xrange = xrange, .yrange = yrange };

    scalar_t h1 = get_range_step(&(ranges.xrange));
    scalar_t h2 = get_range_step(&(ranges.yrange));

    Range2 A_ranges = {
        .xrange = {
            .count = config->x_grid_count,
            .start = problem->area.x1 - h1 / 2,
            .end = problem->area.x2 - h1 / 2
        },
        .yrange = yrange
    };

    Range2 B_ranges = {
        .xrange = xrange,
        .yrange = {
            .count = config->y_grid_count,
            .start = problem->area.y1 - h2 / 2,
            .end = problem->area.y2 - h2 / 2
        }
    };

    Operator *op = (Operator *)malloc(sizeof(Operator));

    op->A = get_matrix_from_func2(problem->funcs.k, &(A_ranges));
    op->B = get_matrix_from_func2(problem->funcs.k, &(B_ranges));
    op->Q = get_matrix_from_func2(problem->funcs.q, &(ranges));
    op->F = get_matrix_from_func2(problem->funcs.f, &(ranges));

    op->h1 = h1;
    op->h2 = h2;

    init_boundary(op, problem);
    return op;
}



void delete_operator(Operator *op)
{
    delete_matrix(op->A);
    delete_matrix(op->B);
    delete_matrix(op->Q);
    delete_matrix(op->F);

    free(op);
}