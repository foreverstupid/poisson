#include "solve.h"

/*
 * Performs a single iteration of the method of the smallest residuals
 * (u -> v).
 */
static void perform_iteration(
    Matrix *v,
    const Operator *op,
    const Matrix *u)
{
    scalar_t tau;           /* iterational parameter */
    apply(v, op, u);
    sub(v, op->F);          /* residual */

    Matrix *Ar = copy_matrix(v);    /* operator applied to residual */
    apply(Ar, op, v);

    tau =
        dot_product(Ar, v, op->h1, op->h2) /
        get_squared_norm(Ar, op->h1, op->h2);

    multiply(v, tau);
    linear_combination(v, u, -tau, v);
    delete_matrix(Ar);
}

/*
 * Performs iteration process until reaching the given error
 * between iterations.
 */
static Matrix *find_solution(const Operator *op, scalar_t eps)
{
    Matrix *u = copy_matrix(op->F);
    Matrix *v = new_matrix(u->nx, u->ny);
    Matrix *tmp;
    scalar_t curr_eps;

    do
    {
        perform_iteration(v, op, u);
        curr_eps = get_difference_cnorm(u, v);

        tmp = u;
        u = v;
        v = tmp;
    }
    while (curr_eps > eps);

    delete_matrix(v);
    return u;
}



/*
 * Performs correction of the right part in the corners of the area
 * accorfing to the boundary conditions.
 */
static void corners_correction(const Operator *op)
{
    int M = op->F->nx - 1;
    int N = op->F->ny - 1;
    scalar_t mul = 2.0 / op->h1 + 2.0 / op->h2;
    scalar_t tmp;

    if (op->left_type != first && op->bottom_type != first)
    {
        tmp = at(op->F, 0, 0) + mul * at(op->PhiL, 0, 0);
        set(op->F, 0, 0, tmp);
    }

    if (op->right_type != first && op->bottom_type != first)
    {
        tmp = at(op->F, M, 0) + mul * at(op->PhiR, 0, 0);
        set(op->F, M, 0, tmp);
    }

    if (op->right_type != first && op->top_type != first)
    {
        tmp = at(op->F, M, N) + mul * at(op->PhiR, 0, N);
        set(op->F, M, N, tmp);
    }

    if (op->left_type != first && op->top_type != first)
    {
        tmp = at(op->F, 0, N) + mul * at(op->PhiL, 0, N);
        set(op->F, 0, N, tmp);
    }
}



/*
 * Performs correction of the right part according to the
 * boundary conditions.
 */
static void correct_right_part(const Operator *op)
{
    int i;
    int j;
    int M = op->F->nx - 1;
    int N = op->F->ny - 1;
    scalar_t tmp;

    if (op->left_type != first)
    {
        for (j = 1; j < N; j++)
        {
            tmp = 2.0 / op->h2 * at(op->PhiL, 0, j);
            set(op->F, 0, j, at(op->F, 0, j) + tmp);
        }
    }

    if (op->right_type != first)
    {
        for (j = 1; j < N; j++)
        {
            tmp = 2.0 / op->h2 * at(op->PhiR, M, j);
            set(op->F, M, j, at(op->F, M, j) + tmp);
        }
    }

    if (op->bottom_type != first)
    {
        for (i = 1; i < M; i++)
        {
            tmp = 2.0 / op->h1 * at(op->PhiB, i, 0);
            set(op->F, i, 0, at(op->F, i, 0) + tmp);
        }
    }

    if (op->top_type != first)
    {
        for (i = 1; i < M; i++)
        {
            tmp = 2.0 / op->h1 * at(op->PhiT, i, 0);
            set(op->F, i, N, at(op->F, i, N) + tmp);
        }
    }

    corners_correction(op);
}

/*
 * Creates an operator, generated byt he equation.
 */
static Operator *create_operator(
    const Problem *problem,
    const SolvingInfo *config)
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

    scalar_t h1 = get_step(&(ranges.xrange));
    scalar_t h2 = get_step(&(ranges.yrange));

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

    op->PhiL = get_column_from_func(problem->boundary.left.phi, &(yrange));
    op->left_type = problem->boundary.left.type;

    op->PhiR = get_column_from_func(problem->boundary.right.phi, &(yrange));
    op->right_type = problem->boundary.right.type;

    op->PhiB = get_row_from_func(problem->boundary.bottom.phi, &(xrange));
    op->bottom_type = problem->boundary.bottom.type;

    op->PhiT = get_row_from_func(problem->boundary.top.phi, &(xrange));
    op->top_type = problem->boundary.top.type;

    op->h1 = h1;
    op->h2 = h2;

    correct_right_part(op);
    return op;
}

static void delete_operator(Operator *op)
{
    delete_matrix(op->A);
    delete_matrix(op->B);
    delete_matrix(op->Q);
    delete_matrix(op->F);
    delete_matrix(op->PhiL);
    delete_matrix(op->PhiR);
    delete_matrix(op->PhiB);
    delete_matrix(op->PhiT);

    free(op);
}

Matrix *solve(const Problem *problem, const SolvingInfo *config)
{
    Operator *op = create_operator(problem, config);
    Matrix *u = find_solution(op, config->eps);
    delete_operator(op);
    return u;
}