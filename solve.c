#include "solve.h"

/*
 * Performs a single iteration (u -> w).
 */
static void perform_iteration(
    const Matrix *u,
    Matrix *w,
    const Matrix *K,
    const Matrix *Q,
    const Matrix *F,
    const Matrix *Phi,
    const Problem *problem)
{

}



/*
 * Performs iteration process getting the solution in the first
 * argument.
 */
static void iteration_process(
    Matrix *u,
    const Matrix *K,
    const Matrix *Q,
    const Matrix *F,
    const Matrix *Phi,
    const Problem *problem,
    scalar_t eps)
{
    Matrix *buf;
    Matrix *w = new_matrix(
        problem->ranges.yrange.count,
        problem->ranges.xrange.count);

    scalar_t current_delta;

    do
    {
        perform_iteration(u, w, K, Q, F, Phi, problem);
        current_delta = get_difference_cnorm(w, u);
        buf = w;
        w = u;
        u = buf;
    }
    while (current_delta > eps);

    u = w;
}



Matrix *solve(const Problem *problem, scalar_t eps)
{
    Matrix *K = get_matrix_from_func2(problem->funcs.k, &(problem->ranges));
    Matrix *Q = get_matrix_from_func2(problem->funcs.q, &(problem->ranges));
    Matrix *F = get_matrix_from_func2(problem->funcs.f, &(problem->ranges));
    Matrix *Phi = get_matrix_from_func2(problem->funcs.phi, &(problem->ranges));

    double h1 = get_step(&(problem->ranges.xrange));
    double h2 = get_step(&(problem->ranges.yrange));

    Matrix *u = new_matrix(
        problem->ranges.yrange.count,
        problem->ranges.xrange.count);

    iteration_process(u, K, Q, F, Phi, problem, eps);

    delete_matrix(K);
    delete_matrix(Q);
    delete_matrix(F);
    delete_matrix(Phi);

    return u;
}