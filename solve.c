#include "solve.h"

/*
 * Gets samples of a two variable function.
 */
static Vector *get_samples(
    Func2 func,
    const Range *x_range,
    const Range *y_range)
{
    int i;
    int j;
    scalar_t x = x_range->start;
    scalar_t y = y_range->start;
    scalar_t x_step = get_step(x_range);
    scalar_t y_step = get_step(y_range);

    Vector *result = new_vector(x_range->count * y_range->count);

    for (i = 0; i < x_range->count; i++)
    {
        y = y_range->start;
        for (j = 0; j < y_range->count; j++)
        {
            result->data[i * y_range->count + j] = func(x, y);
            y += y_step;
        }

        x += x_step;
    }

    return result;
}



/*
 * Performs a single iteration (u -> w).
 */
static void perform_iteration(
    const Vector *u,
    Vector *w,
    const Vector *K,
    const Vector *Q,
    const Vector *F,
    const Vector *Phi,
    const Problem *problem)
{

}



/*
 * Performs iteration process getting the solution in the first
 * argument.
 */
static void iteration_process(
    Vector *u,
    const Vector *K,
    const Vector *Q,
    const Vector *F,
    const Vector *Phi,
    const Problem *problem,
    scalar_t eps)
{
    Vector *buf;
    Vector *w = new_vector(problem->x_range.count * problem->y_range.count);
    scalar_t current_delta;

    do
    {
        perform_iteration(u, w, K, Q, F, Phi, problem);
        current_delta = get_difference_norm(w, u);
        buf = w;
        w = u;
        u = buf;
    }
    while (current_delta > eps);

    u = w;
}



Vector *solve(const Problem *problem, scalar_t eps)
{
    Vector *K = get_samples(problem->k, &(problem->x_range), &(problem->y_range));
    Vector *Q = get_samples(problem->q, &(problem->x_range), &(problem->y_range));
    Vector *F = get_samples(problem->f, &(problem->x_range), &(problem->y_range));
    Vector *Phi = get_samples(problem->phi, &(problem->x_range), &(problem->y_range));

    double h1 = get_step(&(problem->x_range));
    double h2 = get_step(&(problem->y_range));

    Vector *u = new_vector(problem->x_range.count * problem->y_range.count);
    iteration_process(u, K, Q, F, Phi, problem, eps);

    delete_vector(K);
    delete_vector(Q);
    delete_vector(Phi);
    delete_vector(F);

    return u;
}