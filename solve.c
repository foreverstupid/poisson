#include "solve.h"
#include "output.h"

/*
 * Performs a single iteration of the method of the smallest residuals
 * (u -> v).
 */
static void perform_iteration(
    Matrix *v,
    const Operator *op,
    const Matrix *u,
    Matrix *buf)
{
    scalar_t tmp;

    scalar_t tau;                   /* iterational parameter */
    apply(v, op, u);
    sub(v, op->F);                  /* residual */
    apply(buf, op, v);              /* operator applied to residual */

    tmp = dot_product(buf, v, op->h1, op->h2);
    if (fabs(tmp) > EPS)
    {
        tau = tmp / get_squared_norm(buf, op->h1, op->h2);
    }
    else
    {
        tau = 0.0;
    }

    linear_combination(v, u, -tau, v);
}

/*
 * Performs iteration process until reaching the given error
 * between iterations uding the given initial function <u>. It stores
 * the solution in the same argument.
 */
void find_solution(const Operator *op, Matrix **u, const SolvingInfo *config)
{
    Matrix *v = new_matrix((*u)->nx, (*u)->ny);
    Matrix *tmp;
    Matrix *buf = new_matrix((*u)->nx, (*u)->ny);
    scalar_t curr_eps;

    int iteration_idx = 0;
    char name[80];

    do
    {
        perform_iteration(v, op, *u, buf);
        curr_eps = get_difference_cnorm(*u, v);

        tmp = *u;
        *u = v;
        v = tmp;

        if (config->iteration_print_frequency > 0 &&
            iteration_idx % config->iteration_print_frequency == 0)
        {
            sprintf(name, "%s/i%d.csv", config->output_dir, iteration_idx);
            write_as_csv(v, name);
        }

        iteration_idx++;
    }
    while (curr_eps > config->eps);

    delete_matrix(v);
    delete_matrix(buf);
}



/*
 * Initializes the first iteration of the solving process.
 */
void init_first_iteration(
    Matrix *u,
    const Operator *op,
    const Problem *p)
{
    int i;
    scalar_t t;

    if (op->left_type == first)
    {
        t = p->area.y1;
        for (i = 0; i < u->ny; i++)
        {
            set(u, 0, i, p->boundary.left.phi(t));
            t += op->h2;
        }
    }

    if (op->right_type == first)
    {
        t = p->area.y1;
        for (i = 0; i < u->ny; i++)
        {
            set(u, u->nx - 1, i, p->boundary.right.phi(t));
            t += op->h2;
        }
    }

    if (op->bottom_type == first)
    {
        t = p->area.x1;
        for (i = 0; i < u->nx; i++)
        {
            set(u, i, 0, p->boundary.bottom.phi(t));
            t += op->h1;
        }
    }

    if (op->top_type == first)
    {
        t = p->area.x1;
        for (i = 0; i < u->nx; i++)
        {
            set(u, i, u->ny - 1, p->boundary.top.phi(t));
            t += op->h1;
        }
    }
}



Matrix *solve(const Problem *problem, const SolvingInfo *config)
{
    Operator *op = new_operator(problem, config);
    Matrix *u = new_matrix(op->F->nx, op->F->ny);
    init_first_iteration(u, op, problem);

    find_solution(op, &u, config);
    delete_operator(op);
    return u;
}