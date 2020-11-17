#include "solve.h"

/*
 * Performs a single iteration of the method of the smallest residuals
 * (v = u - \tau * r). Returns C-norm of the residual on the previous
 * iteration. Thus, we can do at most one extra iteration.
 */
static scalar_t perform_iteration(
    Matrix *v,
    const Operator *op,
    const Matrix *u,
    Matrix *buf,
    const MatrixMask *mask)
{
    scalar_t r_norm;
    scalar_t tau;                   /* iterational parameter */

    apply(v, op, u);
    sub(v, op->F, mask);            /* residual */
    r_norm = get_cnorm(v, mask);

    apply(buf, op, v);              /* operator applied to residual */

    tau =
        dot_product(buf, v, mask, op->h1, op->h2) /
        get_squared_norm(buf, mask, op->h1, op->h2);

    linear_combination(v, u, -tau, v, mask);
    return r_norm;
}



/*
 * Exchanges boundaries info among processors.
 */
static void exchange_boundaries(
    Matrix *u,
    const ProcessInfo *info,
    scalar_t *buf)
{
    MPI_Status status;
    int M = u->nx - 1;
    int N = u->ny - 1;

    if (info->left_neighbour != MPI_PROC_NULL)
    {
        get_column(u, 1, buf);
        MPI_Sendrecv_replace(
            buf, u->ny, MPI_SCALAR,
            info->left_neighbour, 0, info->left_neighbour, 0,
            info->comm, &status);
        
        set_column(u, 0, buf);
    }

    if (info->right_neighbour != MPI_PROC_NULL)
    {
        get_column(u, M - 1, buf);
        MPI_Sendrecv_replace(
            buf, u->ny, MPI_SCALAR,
            info->right_neighbour, 0, info->right_neighbour, 0,
            info->comm, &status);

        set_column(u, M, buf);
    }

    if (info->bottom_neighbour != MPI_PROC_NULL)
    {
        get_row(u, 1, buf);
        MPI_Sendrecv_replace(
            buf, u->nx, MPI_SCALAR,
            info->bottom_neighbour, 0, info->bottom_neighbour, 0,
            info->comm, &status);

        set_row(u, 0, buf);
    }

    if (info->top_neighbour != MPI_PROC_NULL)
    {
        get_row(u, N - 1, buf);
        MPI_Sendrecv_replace(
            buf, u->nx, MPI_SCALAR,
            info->top_neighbour, 0, info->top_neighbour, 0,
            info->comm, &status);

        set_row(u, N, buf);
    }
}



/*
 * Makes single step of the solving process.
 */
static scalar_t make_step(
    Matrix **v,
    Matrix **u,
    Matrix *buf,
    const ProcessInfo *info,
    scalar_t *local_eps)
{
    scalar_t global_eps;
    Matrix *tmp;

    *local_eps = perform_iteration(*v, info->op, *u, buf, &(info->mask));

    tmp = *u;
    *u = *v;
    *v = tmp;

    exchange_boundaries(*u, info, buf->data);

    MPI_Allreduce(
        local_eps, &global_eps, 1, MPI_SCALAR,
        MPI_MAX, info->comm);

    return global_eps;
}



/*
 * Performs iteration using the given initial function. It stores
 * the solution in the same argument. The exit creteria for the
 * iterational method is ||r||_C < eps or the maximum iteration count
 * is reached. Returns 1 if the solution is actually found (i.e. if
 * residual norm is less than eps), and 0 if the max iteration count
 * is reached.
 */
static int find_solution(Matrix **u, const ProcessInfo *info)
{
    int iteration_idx = -1;
    int global_iteration_idx;
    scalar_t global_eps = 0;
    scalar_t local_eps = 0;

    Matrix *v = new_matrix((*u)->nx, (*u)->ny);
    Matrix *buf = new_matrix((*u)->nx, (*u)->ny);

    apply_first_order_boundary(v, info->op);
    exchange_boundaries(*u, info, buf->data);

    do
    {
        iteration_idx++;
        if (info->log.iteration_print_frequency > 0 &&
            iteration_idx > 0 &&
            iteration_idx % info->log.iteration_print_frequency == 0)
        {
            global_iteration_idx = iteration_idx + info->num.init_iteration;
            info->log.write_matrix(*u, &(info->mask), global_iteration_idx);
            info->log.log_message("(%d) -> "SFI, iteration_idx, local_eps);
        }

        if (iteration_idx < info->num.max_iterations)
        {
            global_eps = make_step(&v, u, buf, info, &local_eps);
        }
    }
    while (global_eps >= info->num.eps &&
           iteration_idx < info->num.max_iterations);

    info->log.log_message(
        "Last iteration: %d (globally %d), "
        "local residual norm: "SFI", "
        "global residual norm: "SFI,
        iteration_idx, iteration_idx + info->num.init_iteration,
        local_eps,
        global_eps);

    delete_matrix(v);
    delete_matrix(buf);

    return global_eps < info->num.eps;
}



void solve(const Problem *problem, const SolvingConfig *config)
{
    if (config->mpi.grid_comm == MPI_COMM_NULL)
    {
        return;
    }

    config->log.log_message("Getting process info...");
    ProcessInfo *info = get_process_info(problem, config);
    Matrix *u = new_matrix(info->op->F->nx, info->op->F->ny);

    config->log.read_matrix(u, &(info->mask), config->num.init_iteration);
    apply_first_order_boundary(u, info->op);

    config->log.log_message("Starting solving process...");
    if (find_solution(&(u), info))
    {
        config->log.log_message("Solution is found");
        info->log.write_matrix(u, &(info->mask), -1);
    }
    else
    {
        config->log.log_message("The maximum iteration count is reached");
    }

    delete_matrix(u);
    delete_operator(info->op);
    free(info);
}