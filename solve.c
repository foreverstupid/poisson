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
 * Stores calculated data.
 * Note: if the process has neighbours then its data has shadow edges,
 * so we should prepare the calculated data before using write function.
 */
static void output_data(
    const Matrix *m,
    scalar_t *buf,
    int iteration_idx,
    const ProcessInfo *info)
{
    int i;
    int j;

    int x_start_shift = info->left_neighbour == MPI_PROC_NULL ? 0 : 1;
    int x_end_shift = info->right_neighbour == MPI_PROC_NULL ? 0 : 1;

    int y_start_shift = info->bottom_neighbour == MPI_PROC_NULL ? 0 : 1;
    int y_end_shift = info->top_neighbour == MPI_PROC_NULL ? 0 : 1;

    Matrix tmp = {
        .data = buf,
        .nx = m->nx - x_start_shift - x_end_shift,
        .ny = m->ny - y_start_shift - y_end_shift
    };

    for (j = y_start_shift; j < m->ny - y_end_shift; j++)
    {
        for (i = x_start_shift; i < m->nx - x_end_shift; i++)
        {
            set(&tmp, i - x_start_shift, j - y_start_shift, at(m, i, j));
        }
    }

    info->log.write_matrix(&tmp, iteration_idx);
}



/*
 * Performs iteration using the given initial function. It stores
 * the solution in the same argument. The exit creteria for the
 * iterational method is ||r||_C < eps.
 */
static void find_solution(Matrix **u, const ProcessInfo *info)
{
    Matrix *v = new_matrix((*u)->nx, (*u)->ny);
    Matrix *tmp;
    Matrix *buf = new_matrix((*u)->nx, (*u)->ny);

    int iteration_idx = 0;
    scalar_t curr_eps;
    scalar_t local_eps;

    apply_first_order_boundary(v, info->op);
    exchange_boundaries(*u, info, buf->data);

    do
    {
        if (info->log.iteration_print_frequency > 0 &&
            iteration_idx % info->log.iteration_print_frequency == 0)
        {
            output_data(*u, buf->data, iteration_idx, info);
            if (iteration_idx != 0)
            {
                info->log.log_message("(%d) -> "SFI, iteration_idx, local_eps);
            }
        }

        local_eps = perform_iteration(v, info->op, *u, buf, &(info->mask));

        tmp = *u;
        *u = v;
        v = tmp;

        iteration_idx++;
        exchange_boundaries(*u, info, buf->data);

        MPI_Allreduce(
            &local_eps, &curr_eps, 1, MPI_SCALAR,
            MPI_MAX, info->comm);
    }
    while (curr_eps >= info->eps);

    info->log.log_message(
        "Last iteration: %d, "
        "local difference: "SFI", "
        "global difference: "SFI,
        iteration_idx, local_eps, curr_eps);

    delete_matrix(v);
    delete_matrix(buf);
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
    scalar_t *buf = (scalar_t *)malloc(u->nx * u->ny * sizeof(scalar_t));

    apply_first_order_boundary(u, info->op);

    config->log.log_message("Starting solving process...");
    find_solution(&(u), info);
    config->log.log_message("Solution is found");

    output_data(u, buf, -1, info);

    delete_matrix(u);
    delete_operator(info->op);
    free(info);
    free(buf);
}