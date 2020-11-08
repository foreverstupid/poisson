#include "solve.h"

/*
 * Configuration of the process, solves the subproblem.
 */
typedef struct ProcessInfo
{
    /*
     * The cart communicator.
     */
    MPI_Comm comm;

    /*
     * The rank of the current process in the cart communicator.
     */
    int rank;

    /*
     * The rank of the left process neighbour in the cart communicator.
     * MPI_PROC_NULL if the current process has no left neighbour.
     */
    int left_neighbour;

    /*
     * The rank of the right process neighbour in the cart communicator.
     * MPI_PROC_NULL if the current process has no right neighbour.
     */
    int right_neighbour;

    /*
     * The rank of the bottom process neighbour in the cart communicator.
     * MPI_PROC_NULL if the current process has no bottom neighbour.
     */
    int bottom_neighbour;

    /*
     * The rank of the top process neighbour in the cart communicator.
     * MPI_PROC_NULL if the current process has no top neighbour.
     */
    int top_neighbour;

    /*
     * Operator of the subproblem.
     */
    Operator *op;

    /*
     * Log configuration. It is defined only for the main process.
     */
    LogConfig log;

    /*
     * The mask for matrix operation that prevents changing values of
     * the boundary condition of the first order.
     */
    MatrixMask mask;

    /*
     * Difference between iterations that indicates that the process of
     * solving is over.
     */
    scalar_t eps;
} ProcessInfo;



/*
 * Performs a single iteration of the method of the smallest residuals
 * (u -> v).
 */
static void perform_iteration(
    Matrix *v,
    const Operator *op,
    const Matrix *u,
    Matrix *buf,
    const MatrixMask *mask)
{
    scalar_t tmp;

    scalar_t tau;                   /* iterational parameter */
    apply(v, op, u);
    sub(v, op->F, mask);            /* residual */
    apply(buf, op, v);              /* operator applied to residual */

    tmp = dot_product(buf, v, mask, op->h1, op->h2);
    if (fabs(tmp) > EPS)
    {
        tau = tmp / get_squared_norm(buf, mask, op->h1, op->h2);
    }
    else
    {
        tau = 0.0;
    }

    linear_combination(v, u, -tau, v, mask);
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
 * the solution in the same argument.
 */
static void find_solution(Matrix **u, const ProcessInfo *info)
{
    Matrix *v = new_matrix((*u)->nx, (*u)->ny);
    Matrix *tmp;
    Matrix *buf = new_matrix((*u)->nx, (*u)->ny);

    int iteration_idx = 0;
    scalar_t curr_eps;
    scalar_t local_eps;
    char msg[256];

    do
    {
        if (info->log.iteration_print_frequency > 0 &&
            iteration_idx % info->log.iteration_print_frequency == 0)
        {
            output_data(*u, buf->data, iteration_idx, info);
        }

        perform_iteration(v, info->op, *u, buf, &(info->mask));

        tmp = *u;
        *u = v;
        v = tmp;

        local_eps = get_difference_cnorm(*u, v, &(info->mask));
        iteration_idx++;
        exchange_boundaries(*u, info, buf->data);

        MPI_Allreduce(
            &local_eps, &curr_eps, 1, MPI_SCALAR,
            MPI_MAX, info->comm);
    }
    while (curr_eps >= info->eps);

    sprintf(
        msg,
        "Last iteration: %d, "
        "local difference: "SFI", "
        "global difference: "SFI,
        iteration_idx, local_eps, curr_eps);

    info->log.log_message(msg);

    delete_matrix(v);
    delete_matrix(buf);
}



/*
 * Sets MPI information for the process.
 */
static void set_mpi_config(ProcessInfo *info, const SolvingConfig *config)
{
    info->comm = config->mpi.grid_comm;

    MPI_Cart_shift(
        info->comm,
        0,
        1,
        &(info->rank),
        &(info->right_neighbour));

    MPI_Cart_shift(
        info->comm,
        0,
        -1,
        &(info->rank),
        &(info->left_neighbour));

    MPI_Cart_shift(
        info->comm,
        1,
        1,
        &(info->rank),
        &(info->top_neighbour));

    MPI_Cart_shift(
        info->comm,
        1,
        -1,
        &(info->rank),
        &(info->bottom_neighbour));
}



/*
 * Help function for inner process boundary confdition.
 */
static scalar_t init_boundary(scalar_t t)
{
    return 0.0;
}



/*
 * Gets the part of the area for the process.
 */
static void get_part_problem(
    Problem *part,
    NumericalConfig *num,
    const Problem *global,
    const SolvingConfig *config)
{
    BoundaryCondition initial = {
        .phi = init_boundary,
        .type = first
    };

    int x_cnt = ceil(config->num.x_grid_count / config->mpi.x_proc_count);
    int y_cnt = ceil(config->num.y_grid_count / config->mpi.y_proc_count);

    int x_start =
        config->mpi.x_proc_idx == 0
        ? 0
        : x_cnt * config->mpi.x_proc_idx - 1;

    int y_start =
        config->mpi.y_proc_idx == 0
        ? 0
        : y_cnt * config->mpi.y_proc_idx - 1;

    int x_end =
        config->mpi.x_proc_idx == config->mpi.x_proc_count - 1
        ? config->num.x_grid_count - 1
        : x_cnt * (config->mpi.x_proc_idx + 1);

    int y_end =
        config->mpi.y_proc_idx == config->mpi.y_proc_count - 1
        ? config->num.y_grid_count - 1
        : y_cnt * (config->mpi.y_proc_idx + 1);

    scalar_t hx = get_step(
        global->area.x1,
        global->area.x2,
        config->num.x_grid_count);

    scalar_t hy = get_step(
        global->area.y1,
        global->area.y2,
        config->num.y_grid_count);

    part->area.x1 = hx * x_start;
    part->area.x2 = hx * x_end;
    part->area.y1 = hy * y_start;
    part->area.y2 = hy * y_end;

    part->boundary.left =
        config->mpi.x_proc_idx == 0
        ? global->boundary.left
        : initial;

    part->boundary.right =
        config->mpi.x_proc_idx == config->mpi.x_proc_count - 1
        ? global->boundary.right
        : initial;

    part->boundary.bottom =
        config->mpi.y_proc_idx == 0
        ? global->boundary.bottom
        : initial;

    part->boundary.top =
        config->mpi.y_proc_idx == config->mpi.y_proc_count - 1
        ? global->boundary.top
        : initial;

    part->funcs = global->funcs;

    num->eps = config->num.eps;
    num->x_grid_count = x_end - x_start + 1;
    num->y_grid_count = y_end - y_start + 1;

    char msg[128];
    sprintf(
        msg, "Area: ["SF", "SF"] x ["SF", "SF"]",
        part->area.x1, part->area.x2, part->area.y1, part->area.y2);

    config->log.log_message(msg);
}



/*
 * Sets the dot product calculating configuration.
 */
static void set_dot_config(ProcessInfo *proc)
{
    proc->mask.x0 = proc->op->left == first ? 1 : 0;
    proc->mask.y0 = proc->op->bottom == first ? 1 : 0;

    proc->mask.x1 =
        proc->op->right == first
        ? proc->op->F->nx - 2
        : proc->op->F->nx - 1;

    proc->mask.y1 =
        proc->op->top == first
        ? proc->op->F->ny - 2
        : proc->op->F->ny - 1;
}



/*
 * Initializes all needed for solution process information.
 */
static ProcessInfo *get_processor_info(
    const Problem *global_problem,
    const SolvingConfig *config)
{
    Problem part_problem;
    NumericalConfig part_num;
    ProcessInfo *proc = (ProcessInfo *)malloc(sizeof(ProcessInfo));

    get_part_problem(&part_problem, &part_num, global_problem, config);
    proc->op = new_operator(&part_problem, &part_num);
    proc->log = config->log;
    proc->eps = config->num.eps;

    set_dot_config(proc);
    set_mpi_config(proc, config);

    return proc;
}



void solve(const Problem *problem, const SolvingConfig *config)
{
    if (config->mpi.grid_comm == MPI_COMM_NULL)
    {
        return;
    }

    config->log.log_message("Getting process info...");
    ProcessInfo *info = get_processor_info(problem, config);
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