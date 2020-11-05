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
} ProcessInfo;



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
 * Exchanges boundaries info among processors.
 */
static void exchange_boundaries(Matrix *u, const ProcessInfo *info)
{
    MPI_Status status;
    scalar_t *buf;
    int M = u->nx - 1;
    int N = u->ny - 1;

    buf = get_column(u, 0);
    MPI_Sendrecv_replace(
        buf, u->ny, MPI_SCALAR_TYPE,
        info->left_neighbour, 0, info->left_neighbour, 0,
        info->comm, &status);

    set_column(info->op->F, 0, buf);
    free(buf);

    buf = get_column(u, M);
    MPI_Sendrecv_replace(
        buf, u->ny, MPI_SCALAR_TYPE,
        info->right_neighbour, 0, info->right_neighbour, 0,
        info->comm, &status);

    set_column(info->op->F, M, buf);
    free(buf);

    buf = get_row(u, 0);
    MPI_Sendrecv_replace(
        buf, u->nx, MPI_SCALAR_TYPE,
        info->bottom_neighbour, 0, info->bottom_neighbour, 0,
        info->comm, &status);

    set_row(info->op->F, 0, buf);
    free(buf);

    buf = get_row(u, N);
    MPI_Sendrecv_replace(
        buf, u->nx, MPI_SCALAR_TYPE,
        info->top_neighbour, 0, info->top_neighbour, 0,
        info->comm, &status);

    set_row(info->op->F, N, buf);
    free(buf);
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
    int max_iter = 1000;

    do
    {
        perform_iteration(v, info->op, *u, buf);

        tmp = *u;
        *u = v;
        v = tmp;

        if (info->log.iteration_print_frequency > 0 &&
            iteration_idx % info->log.iteration_print_frequency == 0)
        {
            info->log.write_matrix(*u, iteration_idx);
        }

        iteration_idx++;
        exchange_boundaries(*u, info);
    }
    while (iteration_idx < max_iter);

    delete_matrix(v);
    delete_matrix(buf);
}



/*
 * Sets MPI information for the process.
 */
static void set_mpi_info(ProcessInfo *info, const SolvingConfig *config)
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

    int x_start = x_cnt * config->mpi.x_proc_idx;
    int y_start = y_cnt * config->mpi.y_proc_idx;

    /* last row/column can have less nodes due to unable to divide
       equally */
    int x_end =
        config->mpi.x_proc_idx == config->mpi.x_proc_count - 1
        ? config->num.x_grid_count - 1
        : x_cnt * (config->mpi.x_proc_idx + 1) - 1;

    int y_end =
        config->mpi.y_proc_idx == config->mpi.y_proc_count - 1
        ? config->num.y_grid_count - 1
        : y_cnt * (config->mpi.y_proc_idx + 1) - 1;

    scalar_t hx =
        (global->area.x2 - global->area.x1) /
        (config->num.x_grid_count - 1);

    scalar_t hy =
        (global->area.y2 - global->area.y1) /
        (config->num.y_grid_count - 1);

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

    set_mpi_info(proc, config);
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
    Matrix *part = copy_matrix(info->op->F);

    config->log.log_message("Starting solving process...");
    find_solution(&(part), info);
    config->log.log_message("Solution is found");

    config->log.write_matrix(part, -1);
    delete_matrix(part);
    delete_operator(info->op);
    free(info);
}