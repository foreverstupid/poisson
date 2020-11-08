#include "process.h"

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



ProcessInfo *get_process_info(
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