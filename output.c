#include "output.h"

/*
 * The path for matrix writing.
 */
char *path = NULL;

/*
 * The length of the matrix output base path.
 */
int base_path_length;

/*
 * X-coordinate of the current process.
 */
int x_proc_coord;

/*
 * Y-coordinate of the current process.
 */
int y_proc_coord;

/*
 * Whether the module is initialized or not.
 */
int initialized = 0;

/* help variables for not output shadow edges */

int x_start_shift;
int x_end_shift;
int y_start_shift;
int y_end_shift;



/*
 * Gets the length of the given string.
 */
static int get_length(const char *str)
{
    int i;
    for (i = 0; str[i]; i++)
    {
    }

    return i;
}



/*
 * Creates the given directory recoursively if it dowsn't exist.
 */
static InitResult recoursive_mkdir(char *path)
{
    int i;
    struct stat s = { 0 };

    if (stat(path, &s) == 0)
    {
        return success;
    }

    for (i = 0; path[i]; i++)
    {
        if (path[i] == '/')
        {
            path[i] = 0;
            if (stat(path, &s) == -1 && mkdir(path, 0777) != 0)
            {
                return creating_path_error;
            }

            path[i] = '/';
        }
    }

    return success;
}



InitResult init_output(
    const char *out_dir,
    int x,
    int y,
    int x_count,
    int y_count)
{
    if (initialized )
    {
        return already_init;
    }

    InitResult res;
    int len = get_length(out_dir);
    path = (char *)malloc((len + 80) * sizeof(char));

    sprintf(path, "%s/%d-%d/", out_dir, x, y);
    base_path_length = get_length(path);

    res = recoursive_mkdir(path);
    if (res == success)
    {
        initialized = 1;
        x_proc_coord = x;
        y_proc_coord = y;
        x_start_shift = x == 0 ? 0 : 1;
        x_end_shift = x == x_count - 1 ? 0 : 1;
        y_start_shift = y == 0 ? 0 : 1;
        y_end_shift = y == y_count - 1 ? 0 : 1;
    }
    else
    {
        dispose_output();
    }

    return res;
}



void log_info(const char *message)
{
    printf("[%d:%d] %s\n", x_proc_coord, y_proc_coord, message);
}



void write_matrix(const Matrix *m, int iteration_idx)
{
    int i;
    int j;

    if (iteration_idx >= 0)
    {
        sprintf(path + base_path_length, "u_%d.csv", iteration_idx);
    }
    else
    {
        sprintf(path + base_path_length, "solution.csv");
    }

    FILE *out = fopen(path, "w");

    for (j = y_start_shift; j < m->ny - y_end_shift; j++)
    {
        fprintf(out, SF, at(m, x_start_shift, j));
        for (i = x_start_shift + 1; i < m->nx - x_end_shift; i++)
        {
            fprintf(out, " " SF, at(m, i, j));
        }

        putc('\n', out);
    }

    fclose(out);
}



void dispose_output()
{
    if (path != NULL)
    {
        free(path);
        path = NULL;
    }
}