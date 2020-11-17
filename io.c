#include "io.h"

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



IOInitResult init_io(const char *out_dir, int x, int y)
{
    if (initialized )
    {
        return already_init;
    }

    int len = get_length(out_dir);
    struct stat s = { 0 };
    path = (char *)malloc((len + 80) * sizeof(char));

    sprintf(path, "%s", out_dir);
    if (path[len - 1] == '/')
    {
        path[len - 1] = 0;
        len--;
    }

    sprintf(path + len, "/%d-%d/", x, y);
    base_path_length = get_length(path);

    if (stat(path, &s) == -1 && mkdir(path, ACCESSPERMS) != 0)
    {
        fprintf(stderr, "ERROR (%s): %s\n", path, strerror(errno));
        dispose_output();
        return creating_path_error;
    }

    initialized = 1;
    x_proc_coord = x;
    y_proc_coord = y;

    return success;
}



void log_info(const char *format, ...)
{
    va_list args;
    va_start(args, format);
    printf("[%d:%d] ", x_proc_coord, y_proc_coord);
    vprintf(format, args);
    putchar('\n');
}



void write_matrix(
    const Matrix *m,
    const MatrixMask *mask,
    int iteration_idx)
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

    for (j = mask->y0; j <= mask->y1; j++)
    {
        fprintf(out, SF, at(m, mask->x0, j));
        for (i = mask->x0 + 1; i <= mask->x1; i++)
        {
            fprintf(out, " " SF, at(m, i, j));
        }

        putc('\n', out);
    }

    fclose(out);
}



void read_matrix(
    Matrix *m,
    const MatrixMask *mask,
    int iteration_idx)
{
    int i;
    int j;

    if (iteration_idx <= 0)
    {
        return;
    }

    sprintf(path + base_path_length, "u_%d.csv", iteration_idx);
    FILE *out = fopen(path, "r");
    scalar_t tmp;

    for (j = mask->y0; j <= mask->y1; j++)
    {
        for (i = mask->x0; i <= mask->x1; i++)
        {
            if (fscanf(out, SF, &tmp) != 1)
            {
                log_info(
                    "ERROR: cannot read element (%d, %d) of iteration %d",
                    i, j, iteration_idx);
            }

            set(m, i, j, tmp);
        }
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
