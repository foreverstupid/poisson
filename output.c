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



InitResult init_output(const char *out_dir, int x, int y)
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

    for (j = 0; j < m->ny; j++)
    {
        fprintf(out, SF, at(m, 0, j));
        for (i = 1; i < m->nx; i++)
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
