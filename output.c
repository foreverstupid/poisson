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



InitResult init_output(const char *out_dir, int x, int y)
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