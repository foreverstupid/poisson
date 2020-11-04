#include "output.h"

void write_as_csv(const Matrix *m, const char *file_name)
{
    int i;
    int j;

    FILE *out = fopen(file_name, "w");
    for (j = 0; j < m->ny; j++)
    {
        fprintf(out, S_FORMAT, at(m, 0, j));
        for (i = 1; i < m->nx; i++)
        {
            fprintf(out, ", " S_FORMAT, at(m, i, j));
        }

        putc('\n', out);
    }

    fclose(out);
}