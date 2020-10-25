#include "matrix.h"

Matrix *new_matrix(int nrows, int ncols)
{
    Matrix *result = (Matrix *)malloc(sizeof(Matrix));

    result->nrows = nrows;
    result->ncols = ncols;
    result->data = (scalar_t *)malloc(nrows * ncols * sizeof(scalar_t));

    return result;
}



Matrix *get_matrix_from_func2(Func2 func, const Range2 *ranges)
{
    int i;
    int j;
    scalar_t x;
    scalar_t y = ranges->yrange.start;
    scalar_t xstep = get_step(&(ranges->xrange));
    scalar_t ystep = get_step(&(ranges->yrange));

    Matrix *result = new_matrix(ranges->yrange.count, ranges->xrange.count);
    scalar_t* data = result->data;

    for (i = 0; i < ranges->yrange.count; i++)
    {
        x = ranges->xrange.start;
        for (j = 0; j < ranges->xrange.count; j++)
        {
            data[j] = func(x, y);
            x += xstep;
        }

        y += ystep;
        data += ranges->xrange.count;
    }

    return result;
}



Matrix *copy(const Matrix *matrix)
{
    Matrix *result = new_matrix(matrix->nrows, matrix->ncols);
    int i;

    for (i = 0; i < matrix->nrows * matrix->ncols; i++)
    {
        result->data[i] = matrix->data[i];
    }

    return result;
}



void delete_matrix(Matrix *matrix)
{
    free(matrix->data);
    free(matrix);
}



scalar_t dot_product(const Matrix *m1, const Matrix *m2)
{
    int i;
    double res = m1->data[0] * m2->data[0];

    for (i = 1; i < m1->nrows * m1->ncols; i++)
    {
        res += m1->data[i] * m2->data[i];
    }

    return res;
}



void multiply(Matrix *matrix, scalar_t f)
{
    int i;
    for (i = 0; i < matrix->nrows * matrix->ncols; i++)
    {
        matrix->data[i] *= f;
    }
}



void add(Matrix *m1, const Matrix *m2)
{
    int i;
    for (i = 0; i < m1->nrows * m2->ncols; i++)
    {
        m1->data[i] += m2->data[i];
    }
}



void sub(Matrix *m1, const Matrix *m2)
{
    int i;
    for (i = 0; i < m1->nrows * m1->ncols; i++)
    {
        m1->data[i] -= m2->data[i];
    }
}



scalar_t get_difference_cnorm(const Matrix *m1, const Matrix *m2)
{
    scalar_t cnorm = fabs(m1->data[0] - m2->data[0]);
    int i;

    for (i = 1; i < m1->nrows * m1->ncols; i++)
    {
        cnorm = fmax(cnorm, fabs(m1->data[i] - m2->data[i]));
    }

    return cnorm;
}