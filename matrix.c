#include "matrix.h"

Matrix *new_matrix(int nx, int ny)
{
    Matrix *result = (Matrix *)malloc(sizeof(Matrix));

    result->nx = nx;
    result->ny = ny;
    result->data = (scalar_t *)malloc(nx * ny * sizeof(scalar_t));

    return result;
}



Matrix *get_matrix_from_func2(Func2 func, const Range2 *ranges)
{
    int i = 0;
    int j;
    scalar_t x;
    scalar_t y = ranges->yrange.start;
    scalar_t xstep = get_step(&(ranges->xrange));
    scalar_t ystep = get_step(&(ranges->yrange));

    Matrix *result = new_matrix(ranges->xrange.count, ranges->yrange.count);
    scalar_t* data = result->data;

    for (j = 0; j < ranges->yrange.count; j++)
    {
        x = ranges->xrange.start;
        for (i = 0; i < ranges->xrange.count; i++)
        {
            data[i] = func(x, y);
            x += xstep;
        }

        y += ystep;
        data += ranges->xrange.count;
    }

    return result;
}



Matrix *get_column_from_func(Func func, const Range *range)
{
    Matrix *res = new_matrix(1, range->count);
    scalar_t y = range->start;
    scalar_t step = get_step(range);
    int j;

    for (j = 0; j < res->ny; j++)
    {
        res->data[j] = func(y);
        y += step;
    }

    return res;
}



Matrix *get_row_from_func(Func func, const Range *range)
{
    Matrix *res = new_matrix(range->count, 1);
    scalar_t x = range->start;
    scalar_t step = get_step(range);
    int i;

    for (i = 0; i < res->nx; i++)
    {
        res->data[i] = func(x);
        x += step;
    }

    return res;
}



Matrix *copy_matrix(const Matrix *matrix)
{
    Matrix *result = new_matrix(matrix->nx, matrix->ny);
    int i;

    for (i = 0; i < matrix->nx * matrix->ny; i++)
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



scalar_t dot_product(
    const Matrix *m1,
    const Matrix *m2,
    scalar_t h1,
    scalar_t h2)
{
    int i;
    int j;
    scalar_t res = 0.0;
    scalar_t rx;
    scalar_t ry;
    int M = m1->nx - 1;
    int N = m1->ny - 1;

    for (j = 0; j < m1->ny; j++)
    {
        ry = j == 0 || j == N ? 0.5 * h2 : h2;
        for (i = 0; i < m1->nx; i++)
        {
            rx = i == 0 || i == M ? 0.5 * h1 : h1;
            res += at(m1, i, j) * at(m2, i, j) * rx * ry;
        }
    }

    return res;
}



void multiply(Matrix *matrix, scalar_t f)
{
    int i;
    for (i = 0; i < matrix->nx * matrix->ny; i++)
    {
        matrix->data[i] *= f;
    }
}



void linear_combination(
    Matrix *res,
    const Matrix *m1,
    scalar_t t,
    const Matrix *m2)
{
    int i;
    for (i = 0; i < m1->nx * m1->ny; i++)
    {
        res->data[i] = m1->data[i] + t * m2->data[i];
    }
}



void sub(Matrix *m1, const Matrix *m2)
{
    int i;
    for (i = 0; i < m1->nx * m1->ny; i++)
    {
        m1->data[i] -= m2->data[i];
    }
}



scalar_t get_difference_cnorm(const Matrix *m1, const Matrix *m2)
{
    scalar_t cnorm = fabs(m1->data[0] - m2->data[0]);
    int i;

    for (i = 1; i < m1->nx * m1->ny; i++)
    {
        cnorm = fmax(cnorm, fabs(m1->data[i] - m2->data[i]));
    }

    return cnorm;
}