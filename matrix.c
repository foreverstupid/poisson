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
    scalar_t xstep = get_range_step(&(ranges->xrange));
    scalar_t ystep = get_range_step(&(ranges->yrange));

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
    scalar_t step = get_range_step(range);
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
    scalar_t step = get_range_step(range);
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
    const MatrixMask *mask,
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

    for (j = mask->y0; j <= mask->y1; j++)
    {
        ry = j == 0 || j == N ? 0.5 * h2 : h2;
        for (i = mask->x0; i <= mask->x1; i++)
        {
            rx = i == 0 || i == M ? 0.5 * h1 : h1;
            res += at(m1, i, j) * at(m2, i, j) * rx * ry;
        }
    }

    return res;
}



void multiply(Matrix *matrix, scalar_t f, const MatrixMask *mask)
{
    int i;
    int j;
    scalar_t tmp;

    for (j = mask->y0; j <= mask->y1; j++)
    {
        for (i = mask->x0; i <= mask->x1; i++)
        {
            tmp = at(matrix, i, j) * f;
            set(matrix, i, j, tmp);
        }
    }
}



void sub(Matrix *m1, const Matrix *m2, const MatrixMask *mask)
{
    int i;
    int j;
    scalar_t tmp;

    for (j = mask->y0; j <= mask->y1; j++)
    {
        for (i = mask->x0; i <= mask->x1; i++)
        {
            tmp = at(m1, i, j) - at(m2, i, j);
            set(m1, i, j, tmp);
        }
    }
}



void linear_combination(
    Matrix *res,
    const Matrix *m1,
    scalar_t t,
    const Matrix *m2,
    const MatrixMask *mask)
{
    int i;
    int j;
    scalar_t tmp;

    for (j = mask->y0; j <= mask->y1; j++)
    {
        for (i = mask->x0; i <= mask->x1; i++)
        {
            tmp = at(m1, i, j) + t * at(m2, i, j);
            set(res, i, j, tmp);
        }
    }
}



scalar_t get_cnorm(const Matrix *m, const MatrixMask *mask)
{
    scalar_t cnorm = 0.0;
    int i;
    int j;

    for (j = mask->y0; j <= mask->y1; j++)
    {
        for (i = mask->x0; i <= mask->x1; i++)
        {
            cnorm = fmax(cnorm, fabs(at(m, i, j)));
        }
    }

    return cnorm;
}



/*
 * Gets the given elements forn the data.
 */
static void get_elements(
    const scalar_t *data,
    int start,
    int step,
    int cnt,
    scalar_t *buf)
{
    int i;
    int j;

    for (i = start, j = 0; j < cnt; i += step, j++)
    {
        buf[j] = data[i];
    }
}



/*
 * Sets the given elements to the data.
 */
static void set_elements(
    scalar_t *data,
    int start,
    int step,
    int cnt,
    const scalar_t *buf)
{
    int i;
    int j;

    for (i = start, j = 0; j < cnt; i += step, j++)
    {
        data[i] = buf[j];
    }
}



void get_row(const Matrix *m, int i, scalar_t *row)
{
    get_elements(m->data, i * m->nx, 1, m->nx, row);
}



void set_row(Matrix *m, int i, scalar_t *row)
{
    set_elements(m->data, i * m->nx, 1, m->nx, row);
}



void get_column(const Matrix *m, int i, scalar_t *column)
{
    get_elements(m->data, i, m->nx, m->ny, column);
}



void set_column(Matrix *m, int i, scalar_t *column)
{
    set_elements(m->data, i, m->nx, m->ny, column);
}