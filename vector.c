#include "vector.h"

Vector *new_vector(int length)
{
    Vector *result = (Vector *)malloc(sizeof(Vector));

    result->length = length;
    result->data = (scalar_t *)malloc(length * sizeof(scalar_t));

    return result;
}



Vector *get_vector_from_func(
    Func func,
    scalar_t start,
    scalar_t end,
    int length)
{
    Vector *result = new_vector(length);
    int i;
    scalar_t x = start;
    scalar_t step = (end - start) / (length - 1);

    for (i = 0; i < length; i++)
    {
        result->data[i] = func(x);
        x += step;
    }

    return result;
}



Vector *copy(Vector *vector)
{
    Vector *result = new_vector(vector->length);
    result->length = vector->length;
    int i;

    for (i = 0; i < vector->length; i++)
    {
        result->data[i] = vector->data[i];
    }

    return result;
}



void delete_vector(Vector *vector)
{
    free(vector->data);
    free(vector);
}



scalar_t dot_product(Vector *v1, Vector *v2)
{
    int i;
    double res = v1->data[0] * v2->data[0];

    for (i = 1; i < v1->length; i++)
    {
        res += v1->data[i] * v2->data[i];
    }

    return res;
}



scalar_t get_norm(Vector *vector)
{
    return sqrt(dot_product(vector, vector));
}



void multiply(Vector *vector, scalar_t f)
{
    int i;
    for (i = 0; i < vector->length; i++)
    {
        vector->data[i] *= f;
    }
}



void add(Vector *v1, Vector *v2)
{
    int i;
    for (i = 0; i < v1->length; i++)
    {
        v1->data[i] += v2->data[i];
    }
}



void sub(Vector *v1, Vector *v2)
{
    int i;
    for (i = 0; i < v1->length; i++)
    {
        v1->data[i] -= v2->data[i];
    }
}