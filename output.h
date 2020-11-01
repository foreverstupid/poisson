#ifndef HPC_POISSON_OTPUT_MODULE_H
#define HPC_POISSON_OTPUT_MODULE_H

#include "definitions.h"
#include "matrix.h"
#include "stdio.h"

/*
 * Stores the matrix into the CSV format.
 */
void write_as_csv(const Matrix *m, const char *file_name);

#endif