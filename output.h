#ifndef HPC_POISSON_OTPUT_MODULE_H
#define HPC_POISSON_OTPUT_MODULE_H

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdarg.h>
#include <errno.h>
#include <string.h>

#include "definitions.h"
#include "matrix.h"

typedef enum InitResult{
    success = 0,
    already_init = 1,
    creating_path_error = 2
} InitResult;

/*
 * Initializes the output information from the path to the output
 * directory and XY-coordinates of the process in the virtual
 * topology. For using the module you should successfully call this
 * function at least once. After that all others invocations are ignored.
 */
InitResult init_output(const char *out_dir, int x, int y);

/*
 * Logs the information message.
 */
void log_info(const char *format, ...);

/*
 * Stores the matrix. Negative iteration index for the final solution.
 */
void write_matrix(const Matrix *m, int iteration_idx);

/*
 * Disposes the resources, that was taken by the module.
 */
void dispose_output();

#endif
