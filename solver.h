#pragma once

#include "px_func.h"
#include "qx_func.h"
#include "fxy_func.h"

typedef unsigned int uint;

double total_func(double x, double y, double z);
void solving_step(uint index, double *x, double *y, double *z, double h);