#pragma once

#include "px_func.h"
#include "qx_func.h"
#include "fxy_func.h"

typedef unsigned int uint;

double right_expr(double x, double z1, double z2);
double calculate_const(double x, double y, double z);

void z1_step(uint idx, double *z1, double x, double z2, double c);
void z2_step(uint idx, double *z2, double *x, double *z1, double h);