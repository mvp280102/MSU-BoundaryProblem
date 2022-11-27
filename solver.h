#pragma once

#include "px_func.h"
#include "qx_func.h"
#include "fxy_func.h"

#include <cstdio>
#include <cmath>

#define ACC_ORDER 5

typedef unsigned int uint;

struct BoundaryData
{
	uint intervals;

	double arg_a;
	double arg_b;
	double func_a;
	double func_b;
};

double vector_distance(uint len, double *x, double *y1, double *y2);

double numerical_derivative(BoundaryData data, double arg, double *res, double (*func)(BoundaryData, double, double*));
double newton_step(BoundaryData data, double arg, double *res, double (*func)(BoundaryData, double, double*));

double right_expr_z1(double x, double z1, double z2);
double right_expr_z2(double x, double z2, double z1);

void rk_step(uint idx, double *arg, double *cur, double *other, double step, double (*expr)(double, double, double));
void adams_step(uint idx, double *arg, double *cur, double *other, double step, double (*expr)(double, double, double));

double boundary_solve(BoundaryData data, double der_a, double *result);