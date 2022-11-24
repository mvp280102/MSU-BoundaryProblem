#include "solver.h"

double right_expr(double x, double z1, double z2)
{
	return -1 * px_func(x) * z2 - qx_func(x) * z1 + fxy_func(x, z1);
}

double calculate_const(double x, double y, double z)
{
	return y - z * x;
}

void z1_step(uint idx, double *z1, double x, double z2, double c)
{
	z1[idx] = z2 * x + c;
}

void z2_step(uint idx, double *z2, double *x, double *z1, double h)
{
	// Predictor step:
	z2[idx] = z2[idx - 1] + h * (1901.0 / 720 * right_expr(x[idx - 1], z2[idx - 1], z1[idx - 1]) -
								 1387.0 / 360 * right_expr(x[idx - 2], z2[idx - 2], z1[idx - 2]) +
								 109.0  / 30  * right_expr(x[idx - 3], z2[idx - 3], z1[idx - 3]) -
								 637.0  / 360 * right_expr(x[idx - 4], z2[idx - 4], z1[idx - 4]) +
								 251.0  / 720 * right_expr(x[idx - 5], z2[idx - 5], z1[idx - 5]));

	// Corrector step:
	z2[idx] = z2[idx - 1] + h / 720 * (251 * right_expr(x[idx], z2[idx], z1[idx]) +
									   646 * right_expr(x[idx - 1], z2[idx - 1], z1[idx - 1]) -
									   264 * right_expr(x[idx - 2], z2[idx - 2], z1[idx - 2]) +
									   106 * right_expr(x[idx - 3], z2[idx - 3], z1[idx - 3]) -
									   19  * right_expr(x[idx - 4], z2[idx - 4], z1[idx - 4]));
}