#include "solver.h"

double right_expr_z1(double x, double z1, double z2)
{
	return z2;
}

double right_expr_z2(double x, double z2, double z1)
{
	return -1 * px_func(x) * z2 - qx_func(x) * z1 + fxy_func(x, z1);
}

void rk_step(uint idx, double *arg, double *cur, double *other, double step, double (*expr)(double, double, double))
{
	double k1 = expr(arg[idx - 1], cur[idx - 1], other[idx - 1]);
	double k2 = expr(arg[idx - 1] + step / 2, cur[idx - 1] + k1 * step / 2, other[idx - 1]);
	double k3 = expr(arg[idx - 1] + step / 2, cur[idx - 1] + k2 * step / 2, other[idx - 1]);
	double k4 = expr(arg[idx - 1] + step, cur[idx - 1] + k3 * step, other[idx - 1]);

	cur[idx] = cur[idx - 1] + step / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
}

void adams_step(uint idx, double *arg, double *cur, double *other, double step, double (*expr)(double, double, double))
{
	// Predictor step:
	cur[idx] = cur[idx - 1] + step * (1901.0 / 720 * expr(arg[idx - 1], cur[idx - 1], other[idx - 1]) -
								      1387.0 / 360 * expr(arg[idx - 2], cur[idx - 2], other[idx - 2]) +
								      109.0  / 30  * expr(arg[idx - 3], cur[idx - 3], other[idx - 3]) -
								      637.0  / 360 * expr(arg[idx - 4], cur[idx - 4], other[idx - 4]) +
								      251.0  / 720 * expr(arg[idx - 5], cur[idx - 5], other[idx - 5]));

	// Corrector step:
	cur[idx] = cur[idx - 1] + step / 720 * (251 * expr(arg[idx], cur[idx], other[idx - 1]) +
											646 * expr(arg[idx - 1], cur[idx - 1], other[idx - 1]) -
									        264 * expr(arg[idx - 2], cur[idx - 2], other[idx - 2]) +
									        106 * expr(arg[idx - 3], cur[idx - 3], other[idx - 3]) -
									        19  * expr(arg[idx - 4], cur[idx - 4], other[idx - 4]));
}