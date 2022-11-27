#include "solver.h"

double vector_distance(uint len, double *y1, double *y2)
{
	double res = 0;

	for (uint i = 0, j = 0; i < len; i += 2, ++j)
		res += pow(y1[i] - y2[j], 2);

	return sqrt(res) / len;
}

double numerical_derivative(BoundaryData *data, double arg, double *res, double (*func)(BoundaryData*, double, double*))
{
	double step = (data->arg_b - data->arg_a) / data->intervals;
	return (func(data, arg + step, res) - func(data, arg - step, res)) / (2 * step);
}

double newton_step(BoundaryData *data, double arg, double *res, double (*func)(BoundaryData*, double, double*))
{
	return arg - func(data, arg, res) / numerical_derivative(data, arg, res, func);
}

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
	double predictor_cfs[ACC_ORDER] = { 1901.0 / 720, 1387.0 / 360, 109.0 / 30, 637.0 / 360, 251.0 / 720 },
		   corrector_cfs[ACC_ORDER] = { 251.0 / 720, 646.0 / 720, 264.0 / 720, 106.0 / 720, 19.0 / 720 };

	double predictor_value = 0,
		   corrector_value = 0;

	// Predictor step:

	for (uint i = 0; i < ACC_ORDER; ++i)
		predictor_value += predictor_cfs[i] * expr(arg[idx - i - 1], cur[idx - i - 1], other[idx - i - 1]) * (i & 1 ? -1 : 1);

	cur[idx] = cur[idx - 1] + step * predictor_value;

	// Corrector step:

	for (uint i = 0; i < ACC_ORDER; ++i)
		corrector_value += corrector_cfs[i] * expr(arg[idx - i], cur[idx - i], other[!i ? idx - 1 : idx - i]) * (!i || i & 1 ? 1 : -1);

	cur[idx] = cur[idx - 1] + step * corrector_value;
}

double rk_adams_solve(BoundaryData *data, double der_a, double *res)
{
	double step = (data->arg_b - data->arg_a) / data->intervals;

	auto *x_array = (double*)malloc(sizeof(double) * (data->intervals + 1));
	x_array[0] = data->arg_a;
	x_array[data->intervals] = data->arg_b;

	auto *z_array = (double*)malloc(sizeof(double) * (data->intervals + 1));
	z_array[0] = der_a;

	res[0] = data->func_a;

	for (uint i = 1; i < data->intervals; ++i)
		x_array[i] = x_array[i - 1] + step;

	for (uint i = 1; i < ACC_ORDER; ++i)
	{
		rk_step(i, x_array, res, z_array, step, right_expr_z1);
		rk_step(i, x_array, z_array, res, step, right_expr_z2);
	}

	for (uint i = ACC_ORDER; i < data->intervals + 1; ++i)
	{
		adams_step(i, x_array, res, z_array, step, right_expr_z1);
		adams_step(i, x_array, z_array, res, step, right_expr_z2);
	}

	free(z_array);
	free(x_array);

	return res[data->intervals] - data->func_b;
}