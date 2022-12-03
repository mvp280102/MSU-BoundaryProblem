#include "solver.h"


//// ФУНКЦИИ СИСТЕМЫ ДИФФЕРЕНЦИАЛЬНЫХ УРАВНЕНИЙ, К КОТОРОЙ СВОДИТСЯ ИСХОДНАЯ КРАЕВАЯ ЗАДАЧА:

// Возвращает значение первого уравнения системы на текущем шаге.
double right_expr_z1(double x, double z1, double z2)
{
	return z2;
}

// Возвращает значение второго уравнения системы на текущем шаге.
double right_expr_z2(double x, double z2, double z1)
{
	return -1 * px_func(x, z1) * z2 - qx_func(x, z1) * z1 + fxy_func(x, z1);
}


//// ФУНКЦИИ ЧИСЛЕННЫХ МЕТОДОВ РЕШЕНИЯ СИСТЕМЫ ДИФФЕРЕНЦИАЛЬНЫХ УРАВНЕНИЙ:

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

	for (uint i = 0; i < ACC_ORDER; ++i)
		predictor_value += predictor_cfs[i] * expr(arg[idx - i - 1], cur[idx - i - 1], other[idx - i - 1]) * (i & 1 ? -1 : 1);

	cur[idx] = cur[idx - 1] + step * predictor_value;

	for (uint i = 0; i < ACC_ORDER; ++i)
		corrector_value += corrector_cfs[i] * expr(arg[idx - i], cur[idx - i], other[!i ? idx - 1 : idx - i]) * (!i || i & 1 ? 1 : -1);

	cur[idx] = cur[idx - 1] + step * corrector_value;
}

double rk_adams_solve(BoundaryData *data, double der_a, double *res)
{
	auto *x_array = (double*)malloc(sizeof(double) * (data->intervals + 1));
	auto step = step_fill(data->intervals + 1, x_array, data->arg_a, data->arg_b);

	auto *z_array = (double*)malloc(sizeof(double) * (data->intervals + 1));
	z_array[0] = der_a;

	res[0] = data->func_a;

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


//// ФУНКЦИИ ЧИСЛЕННЫХ МЕТОДОВ РЕШЕНИЯ ВСПОМОГАТЕЛЬНОГО УРАВНЕНИЯ ДЛЯ ВЫЧИСЛЕНИЯ ЗНАЧЕНИЯ ПРОИЗВОДНОЙ НА ЛЕВОМ КОНЦЕ:

// Возвращает численное значение производной переданной функции в заданной точке.
double numerical_derivative(BoundaryData *data, double arg, double *res, double (*func)(BoundaryData*, double, double*))
{
	double step = (data->arg_b - data->arg_a) / data->intervals;
	return (func(data, arg + step, res) - func(data, arg - step, res)) / (2 * step);
}

// Выполняет шаг метода Ньютона решения уравнения вида f(x) = 0.
// Возвращает значение следующего приближения аргумента функции.
double newton_step(BoundaryData *data, double arg, double *res, double (*func)(BoundaryData*, double, double*))
{
	return arg - func(data, arg, res) / numerical_derivative(data, arg, res, func);
}

// Реализует метод Ньютона решения уравнения вида f(x) = 0 с заданной точностью.
// После выполнения массив значений функций заполнен в соответствии с данными задачи.
void newton_solve(BoundaryData *data, double arg, double *res, double eps)
{
	double der_a, diff = rk_adams_solve(data, arg, res);

	while (fabs(diff) > eps)
	{
		der_a = newton_step(data, diff, res, rk_adams_solve);
		diff = rk_adams_solve(data, der_a, res);
	}
}


//// ФУНКЦИИ МЕТОДОВ АПОСТЕРИОРНОЙ ОЦЕНКИ ПОГРЕШНОСТИ ПО ПРАВИЛУ РУНГЕ:

//Возвращает евклидово расстояние между двумя векторами.
// Подразумевается, что второй массив в два раза длиннее первого.
// При вычислениях используется каждая точка первого массива и каждая вторая точка второго.
double vector_distance(uint len, double *y1, double *y2)
{
	double res = 0;

	for (uint i = 0; i < len; ++i)
		res += pow(y1[i] - y2[2 * i], 2);

	return sqrt(res) / len;
}

double runge_error_step(BoundaryData *data_n, double arg, double *res_n, double *res_2n, double eps)
{
	BoundaryData data_2n = *data_n;
	data_2n.intervals *= 2;

	newton_solve(data_n, arg, res_n, eps);
	newton_solve(&data_2n, arg, res_2n, eps);

	return vector_distance(data_n->intervals + 1, res_n, res_2n) / ((1 << ACC_ORDER) - 1);
}

void runge_error_solve(BoundaryData *data, double arg, double **res, double eps)
{
	auto *res_doubled = (double*)malloc(sizeof(double) * (data->intervals * 2 + 1));

	while (runge_error_step(data, arg, *res, res_doubled, eps) > eps)
	{
		data->intervals *= 2;
		*res = (double*)realloc(*res, sizeof(double) * (data->intervals + 1));
		memcpy(*res, res_doubled, sizeof(double) * data->intervals + 1);
		res_doubled = (double*)realloc(res_doubled, sizeof(double) * (data->intervals * 2 + 1));
	}

	free(res_doubled);
}


//// ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ:

// Заполняет массив значениями с равномерным шагом в заданном диапазоне.
// Возвращает разность между двумя последовательными значениями в массиве.
double step_fill(uint len, double *arr, double start, double stop)
{
	double step = (stop - start) / (len - 1);

	arr[0] = start;
	arr[len - 1] = stop;

	for (uint i = 1; i < len - 1; ++i)
		arr[i] = arr[i - 1] + step;

	return step;
}