#include "solver.h"

// z - результат второго уравнения
double total_func(double x, double y, double z)
{
	return -1 * px_func(x) * y - qx_func(x) * z + fxy_func(x, z);
}

void solving_step(uint index, double *x, double *y, double *z, double h)
{
	// Шаг предиктора (метод Адамса-Башфорта).
	y[index] = y[index - 1] + h * (1901.0 / 720 * total_func(x[index - 1], y[index - 1], z[index - 1]) -
								   1387.0 / 360 * total_func(x[index - 2], y[index - 2], z[index - 2]) +
								   109.0 / 30 * total_func(x[index - 3], y[index - 3], z[index - 3]) -
								   637.0 / 360 * total_func(x[index - 4], y[index - 4], z[index - 4]) +
								   251.0 / 720 * total_func(x[index - 5], y[index - 5], z[index - 5]));

	// Шаг корректора (методы Адамса-Мултона).
	y[index] = y[index - 1] + h / 720 * (251 * total_func(x[index], y[index], z[index]) +
								         646 * total_func(x[index - 1], y[index - 1], z[index - 1]) -
										 264 * total_func(x[index - 2], y[index - 2], z[index - 2]) +
										 106 * total_func(x[index - 3], y[index - 3], z[index - 3]) -
										 19 * total_func(x[index - 4], y[index - 4], z[index - 4]));
}