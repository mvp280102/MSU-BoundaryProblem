#pragma once

#include "px_func.h"
#include "qx_func.h"
#include "fxy_func.h"

#define ACC_ORDER 5

// Начальные данные для решения краевой задачи.
struct BoundaryData
{
	uint intervals;								// Количество точек деления отрезка.

	double arg_a;								// Координата левого конца отрезка.
	double arg_b;								// Координата правого конца отрезка.
	double func_a;								// Значение функции на левом конце отрезка.
	double func_b;								// Значение функции на правом конце отрезка.
};


//// ФУНКЦИИ СИСТЕМЫ ДИФФЕРЕНЦИАЛЬНЫХ УРАВНЕНИЙ, К КОТОРОЙ СВОДИТСЯ ИСХОДНАЯ КРАЕВАЯ ЗАДАЧА:

/*
 * ОПИСАНИЕ:
 * Возвращает значение первого уравнения системы на текущем шаге.
 *
 * ПАРАМЕТРЫ:
 * double x - значение аргумента
 * double z1 - значение первого уравнения системы на предыдущем шаге
 * double z2 - значение второго уравнения системы на предыдущем шаге
 */
double right_expr_z1(double x, double z1, double z2);

/*
 * ОПИСАНИЕ:
 * Возвращает значение второго уравнения системы на текущем шаге.
 *
 * ПАРАМЕТРЫ:
 * double x - значение аргумента
 * double z2 - значение второго уравнения системы на предыдущем шаге
 * double z1 - значение первого уравнения системы на предыдущем шаге
 */
double right_expr_z2(double x, double z2, double z1);


//// ФУНКЦИИ ЧИСЛЕННЫХ МЕТОДОВ РЕШЕНИЯ СИСТЕМЫ ДИФФЕРЕНЦИАЛЬНЫХ УРАВНЕНИЙ:
void rk_step(uint idx, double *arg, double *cur, double *other, double step, double (*expr)(double, double, double));
void adams_step(uint idx, double *arg, double *cur, double *other, double step, double (*expr)(double, double, double));
double rk_adams_solve(BoundaryData *data, double der_a, double *res);


//// ФУНКЦИИ ЧИСЛЕННЫХ МЕТОДОВ РЕШЕНИЯ ВСПОМОГАТЕЛЬНОГО УРАВНЕНИЯ ДЛЯ ВЫЧИСЛЕНИЯ ЗНАЧЕНИЯ ПРОИЗВОДНОЙ НА ЛЕВОМ КОНЦЕ:
double numerical_derivative(BoundaryData *data, double arg, double *res, double (*func)(BoundaryData*, double, double*));
double newton_step(BoundaryData *data, double arg, double *res, double (*func)(BoundaryData*, double, double*));
void newton_solve(BoundaryData *data, double arg, double *res, double eps);


//// ФУНКЦИИ МЕТОДОВ АПОСТЕРИОРНОЙ ОЦЕНКИ ПОГРЕШНОСТИ ПО ПРАВИЛУ РУНГЕ:
double vector_distance(uint len, double *y1, double *y2);
double runge_error_step(BoundaryData *data_n, double arg, double *res_n, double *res_2n, double eps);
void runge_error_solve(BoundaryData *data, double arg, double **res, double eps);


//// ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ:
double step_fill(uint len, double *arr, double start, double stop);
void array_output(uint len, double *arr);
