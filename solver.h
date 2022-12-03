#pragma once

#include "pxy_func.h"
#include "qxy_func.h"
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

/*
 * ОПИСАНИЕ:
 * Возвращает численное значение производной переданной функции в заданной точке.
 *
 * ПАРАМЕТРЫ:
 * BoundaryData *data - данные для решения основной задачи
 * double arg - значение аргумента функции
 * double *res - массив значений функции
 * double (*func) - указатель на функцию, для которой вычисляется производная
 */
double numerical_derivative(BoundaryData *data, double arg, double *res, double (*func)(BoundaryData*, double, double*));

/*
 * ОПИСАНИЕ:
 * Выполняет шаг метода Ньютона решения уравнения вида f(x) = 0.
 * Возвращает значение следующего приближения аргумента функции.
 *
 * ПАРАМЕТРЫ:
 * BoundaryData *data - данные для решения основной задачи
 * double arg - текущее значение приближения аргумента функции
 * double *res - массив значений функции
 * double (*func) - указатель на функцию, для которой вычисляется производная
 */
double newton_step(BoundaryData *data, double arg, double *res, double (*func)(BoundaryData*, double, double*));

/*
 * ОПИСАНИЕ:
 * Реализует метод Ньютона решения уравнения вида f(x) = 0 с заданной точностью.
 * После выполнения массив значений функций заполнен в соответствии с данными задачи.
 *
 * ПАРАМЕТРЫ:
 * BoundaryData *data - данные для решения основной задачи
 * double arg - начальное приближение аргумента функции
 * double *res - массив значений функции
 * double eps - точность
 */
void newton_solve(BoundaryData *data, double arg, double *res, double eps);


//// ФУНКЦИИ МЕТОДОВ АПОСТЕРИОРНОЙ ОЦЕНКИ ПОГРЕШНОСТИ ПО ПРАВИЛУ РУНГЕ:

/*
 * ОПИСАНИЕ:
 * Возвращает евклидово расстояние между двумя векторами.
 * Подразумевается, что второй массив в два раза длиннее первого.
 * При вычислениях используется каждая точка первого массива и каждая вторая точка второго.
 *
 * ПАРАМЕТРЫ:
 * uint len - длина первого массива
 * double *y1 - первый массив
 * double *y2 - второй массив
 */
double vector_distance(uint len, double *y1, double *y2);
double runge_error_step(BoundaryData *data_n, double arg, double *res_n, double *res_2n, double eps);
void runge_error_solve(BoundaryData *data, double arg, double **res, double eps);


//// ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ:

/*
 * ОПИСАНИЕ:
 * Заполняет массив значениями с равномерным шагом в заданном диапазоне.
 * Возвращает разность между двумя последовательными значениями в массиве.
 *
 * ПАРАМЕТРЫ:
 * uint len - длина заполняемого массива
 * double *arr - заполняемый массив
 * double start - левый конец диапазона
 * double stop - правый конец правый конец
 */
double step_fill(uint len, double *arr, double start, double stop);