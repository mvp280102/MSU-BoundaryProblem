#pragma once

#include "func.h"

#define ACC_ORDER 5

// Начальные данные для решения краевой задачи.
typedef struct BoundaryData
{
	size_t intervals;								// Количество узлов отрезка.

	double arg_a;								// Координата левого конца отрезка.
	double arg_b;								// Координата правого конца отрезка.
	double func_a;								// Значение функции на левом конце отрезка.
	double func_b;								// Значение функции на правом конце отрезка.
} BoundaryData;


//// ФУНКЦИИ СИСТЕМЫ ДИФФЕРЕНЦИАЛЬНЫХ УРАВНЕНИЙ (подробнее - readme.md), К КОТОРОЙ СВОДИТСЯ ИСХОДНАЯ КРАЕВАЯ ЗАДАЧА:

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


//// ФУНКЦИИ ЧИСЛЕННЫХ МЕТОДОВ РЕШЕНИЯ СИСТЕМЫ ДИФФЕРЕНЦИАЛЬНЫХ УРАВНЕНИЙ (подробнее - readme.md):

/*
 * ОПИСАНИЕ:
 * Выполняет шаг метода Рунге-Кутты 4-го порядка.
 * После выполнения элемент массива cur с индексом idx содержит вычисленное значение производной.
 *
 * ПАРАМЕТРЫ:
 * uint idx - индекс узла в массиве узлов и искомого значения уравнения в этом узле в массиве значений уравнения
 * double *arg - массив узлов
 * double *cur - массив значений решаемого уравнения системы
 * double *other - массив значений другого уравнения системы
 * double step - разность двух последовательных узлов
 * double (*expr) - правая часть решаемого уравнения
 */
void rk_step(size_t idx, double *arg, double *cur, double *other, double step, double (*expr)(double, double, double));

/*
 * ОПИСАНИЕ:
 * Выполняет шаг метода Адамса-Мултона 5-го порядка.
 * После выполнения элемент массива cur с индексом idx содержит вычисленное значение производной.
 *
 * ПАРАМЕТРЫ:
 * uint idx - индекс узла в массиве узлов и искомого значения уравнения в этом узле в массиве значений уравнения
 * double *arg - массив узлов
 * double *cur - массив значений решаемого уравнения системы
 * double *other - массив значений другого уравнения системы
 * double step - разность двух последовательных узлов
 * double (*expr) - правая часть решаемого уравнения
 */
void adams_step(size_t idx, double *arg, double *cur, double *other, double step, double (*expr)(double, double, double));

/*
 * ОПИСАНИЕ:
 * Решает систему уравнений методами Рунге-Кутты и Адамса-Мултона при заданном значении производной на левом конце.
 * Метод Рунге-Кутты используется для вычисления значений уравнений в узлах 2, 3, 4.
 * Для вычисления значений уравнений в остальных узлах используется метод Адамса-Мултона.
 * После выполнения массив значений функции res заполнен в соответствии с данными задачи и заданным значением производной.
 * Возвращает разность между полученным и заданным значениями искомой функции на правом конце отрезка.
 *
 * ПАРАМЕТРЫ:
 * BoundaryData *data - данные для решения основной задачи
 * double der_a - значение производной искомой функции на левом конце отрезка
 * double *res - массив значений искомой функции
 */
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
 * После выполнения массив значений функции res заполнен в соответствии с данными задачи.
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
 * Подразумевается, что второй вектор в два раза длиннее первого.
 * При вычислениях используется каждая координата первого вектора и каждая вторая координата второго.
 *
 * ПАРАМЕТРЫ:
 * uint len - длина первого вектора
 * double *y1 - первый вектор
 * double *y2 - второй вектор
 */
double vector_distance(size_t len, double *y1, double *y2);

/*
 * ОПИСАНИЕ:
 * Выполняет шаг метода Рунге апостериорной оценки погрешности.
 * Возвращает значение погрешности при количестве узлов, в два раза большем исходного.
 *
 * ПАРАМЕТРЫ:
 * BoundaryData *data - данные для решения основной задачи
 * double arg - начальное приближение для метода Ньютона
 * double *res_n - массив значений искомой функции при исходном количестве узлов
 * double *res_2n - массив значений искомой функции при удвоенном количестве узлов
 * double eps - точность
 */
double runge_error_step(BoundaryData *data_n, double arg, double *res_n, double *res_2n, double eps);

/*
 * Реализует метод Рунге апостериорной оценки погрешности.
 * После выполнения массив значений искомой функции res увеличен до размера, равного количеству узлов,
 * необходимому для достижения заданной точности, и заполнен в соответствии с данными задачи.
 *
 * ПАРАМЕТРЫ:
 * BoundaryData *data - данные для решения основной задачи
 * double arg - начальное приближение для метода Ньютона
 * double **res - указатель на массив значений искомой функции
 * double eps - точность
 */
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
double step_fill(size_t len, double *arr, double start, double stop);

/*
 * ОПИСАНИЕ:
 * Заполняет массив значениями функции на заданных в другом массиве значениях аргумента.
 * Подразумевается, что массивы значений аргумента и функции имеют одинаковую длину.
 * Если вместо указателя на функцию передается NULL, заполнения не происходит.
 *
 * ПАРАМЕТРЫ:
 * uint len - длина используемых массивов
 * double *arg - массив значений аргумента
 * double *res - массив значений функции
 * double (*func) - указатель на функцию
 */
void func_fill(size_t len, double *arg, double *res, double (*func)(double));