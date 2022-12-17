#pragma once

#include "headers.h"


/*
 * ОПИСАНИЕ:
 * Возвращает значение точного решения краевой задачи.
 *
 * ПАРАМЕТРЫ:
 * double x - аргумент
 */
double yx_func(double x);

/*
 * ОПИСАНИЕ:
 * Возвращает значение коэффициента p(x).
 *
 * ПАРАМЕТРЫ:
 * double x - аргумент
 */
double px_func(double x);

/*
 * ОПИСАНИЕ:
 * Возвращает значение коэффициента q(x).
 *
 * ПАРАМЕТРЫ:
 * double x - аргумент
 */
double qx_func(double x);

/*
 * ОПИСАНИЕ:
 * Возвращает значение правой части f(x, y).
 *
 * ПАРАМЕТРЫ:
 * double x - аргумент x
 * double y - аргумент y
 */
double fxy_func(double x, double y);