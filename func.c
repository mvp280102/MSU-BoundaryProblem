#include "func.h"


// Возвращает значение точного решения краевой задачи.
double yx_func(double x)
{
    return 0;
}

// Возвращает значение коэффициента p(x).
double px_func(double x)
{
    return -1 * (x + 4) / (2 * x);
}

// Возвращает значение коэффициента q(x).
double qx_func(double x)
{
    return (1 + 4 / x) / (2 * x);
}

// Возвращает значение правой части f(x, y).
double fxy_func(double x, double y)
{
    return pow(x, 3) / (2 * x);
}