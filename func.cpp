#include "func.h"

double yx_func(double x)
{
    return pow(sin(x), 3);
}

double px_func(double x)
{
	return -1 * (x + 4) / (2 * x);
}

double qx_func(double x)
{
	return (1 + 4 / x) / (2 * x);
}

double fxy_func(double x, double y)
{
	return pow(x, 3) / (2 * x);
}