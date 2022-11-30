#include "fxy_func.h"

double fxy_func(double x, double y)
{
	return 6 * cos(x) * sin(x) * (cos(x) - y * sin(x));
}