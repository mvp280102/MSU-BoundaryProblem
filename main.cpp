#include "px_func.h"
#include "qx_func.h"
#include "fxy_func.h"

#include <cstdio>
#include <cstdlib>

typedef unsigned int uint;

int main(int argc, char *argv[])
{
	double x0, x1, y0, y1, h, eps;

	FILE *in_file = fopen(argv[1], "r"),
			*out_file = fopen(argv[2], "w");

	if (in_file != nullptr)
	{
		fscanf_s(in_file, "%lf %lf %lf %lf %lf %lf", &x0, &x1, &y0, &y1, &h, &eps);
	}
	else
	{
		printf("Input file opening error!\nAborting program!\n");
		exit(1);
	}

	return 0;
}
