#include "solver.h"

#include <cstdio>
#include <cstdlib>

void show_array(uint len, double *arr)
{
	for (uint i = 0; i < len; ++i)
		printf("%f\n", arr[i]);
}

int main(int argc, char *argv[])
{
	double eps, *y_array;

	BoundaryData data = {};

	FILE *in_file = fopen(argv[1], "r"),
		 *out_file = fopen(argv[2], "w");

	if (in_file && out_file)
	{
		fscanf_s(in_file, "%lf %lf %lf %lf %lf %d", &data.arg_a, &data.arg_b, &data.func_a, &data.func_b, &eps, &data.intervals);

		y_array = (double*)malloc(sizeof(double) * (data.intervals + 1));

		runge_error_solve(&data, -1, &y_array, eps);

		printf("len = %d\n", data.intervals);
		printf("y_array:\n");
		show_array(data.intervals + 1, y_array);

		fprintf(out_file, "%d\n", data.intervals);

		for (uint i = 0; i < data.intervals + 1; ++i)
			fprintf(out_file, "%lf\n", y_array[i]);

		free(y_array);

		fclose(out_file);
		fclose(in_file);
	}
	else
	{
		printf("File opening error!\nAborting program!\n");
		exit(1);
	}

	return 0;
}
