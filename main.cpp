#include "solver.h"

#include <cstdio>
#include <cstdlib>

void show_array(uint len, double *arr)
{
	for (uint i = 0; i < len; ++i)
		printf("%f ", arr[i]);

	printf("\n");
}

int main(int argc, char *argv[])
{
	double eps, der_a, diff, *y_array;

	BoundaryData data = {};

	FILE *in_file = fopen(argv[1], "r"),
		 *out_file = fopen(argv[2], "w");

	if (in_file && out_file)
	{
		fscanf_s(in_file, "%lf %lf %lf %lf %lf %d", &data.arg_a, &data.arg_b, &data.func_a, &data.func_b, &eps, &data.intervals);

		y_array = (double*)malloc(sizeof(double) * (data.intervals + 1));

		diff = boundary_solve(data, -1, y_array);

		while (fabs(diff) > eps)
		{
			der_a = newton_step(data, diff, y_array, boundary_solve);
			diff = boundary_solve(data, der_a, y_array);
		}

		printf("der_a = %lf\n\n", der_a);

		printf("y_array:\n");
		show_array(data.intervals + 1, y_array);

		fprintf(out_file, "%lf %lf %d\n", data.arg_a, data.arg_b, data.intervals);

		for (uint i = 0; i < data.intervals + 1; ++i)
			fprintf(out_file, "%lf\n", y_array[i]);

		system("python plot.py");

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
