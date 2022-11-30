#include "boundary_solver.h"
#include "plot_builder.h"

int main(int argc, char *argv[])
{
	double eps, *x_array, *y_array;

	BoundaryData data = {};

	FILE *in_file = fopen(argv[1], "r"),
		 *out_file = fopen(argv[2], "w");

	if (in_file && out_file)
	{
		fscanf_s(in_file, "%lf %lf %lf %lf %lf %d", &data.arg_a, &data.arg_b, &data.func_a, &data.func_b, &eps, &data.intervals);

		y_array = (double*)malloc(sizeof(double) * (data.intervals + 1));

		runge_error_solve(&data, -1, &y_array, eps);

		x_array = (double*)malloc(sizeof(double) * (data.intervals + 1));
		step_fill(data.intervals + 1, x_array, data.arg_a, data.arg_b);

		fprintf(out_file, "%d\n", data.intervals);

		for (uint i = 0; i < data.intervals + 1; ++i)
			fprintf(out_file, "%lf\n", y_array[i]);

		plot_build(1024, data.intervals + 1, x_array, y_array, (char*)"plot.png");

		free(x_array);
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
