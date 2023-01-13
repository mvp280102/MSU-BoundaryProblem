#include "solver.h"
#include "plotter.h"


int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        printf(ARGS_ERROR);
        printf("Aborting application!\n");
        exit(1);
    }

    size_t init_intervals;                      // Начальное количество узлов отрезка.

    double eps,                                 // Точность.
           *x_array,                            // Массив узлов.
           *y_approx,                           // Массив приближенных значений искомой функции.
           *y_exact;                            // Массив точных значений искомой функции.

    BoundaryData data = {};                     // Начальные данные для решения задачи.

    FILE *in_file = fopen(argv[1], "r"),        // Входной файл.
         *out_file = fopen(argv[2], "w");       // Выходной файл.

    if (!in_file || !out_file)
    {
        printf(FILES_ERROR);
        printf("Aborting application!\n");
        exit(1);
    }

    fscanf_s(in_file, "%lf %lf %lf %lf %lf %d", &data.arg_a, &data.arg_b, &data.func_a, &data.func_b, &eps, &data.intervals);

    init_intervals = data.intervals;

    y_approx = (double*)malloc(sizeof(double) * (data.intervals + 1));

    runge_error_solve(&data, -1, &y_approx, eps);

    x_array = (double*)malloc(sizeof(double) * (data.intervals + 1));
    step_fill(data.intervals + 1, x_array, data.arg_a, data.arg_b);

    y_exact = (double*)malloc(sizeof(double) * (data.intervals + 1));
    func_fill(data.intervals + 1, x_array, y_exact, NULL);

    fprintf(out_file, "%zu\n", init_intervals);

    for (size_t i = 0; i < data.intervals + 1; i += (data.intervals / init_intervals))
        fprintf(out_file, "%lf\n", y_approx[i]);

    plot_build(1024, argv[3], data.intervals + 1, x_array, y_approx, NULL);

    free(x_array);
    free(y_exact);
    free(y_approx);

    fclose(in_file);
    fclose(out_file);

    return 0;
}
