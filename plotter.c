#include "plotter.h"


// Строит графики одной или двух функций на одной координатной плоскости, и выводит изображение в PNG-файл.
void plot_build(size_t size, char* filename, size_t len, double *arg, double *func1, double *func2)
{
    size_t series_len = 1;

    StringReference *error_message = CreateStringReference((wchar_t*)"Plot building error!\n", 21);
	RGBABitmapImageReference *canvas_reference = CreateRGBABitmapImageReference();
    ScatterPlotSeries **series, *series1, *series2;
    ScatterPlotSettings *settings;

    series1 = GetDefaultScatterPlotSeriesSettings();
    series1->xs = arg;
    series1->xsLength = len;
    series1->ys = func1;
    series1->ysLength = len;
    series1->linearInterpolation = false;
    series1->color = CreateRGBColor(1, 0, 0);

    if (func2 != NULL)
    {
        series_len = 2;

        series2 = GetDefaultScatterPlotSeriesSettings();
        series2->xs = arg;
        series2->xsLength = len;
        series2->ys = func2;
        series2->ysLength = len;
        series2->linearInterpolation = false;
        series2->color = CreateRGBColor(0, 0, 1);
    }

    settings =GetDefaultScatterPlotSettings();
    settings->width = (double)size;
    settings->height = (double)size;
    settings->autoPadding = true;
    settings->autoBoundaries = true;

    series = (ScatterPlotSeries**)malloc(sizeof(ScatterPlotSeries*) * series_len);
    series[0] = series1;

    if (series_len == 2)
        series[1] = series2;

    settings->scatterPlotSeries = series;
    settings->scatterPlotSeriesLength = series_len;

	if (DrawScatterPlotFromSettings(canvas_reference, settings, error_message))
	{
		size_t length;
		double *png_data = ConvertToPNG(&length, canvas_reference->image);

		WriteToFile(png_data, length, filename);
		DeleteImage(canvas_reference->image);
	}
	else
	{
		printf("%ls\nAborting application!\n", error_message->string);
		exit(1);
	}

    free(series);
	FreeAllocations();
}