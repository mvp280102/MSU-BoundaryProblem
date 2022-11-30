#include "plotter.h"

/*
 * Строит график функции и выводит его изображение в PNG-файл.
 */
void plot_build(uint size, uint len, double *arg, double *res, char* filename)
{
	RGBABitmapImageReference *canvas_reference = CreateRGBABitmapImageReference();
	StringReference *error_message = CreateStringReference((wchar_t*)"Plot building error!\n", 21);

	if (DrawScatterPlot(canvas_reference, size, size, arg, len, res, len, error_message))
	{
		size_t length;
		double *png_data = ConvertToPNG(&length, canvas_reference->image);

		WriteToFile(png_data, length, filename);
		DeleteImage(canvas_reference->image);
	}
	else
	{
		printf("%ls\nAborting program!\n", error_message->string);
		exit(1);
	}

	FreeAllocations();
}