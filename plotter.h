#pragma once

#include "headers.h"

/*
 * ОПИСАНИЕ:
 * Строит график функции и выводит его изображение в PNG-файл.
 *
 * ПАРАМЕТРЫ:
 * uint size - размер изображения
 * char* filename - имя выходного файла
 * uint len - длина массивов значений аргумента и функции
 * double *arg - массив значений аргумента
 * double *res - массив значений функции
 */
void plot_build(uint size, char* filename, uint len, double *arg, double *res);