#pragma once

#include "headers.h"

/*
 * ОПИСАНИЕ:
 * Строит график функции и выводит его изображение в PNG-файл.
 *
 * ПАРАМЕТРЫ:
 * uint size - размер изображения
 * uint len - длина массивов аргументов и значений функции
 * double *arg - массив аргументов функции
 * double *res - массив значений функции
 * char* filename - имя выходного файла
 */
void plot_build(uint size, uint len, double *arg, double *res, char* filename);