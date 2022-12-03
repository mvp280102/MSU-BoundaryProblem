### ПРИЛОЖЕНИЕ, РЕШАЮЩЕЕ КРАЕВЫЕ ЗАДАЧИ ДЛЯ ЛИНЕЙНЫХ ДИФФЕРЕНЦИАЛЬНЫХ УРАВНЕНИЙ 2-ГО ПОРЯДКА С ПЕРЕМЕННЫМИ КОЭФФИЦИЕНТАМИ

--- 

#### СОДЕРЖАНИЕ:
1. Постановка задачи.
2. Теоретическое описание алгоритма.
3. Взаимодействие с приложением.

---

#### 1. ПОСТАНОВКА ЗАДАЧИ

Дана краевая задача для линейного дифференциального уравнения 2-го порядка с переменными коэффициентами вида
y'' + p(x, y) * y' + q(x, y) * y = f(x, y). Даны также координаты отрезка [x0; x1], значения искомой функции на его
концах вида y(x0) = y0 и y(x1) = y1, желаемая точность решения и начальное количество узлов. Требуется указать значения
неизвестной функции y в узлах отрезка, количество которых позволяет достичь заданной точности, а также построить график
этого приближенного решения.

---

#### 2. ТЕОРЕТИЧЕСКОЕ ОПИСАНИЕ АЛГОРИТМА

Применяется метод сведения линейного дифференциального уравнения второго порядка к системе обыкновенных линейных
дифференциальных уравнений первого порядка путем введения замен y = z1 и y' = z2. Полученная система состоит из двух
уравнений вида z1' = z2 и z2' = -p(x, z1) * z2 - q(x, z1) * z1 + f(x, z1). Для решения этой системы на отрезке [x0; x1]
необходимы начальные условия вида z1(x0) = y1 и z2(x0) = t. Следует отметить, что z2(x0) является значением производной
на левом конце отрезка, что в данном случае не известно. Поэтому в процессе решения исходной задачи это значение
варьируется для того, чтобы значение неизвестной функции на правом конце совпадало с заданным.

Для решения полученной системы дифференциальных уравнений используется неявный метод Адамса-Мултона 5-го порядка. Так
как для вычисления значения каждой неизвестной функции в каждой текущей точке он требует известных значений функции в 5
предыдущих точках, а в самом начале известно лишь одно значение - в точке 1, до запуска метода Адамса-Мултона значения
функции в точках 2, 3, 4 и 5 вычисляются с помощью классического метода Рунге-Кутты 4-го порядка. Далее запускается
собственно метод Адамса-Мултона. Так как система содержит два связанных уравнения, для вычисления значения неизвестной
функции в каждом уравнении системы на текущем шаге используются вычисленные значения неизвестных функций обоих уравнений
на предыдущих шагах. Таким образом, результатом решения исходного уравнения 2-го порядка являются значения неизвестной
функции z1, т. к. замена имеет вид y = z1.

Вычисление правильного значения производной на левом конце отрезка производится с помощью метода Ньютона. В самом начале
задается некое произвольное значение производной - начальное приближение, с использованием которого решается исходная
задача. Далее находится разность между вычисляемым и заданным изначально значением неизвестной функции на правом конце
отрезка. Разность считается функцией, зависящей от параметра t - искомого значения производной. Метод Ньютона
используется для минимизации этой функции, то есть такого значения параметра, при котором модуль разности между
вычисляемым и данным изначально значением искомой функции не будет превосходить заданного значения погрешности.

Метод Ньютона требует знания значения производной функции в каждой текущей точке. Так как для вычисления значения этой
функции необходимо решать систему уравнений вышеописанным методом, и в явном виде указать ее затруднительно, ее
производная вычисляется с помощью методов численного дифференцирования.

После того как правильное значение производной на левом конце вычислено, количество узлов отрезка удваивается в
соответствии с правилом Рунге апостериорной оценки погрешности численного решения обыкновенного дифференциального
уравнения при p = 5 (порядок точности используемого метода Адамса-Мултона). Все вышеописанные действия повторяются на
отрезке, имеющем вдвое большее количество узлов. Затем по правилу Рунге вычисляется погрешность решения при удвоенном
числе узлов. Уравнение решается на удваивающемся на каждом шаге количестве узлов, пока вычисляемая погрешность Рунге не
станет меньше заданного значения.

Более подробные описания функций, используемых в исходном коде приложения, содержатся в его заголовочных файлах.

---

#### 3. ВЗАИМОДЕЙСТВИЕ С ПРИЛОЖЕНИЕМ

Приложение является консольным, ожидая текстовый ввод из файла и производя текстовый и графический вывод в два разных
файла. Имена этих трех файлов должны быть переданы приложению при запуске в качестве аргументов. Входной файл должен
быть предварительно создан. Выходные создаются при необходимости. Во входом файле через пробел задаются: левый и правый
концы отрезка x0 и x1, значения неизвестной функции на левом и правом концах отрезка y0 и y1, желаемая точность и
начальное количество узлов на отрезке. Все числа кроме последнего могут быть вещественными, последнее же - целое
положительное.

Функции p(x, y), q(x, y) и f(x, y), задающие исходное линейное дифференциальное уравнение 2-го порядка, определяются в
файлах pxy_func.(h|cpp), qxy_func.(h|cpp) и fxy_func.(h|cpp) соответственно. Эти файлы включаются в исходный код
приложения до компиляции.

После достижения погрешностью решения задачи заданного значения в выходной текстовый файл выводится начальное количество
узлов на отрезке, а затем значения искомой функции в этих узлах, каждое с новой строки. Также приложение строит график
полученного решения и выводит его в выходной графический файл. Построение графика производится с помощью библиотеки с
открытым исходным кодом pbPlots. Библиотека представлена исключительно исходными файлами, уже содержащимися в проекте,
никаких дополнительных заранее скомпилированных файлов для запуска приложения не требуется.