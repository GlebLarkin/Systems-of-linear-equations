# Библиотека для решения систем линейных уравнений (СЛАУ)

Библиотека предоставляет инструменты для решения систем линейных уравнений (СЛАУ). На данный момент реализованы:
- Алгоритм прогонки для трехдиагональной матрицы
- :(

## Как это использовать то?
### Пример для метода прогонки для трехдиагональной матрицы:
```cxx
#include "Thomas_algorithm.h"

int main()
{
    // Задаем диагонали (пример для матрицы 4x4):
    size_t matrix_size = 4;

    std::vector<double> lower_diag = {3, 2, 0};      // Нижняя диагональ (N-1 элементов)
    std::vector<double> main_diag = {5, 8, -10, 1};  // Главная диагональ (N элементов)
    std::vector<double> upper_diag = {1, 2, 6};      // Верхняя диагональ (N-1 элементов)

    // Тогда получим:
    // { 5   1   0   0 }
    // { 3   8   2   0 }
    // { 0   2  -10  6 }
    // { 0   0   0   1 }


    // Создаем трехдиагональную матрицу
    Tridiagonal_matrix<double> matrix(upper_diag, main_diag, lower_diag, matrix_size);


    // Вектор свободных членов (должен иметь размер = matrix_size)
    std::vector<double> free_col = {3, 8, 5, 1};


    // Решаем систему
    Thomas_algorithm<double> solver(matrix, free_col);


    // Получаем решение
    std::vector<double> solution = solver.Get_solution();

    return 0;
}
```
