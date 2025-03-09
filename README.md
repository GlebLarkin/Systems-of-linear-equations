# Библиотека для решения систем линейных уравнений (СЛАУ)

Библиотека предоставляет инструменты для решения систем линейных уравнений (СЛАУ). На данный момент реализованы:
- Перегрузка основынх операторов для работы с векторами:
    - Сложение/вычитание векторов
    - Умножение вектора на число
    - Скалярное произведение векоров
    - Умножение матрицы на вектор
- Класс плотной матрицы (она умеет перемножаться на вектор)
- Класс CSR матрицы (она тоже умеет перемножаться на вектор, но лучше)
- Класс трехдиагональной матрицы (она почти ничего не умеет)
- Алгоритм прогонки для трехдиагональной матрицы
- QR разложение для CSR матрицы и решение СЛАУ с его помощью
- Итерационные методы решения СЛАУ:
    - Метод простых итераций (в том числе с помощью чебышёвского ускорения)
    - Метод Гаусса-Зейнделя
    - Метод Якоби (не используйте его, пожалуйста)
- :(

## Как это использовать то?
### Пример для работы с векторами:
```cxx
#include <iostream>
#include "vector_operations.h"


//=================Vector operations tests=================
void TestVectorAddition()
{
  std::vector<double> vec1 = {1.0, 2.0, 3.0};
  std::vector<double> vec2 = {4.0, 5.0, 6.0};

  std::vector<double> result = vec1 + vec2; // Тривиальное сложение std::vector<T>

  std::cout << "Results for vector addition: ";
  for (double x : result) std::cout << x << " ";
  std::cout << "\n";
}

void TestVectorSubtraction()
{
  std::vector<double> vec1 = {10.1, 20.2, 30.3};
  std::vector<double> vec2 = {4.4, 5.5, 6.6};

  std::vector<double> result = vec1 - vec2; // Тривиальное вычитание std::vector<T>

  std::cout << "Results for vector subtraction: ";
  for (double x : result) std::cout << x << " ";
  std::cout << "\n";
}

void TestVectorScalarProduct()
{
  std::vector<double> vec1 = {1.0, 2.0, 3.0};
  std::vector<double> vec2 = {4.0, 5.0, 6.0};

  double result = vec1 * vec2; // Тривиальное скалярное произведение std::vector<T>

  std::cout << "Results for vector scalar product: " << result;
  std::cout << "\n";
}

void TestVectorScalarMultiplication()
{
  std::vector<double> vec = {1.0, 2.0, 3.0};
  double scalar = 2.0;

  std::vector<double> result = vec * scalar; // Тривиальное std::vector<T> на число

  std::cout << "Results for vector scalar multiplication: ";
  for (double x : result) std::cout << x << " ";
  std::cout << "\n";
}


int main()
{
  TestVectorAddition();
  TestVectorSubtraction();
  TestVectorScalarProduct();
  TestVectorScalarMultiplication();

  std::cout <<"\n";
  return 0;
}
```

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

    // К элементам матрицы можно обращаться через matrix(i, j)
    // Изменять можно только диагональные элементы, считывать - любые

    // Вектор свободных членов (должен иметь размер = matrix_size)
    std::vector<double> free_col = {3, 8, 5, 1};


    // Решаем систему
    Thomas_algorithm<double> solver(matrix, free_col);


    // Получаем решение
    std::vector<double> solution = solver.Get_solution();

    return 0;
}
```

### Пример работы с плотной (Dense) и разреженной (CSR) матрицами
```cxx
#include <iostream>
#include "CSR_matrix.h"
#include "Dense_matrix.h"

//=================Dense matrix tests=================
void DenseMatrixTest()
{
  std::vector<std::vector<double>> init_matrix =
  {
    {1.0, 2.0, 3.0},
    {4.0, 5.0, 6.0},
    {7.0, 8.0, 9.0}
  };

  DenseMatrix<double> dense_matrix(3, 3, init_matrix); // Создаем матрицу

  std::cout << "Element (0, 0): " << dense_matrix(0, 0) << " || should be 1.0\n"; // Выводим ее элементы с помощью ()
  std::cout << "Element (1, 2): " << dense_matrix(1, 2) << " || should be 6.0\n";
  std::cout << "Element (2, 1): " << dense_matrix(2, 1) << " || should be 8.0\n";
  std::cout << "Element (0, 2): " << dense_matrix(0, 2) << " || should be 3.0\n";

  std::vector<double> vec = {1.0, 1.0, 1.0}; // На этот вектор будем умножать
  std::vector<double> result = dense_matrix * vec;

  std::cout << "Dense matrix multiplication result: ";
  for (double x : result) std::cout << x << " ";
  std::cout << "\n";
}

//=================CSR matrix tests=================
void CSRMatrixTest()
{
  // Создаем DOK матрицу
  std::map<std::pair<size_t, size_t>, double> dok_matrix =
  {
    {{0, 0}, 1.0},
    {{1, 2}, 2.5},
    {{2, 1}, 3.7},
    {{3, 3}, 4.2}
  }; 

  CSR_Matrix<double> csr_matrix(dok_matrix, 5, 4);   // Создаем CSR матрицу, передавая в конструктор DOK матрицу и необходимые размеры

  std::cout << "Element (0, 0): " << csr_matrix(0, 0) << " || should be 1.0\n"; // Выводим ее элементы с помощью ()
  std::cout << "Element (1, 2): " << csr_matrix(1, 2) << " || should be 2.5\n"; // ВАЖНО - нельзя изменять элементы CSR матрицы
  std::cout << "Element (2, 1): " << csr_matrix(2, 1) << " || should be 3.7\n"; // Возможно, я потом это реализую
  std::cout << "Element (3, 3): " << csr_matrix(3, 3) << " || should be 4.2\n";
  std::cout << "Element (1, 1): " << csr_matrix(1, 1) << " || should be 0.0\n";
  std::cout << "Element (4, 4): " << csr_matrix(4, 3) << " || should be 0.0\n"; // Такой элемент мы не передавали, но передали размер - 
                                                                                // у нас одна нулевая строка снизу

  std::vector<double> vec = {1.0, 1.0, 1.0, 1.0}; // На этот вектор будем умножать
  std::vector<double> result = csr_matrix * vec;

  std::cout << "CSR matrix multiplication result: ";
  for (double x : result) std::cout << x << " ";
  std::cout << "\n";
}


int main()
{
  // Dense matrix tests
  DenseMatrixTest();

  std::cout <<"\n";

  // CSR matrix tests
  CSRMatrixTest();

  return 0;
}
```

## Что такое CSR матрица?
CSR матрица - это структура данных, которая удобна для хранения неплотных матриц. Она позволяет быстро перемножать неплотную матрицу на вектор. Вот сравнение скорости умножения в логарифмических координатах:
![Картинку съел Ульяновский автомеханический завод](https://github.com/GlebLarkin/Systems-of-linear-equations/blob/main/data/t(n).png)
Как видно из графиков, ассимтотики умножения обычной и CSR матрицы на вектор одинаковые, но CSR матрица умножает в среднем быстрее. Но если ее плоность равна 1.0, то скорость умножения такая же, как и для плотной матрицы, что в общем (и целом) ожидаемо. 

## Пара слов об итерационных методах
Ниже приведены графики сравнения скорости сходимости реализованных итерационных (именно зависимость невязки от количества итераций). На сходимость метода простых итераций сильно влияет выбор τ, оптимальный выбор которого является нетривиальной задачей. 
Эту задачу решает ускорение Чебышёва, которая, как можно видеть из графиков, действительно уменьшает количество итераций, необходимое для сходиомсти. К сожалению, для работы метода необходимо знать максимальное и минимальное собыственные значения матрицы, задающей СЛАУ. Максимальное собственне значение можно оценить с помощью функции Estimate_max_eigenvalue, а минимальное, к сожалению, оценить тяжело. 
![Изображение было поглощено консерном Алмаз-Антей](https://github.com/GlebLarkin/Systems-of-linear-equations/blob/main/data/all_iteration_methods.png)