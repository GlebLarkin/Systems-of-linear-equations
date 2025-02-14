#include <iostream>
#include "Thomas_algorithm.h"
#include "CSR_matrix.h"
#include "Dense_matrix.h"
#include "vector_operations.h"
#include <utility>

//=================Tridiagonal matrix tests=================
void TestSimpleCase()
{
  std::vector<double> lower_diag = {1, 1};
  std::vector<double> main_diag = {3, 3, 3};
  std::vector<double> upper_diag = {1, 1};
  std::vector<double> free_col = {1, 1, 1};

  std::vector<double> expected = {2.0 / 7, 1.0 / 7, 2.0 / 7};

  Tridiagonal_matrix<double> matrix(upper_diag, main_diag, lower_diag, 3);
  Thomas_algorithm<double> solver(matrix, free_col);

  std::vector<double> result = solver.Get_solution();

  std::cout << "Test 1 Thomas algorithm for tridiagonal matrix results: ";
  for (double x : result) std::cout << x << " ";
  std::cout << "\nExpected results: ";
  for (double x : expected) std::cout << x << " ";
  std::cout << "\n";
}

void TestAlsoSimpleCase()
{
  std::vector<float> lower_diag = {5, 1};
  std::vector<float> main_diag = {6, 8, 4};
  std::vector<float> upper_diag = {2, 2};
  std::vector<float> free_col = {1, 1, 2};

  std::vector<float> expected = {3.0 / 14, -1.0 / 7, 15.0 / 28};

  Tridiagonal_matrix<float> matrix(upper_diag, main_diag, lower_diag, 3);
  Thomas_algorithm<float> solver(matrix, free_col);

  std::vector<float> result = solver.Get_solution();

  std::cout << "Test 2 Thomas algorithm for tridiagonal matrix results: ";
  for (float x : result) std::cout << x << " ";
  std::cout << "\nExpected results: ";
  for (float x : expected) std::cout << x << " ";
  std::cout << "\n";
}

void FailTest()
{
  std::vector<double> lower_diag = {1, 1};
  std::vector<double> main_diag = {4, 0, 6};
  std::vector<double> upper_diag = {2, 2};
  std::vector<double> free_col = {7, 8, 9};

  try
  {
    std::cout << "Test 3 Thomas algorithm for tridiagonal matrix results (it should fail)\n";
    Tridiagonal_matrix<double> matrix(upper_diag, main_diag, lower_diag, 3);
    Thomas_algorithm<double> solver(matrix, free_col);
  }
  catch (const std::exception & e)
  {
    std::cout << e.what();
  }
}

//=================Vector operations tests=================
void TestVectorAddition()
{
  std::vector<double> vec1 = {1.0, 2.0, 3.0};
  std::vector<double> vec2 = {4.0, 5.0, 6.0};
  std::vector<double> expected = {5.0, 7.0, 9.0};

  std::vector<double> result = vec1 + vec2;

  std::cout << "Test 1 for vector addition: ";
  for (double x : result) std::cout << x << " ";

  std::cout << "\nExpected results: ";
  for (double x : expected) std::cout << x << " ";
  std::cout << "\n";
}

void TestVectorSubtraction1()
{
  std::vector<double> vec1 = {10.1, 20.2, 30.3};
  std::vector<double> vec2 = {4.4, 5.5, 6.6};
  std::vector<double> expected = {5.7, 14.7, 23.7};

  std::vector<double> result = vec1 - vec2;

  std::cout << "Test 1 for vector subtraction: ";
  for (double x : result) std::cout << x << " ";

  std::cout << "\nExpected results: ";
  for (double x : expected) std::cout << x << " ";
  std::cout << "\n";
}

void TestVectorSubtraction2()
{
  std::vector<double> vec1 = {10.1, 20.2, 30.3};
  std::vector<double> vec2 = {4.4, 5.5};

  std::cout << "Test 2 for vector subtraction: ";

  try
  {
    std::vector<double> result = vec1 - vec2;
    for (double x : result) std::cout << x << " ";
  }
  catch (const std::runtime_error & e)
  {
    std::cout << e.what();
  }
  std::cout << "\nExpected results: should fail\n";
}

void TestVectorScalarProduct()
{
  std::vector<double> vec1 = {1.0, 2.0, 3.0};
  std::vector<double> vec2 = {4.0, 5.0, 6.0};
  double expected = 32.0;

  double result = vec1 * vec2;

  std::cout << "Test 1 for vector scalar product: " << result;
  std::cout << "\nExpected result: " << expected << "\n";
}

void TestVectorScalarMultiplication()
{
  std::vector<double> vec = {1.0, 2.0, 3.0};
  double scalar = 2.0;
  std::vector<double> expected = {2.0, 4.0, 6.0};

  std::vector<double> result = vec * scalar;

  std::cout << "Test 1 for vector scalar multiplication: ";
  for (double x : result) std::cout << x << " ";

  std::cout << "\nExpected results: ";
  for (double x : expected) std::cout << x << " ";
  std::cout << "\n";
}

//=================Dense matrix tests=================
void DenseMatrixTest()
{
  std::vector<std::vector<double>> init_matrix =
  {
    {1.0, 2.0, 3.0},
    {4.0, 5.0, 6.0},
    {7.0, 8.0, 9.0}
  };

  DenseMatrix<double> dense_matrix(3, 3, init_matrix);

  std::cout << "=======Dense matrix tests=======\n";

  std::cout << "Element (0, 0): " << dense_matrix(0, 0) << " || should be 1.0\n";
  std::cout << "Element (1, 2): " << dense_matrix(1, 2) << " || should be 6.0\n";
  std::cout << "Element (2, 1): " << dense_matrix(2, 1) << " || should be 8.0\n";
  std::cout << "Element (0, 2): " << dense_matrix(0, 2) << " || should be 3.0\n";

  std::vector<double> vec = {1.0, 1.0, 1.0};
  std::vector<double> expected_result = {6.0, 15.0, 24.0};
  std::vector<double> result = dense_matrix * vec;

  std::cout << "Dense matrix multiplication result: ";
  for (double x : result) std::cout << x << " ";
  std::cout << "\nExpected result: ";
  for (double x : expected_result) std::cout << x << " ";
  std::cout << "\n";
}

//=================CSR matrix tests=================
void CSRMatrixTest()
{
  std::map<std::pair<size_t, size_t>, double> dok_matrix =
  {
    {{0, 0}, 1.0},
    {{1, 2}, 2.5},
    {{2, 1}, 3.7},
    {{3, 3}, 4.2}
  };

  CSR_Matrix<double> csr_matrix(dok_matrix, 5, 4);

  std::cout << "=======CSR matrix tests=======\n";

  std::cout << "Element (0, 0): " << csr_matrix(0, 0) << " || should be 1.0\n";
  std::cout << "Element (1, 2): " << csr_matrix(1, 2) << " || should be 2.5\n";
  std::cout << "Element (2, 1): " << csr_matrix(2, 1) << " || should be 3.7\n";
  std::cout << "Element (3, 3): " << csr_matrix(3, 3) << " || should be 4.2\n";
  std::cout << "Element (1, 1): " << csr_matrix(1, 1) << " || should be 0.0\n";
  std::cout << "Element (4, 4): " << csr_matrix(4, 3) << " || should be 0.0\n";

  std::vector<double> vec = {1.0, 1.0, 1.0, 1.0};
  std::vector<double> expected_result = {1.0, 2.5, 3.7, 4.2};
  std::vector<double> result = csr_matrix * vec;

  std::cout << "CSR matrix multiplication result: ";
  for (double x : result) std::cout << x << " ";
  std::cout << "\nExpected result: ";
  for (double x : expected_result) std::cout << x << " ";
  std::cout << "\n";
}





int main()
{
  // Tridiagonal matrix tests
  TestSimpleCase();
  TestAlsoSimpleCase();
  FailTest();

  std::cout <<"\n";

  // Vector operations tests
  TestVectorAddition();
  TestVectorSubtraction1();
  TestVectorSubtraction2();
  TestVectorScalarProduct();
  TestVectorScalarMultiplication();

  std::cout <<"\n";

  // Dense matrix tests
  DenseMatrixTest();

  std::cout <<"\n";

  // CSR matrix tests
  CSRMatrixTest();

  return 0;
}