#include <iostream>
#include "Thomas_algorithm.h"


//я хз кто такие тесты

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
  } catch (const std::exception & e) {}
}


int main() 
{
    TestSimpleCase();
    TestAlsoSimpleCase();
    FailTest();
    return 0;
}