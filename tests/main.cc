#include <iostream>
#include "Thomas_algorithm.h"


//я хз кто такие тесты

void TestSimpleCase()
{
  std::vector<double> lower_diag = {1};
  std::vector<double> main_diag = {4, 5}; 
  std::vector<double> upper_diag = {2}; 
  std::vector<double> free_col = {6, 9}; 

  Tridiagonal_matrix<double> matrix(upper_diag, main_diag, lower_diag, 2);
  Thomas_algorithm<double> solver(matrix, free_col);

  std::vector<double> result = solver.Get_solution();
  std::cout << "Test 1: Simple case result: ";
  for (double x : result) std::cout << x << " ";
  std::cout << std::endl;
}

void TestAlsoSimpleCase() 
{
  std::vector<double> lower_diag = {1, 1};
  std::vector<double> main_diag = {4, 5, 6};
  std::vector<double> upper_diag = {2, 2};
  std::vector<double> free_col = {7, 8, 9};

  Tridiagonal_matrix<double> matrix(upper_diag, main_diag, lower_diag, 3);
  Thomas_algorithm<double> solver(matrix, free_col);

  std::vector<double> result = solver.Get_solution();
  std::cout << "Test 2: Typical case result: ";
  for (double x : result) std::cout << x << " ";
  std::cout << std::endl;
}

void FailTest() 
{
  std::vector<double> lower_diag = {1, 1};
  std::vector<double> main_diag = {4, 0, 6};
  std::vector<double> upper_diag = {2, 2};
  std::vector<double> free_col = {7, 8, 9};

  try 
  {
    Tridiagonal_matrix<double> matrix(upper_diag, main_diag, lower_diag, 3);
    Thomas_algorithm<double> solver(matrix, free_col);
    std::cout << "Test 3: Unexpected success (should fail!)\n";
  } catch (const std::exception &e) {
    std::cout << "Test 3: Caught expected error: " << e.what() << std::endl;
  }
}

void TestLargeSystem() 
{
  size_t n = 100;
  std::vector<double> lower_diag(n - 1, -1);
  std::vector<double> main_diag(n, 4);
  std::vector<double> upper_diag(n - 1, -1);
  std::vector<double> free_col(n, 2);

  Tridiagonal_matrix<double> matrix(upper_diag, main_diag, lower_diag, n);
  Thomas_algorithm<double> solver(matrix, free_col);

  std::vector<double> result = solver.Get_solution();
  std::cout << "Test 4: Large system, first 5 results: ";
  for (size_t i = 0; i < n; i++) std::cout << result[i] << " ";
  std::cout << std::endl;
}

int main() 
{
    TestSimpleCase();
    TestAlsoSimpleCase();
    FailTest();
    TestLargeSystem();
    return 0;
}