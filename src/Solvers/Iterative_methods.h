#ifndef ITERATIVE_METHODS_H
#define ITERATIVE_METHODS_H

#include "./Matrixes/Dense_matrix.h"
#include "./Matrixes/CSR_matrix.h"
#include "./Matrixes/Tridiagonal_matrix.h"
#include "./Operations/Some_matrix_operations.h"
#include "./Operations/vector_operations.h"

#include <limits>
#include <stdexcept>
#include <iostream>

template<class T>
class IterativeMethodsSolver
{
public:
  IterativeMethodsSolver(const CSR_Matrix<T> & A,  
                         const std::vector<T> & b, 
                         const T stop_discrepancy = std::numeric_limits<T>::epsilon() * 1e3) 
    : A_(A), b_(b), stop_discrepancy_(stop_discrepancy), x_(std::vector<T>(b.size(), 0))
  {
    const auto [rows, cols] = A_.Get_matrix_size();
    if (rows != cols)
    {
      throw std::invalid_argument("Holy shit bro. You try to solve system of linear equations with not square matrix");
    }
    
    if (rows != b_.size())
    {
      throw std::invalid_argument("Size of vector and size of matrix do not match for iterative methods solver");
    }
  }
  
  std::vector<T> Jacoby_method()
  {
    const auto [rows, cols] = A_.Get_matrix_size();

    DenseMatrix<T> D_min1(std::vector<std::vector<T>>(rows, std::vector<T>(cols, 0)), rows, cols);
    DenseMatrix<T> D_diag(std::vector<std::vector<T>>(rows, std::vector<T>(cols, 0)), rows, cols);

    for (size_t i = 0; i < rows; ++i)
    {
      T diag = std::abs(A_(i, i));
      T sum_off_diag = 0;

      for (size_t j = 0; j < cols; ++j) { if (i != j) { sum_off_diag += std::abs(A_(i, j)); } }

      if (diag <= sum_off_diag) { std::cout << "Ahtung! Matrix does not have diagonal dominance, Jacobi method is not stable.\n"; break; }
    }

    for (size_t i = 0; i < rows; ++i)
    {
      D_min1(i, i) = 1 / A_(i, i);
      D_diag(i, i) = A_(i, i);
    }

    size_t iteration = 0;

    while (!check_discrepancy() && iteration < max_iterations_) 
    {
      x_ = D_min1 * (b_ - (A_ * x_ - D_diag * x_)); 
      ++iteration;
    }

    return x_;
  }

  std::vector<T> Gauss_Seidel_method()
  {
    const auto [rows, cols] = A_.Get_matrix_size();
    for (size_t i = 0; i < rows; ++i)
    {
      if (A_(i, i) == 0) throw std::invalid_argument("Zero element in init matrix. Gauss Seidel method can not be used.");
    }
  
    size_t iteration = 0;
  
    while (!check_discrepancy() && iteration < max_iterations_)
    {
      for (size_t k = 0; k < rows; ++k)
      {
        T sum1 = 0, sum2 = 0;

        for (size_t j = k + 1; j < cols; ++j) { sum1 += A_(k, j) * x_[j]; }
        for (size_t j = 0; j < k; ++j) { sum2 += A_(k, j) * x_[j]; }
  
        x_[k] = (b_[k] - sum1 - sum2) / A_(k, k);
      }
      ++iteration;
    }
  
    return x_;
  }
  
  std::vector<T> Simple_iteration_method(T tau = 1)
  {
    size_t iteration = 0;

    while (!check_discrepancy() && iteration < max_iterations_) 
    {
      x_ = x_ - (A_ * x_ - b_) * tau;
      ++iteration;
    }

    return x_;
  }

  void ResetSolution() { x_.assign(b_.size(), 0); }
  
private:
  bool check_discrepancy() { return VectorNorm(A_ * x_ - b_) < stop_discrepancy_; }

  const CSR_Matrix<T> A_;
  const std::vector<T> b_;
  const T stop_discrepancy_;  // when ( |Ax - b| < stop_discrepancy_ ) we stop computation
  size_t max_iterations_ = 10000;

  std::vector<T> x_;  // answer

};

#endif