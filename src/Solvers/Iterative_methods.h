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
#include <cmath>
#include <numbers>


template <class T>
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
    if (iteration >= max_iterations_) std::cout << "Iteration number limit had been reached, but the discrepancy is too big. Maybe the method is unstable.";
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
    if (iteration >= max_iterations_) std::cout << "Iteration number limit had been reached, but the discrepancy is too big. Maybe the method is unstable.";
    return x_;
  }
  
  std::vector<T> Simple_iteration_method(T tau = 0.001)
  {
    size_t iteration = 0;

    while (!check_discrepancy() && iteration < max_iterations_) 
    {
      x_ = x_ - (A_ * x_ - b_) * tau;
      ++iteration;
    }
    if (iteration >= max_iterations_) std::cout << "Iteration number limit had been reached, but the discrepancy is too big. Maybe the method is unstable.\n";

    return x_;
  }

  std::vector<T> Chebyshev_simple_iteration_method(T lambda_min, T lambda_max, size_t number_of_iterations = 32) 
  // can be sped up by pre-computing the tau vector without recursion
  {
    if (number_of_iterations == 0) throw std::invalid_argument("Holy moly zero number of iterations is no good");
    auto findExponent = [](T x) -> size_t 
    {
      size_t n = 0;
      T power = 1;
      while (power < x) 
      {
        ++n;
        power *= 2;
      }
      return n;
    };
    number_of_iterations = static_cast<size_t>(pow(2, findExponent(number_of_iterations)));
    
    if (number_of_iterations > max_iterations_) number_of_iterations = max_iterations_;
    auto true_tau_arr = create_true_tau_arr(lambda_min, lambda_max, number_of_iterations);

    size_t iteration = 0;
    while(!check_discrepancy() && iteration < max_iterations_ && iteration < number_of_iterations)
    {
      x_ = x_ - (A_ * x_ - b_) * true_tau_arr[iteration];
      ++iteration;
    }

    if (iteration >= max_iterations_) std::cout << "Iteration number limit had been reached, but the discrepancy is too big. Maybe the method is unstable.";
    if (iteration >= number_of_iterations) return Chebyshev_simple_iteration_method(lambda_min, lambda_max, number_of_iterations * 2);
    return x_;
  }

  void ResetSolution() { x_.assign(x_.size(), 0); }

private:
  bool check_discrepancy() { return VectorNorm(A_ * x_ - b_) < stop_discrepancy_; } const

  // creats the correct order of tau for Chebyshev_simple_iteration_method
  inline std::vector<size_t> get_tau_order(size_t number_of_iterations, std::vector<size_t> init_tau_arr = std::vector<size_t>{0}) const // я могу как то передавать по ссылке init_tau_arr?
  {
    size_t current_tau_size = init_tau_arr.size();

    if (current_tau_size == number_of_iterations) return init_tau_arr;

    std::vector<size_t> tau_arr;
    tau_arr.reserve(2 * current_tau_size);

    for (size_t i = 0; i < current_tau_size; ++i)
    {
      tau_arr.emplace_back(init_tau_arr[i]);
      tau_arr.emplace_back(2 * current_tau_size - 1 - init_tau_arr[i]);
    }

    return get_tau_order(number_of_iterations, tau_arr);
  }

  inline std::vector<T> find_Chebyshev_coefs(size_t number_of_iterations) const
  {
    T cos_1 = std::cos(std::numbers::pi_v<T> / number_of_iterations);
    T sin_1 = std::sin(std::numbers::pi_v<T> / number_of_iterations);
    T t_0 = std::cos(std::numbers::pi_v<T> / (2 * number_of_iterations));

    std::vector<T> ans;
    ans.reserve(number_of_iterations);

    ans.emplace_back(t_0);
    for (size_t i = 0; i < number_of_iterations - 1; ++i) // previous indexes of ans
    {
      T value = std::sqrt(std::max(T(0), 1 - ans[i] * ans[i])) * sin_1;
      ans.emplace_back(ans[i] * cos_1 - value);
    }

    return ans;
  }
  
  inline std::vector<T> create_true_tau_arr(T lambda_min, T lambda_max, size_t number_of_iterations) const
  {
    auto init_tau_arr = find_Chebyshev_coefs(number_of_iterations);
    for (size_t i = 0; i < number_of_iterations; ++i)
    {
      init_tau_arr[i] = 0.5 * (lambda_min + lambda_max) + 0.5 * (lambda_max - lambda_min) * init_tau_arr[i];
    }
  
    std::vector<size_t> tau_order = get_tau_order(number_of_iterations);
  
    std::vector<T> true_tau_arr;
    true_tau_arr.reserve(number_of_iterations);
    
    for (size_t it : tau_order)
    {
      true_tau_arr.emplace_back(1.0 / init_tau_arr[it]);
    }
    
    return true_tau_arr;
  }
  


  const CSR_Matrix<T> A_;
  const std::vector<T> b_;
  const T stop_discrepancy_;  // when ( |Ax - b| < stop_discrepancy_ ) we stop computation
  size_t max_iterations_ = 1000000;

  std::vector<T> x_;  // answer

};

#endif