#ifndef THOMASALGORITHM_H
#define THOMASALGORITHM_H

#include "./Matrixes/Tridiagonal_matrix.h"


template <class T>
class Thomas_algorithm
{
public:
  Thomas_algorithm(const Tridiagonal_matrix<T> & matrix, const std::vector<T> & free_col)
    : matrix_(matrix),
      upper_diag_(matrix.Get_upper_diag()),
      main_diag_(matrix.Get_main_diag()),
      lower_diag_(matrix.Get_lower_diag()),
      free_col_(free_col) 
      {
        Check_if_matrix_correct();
        Solve();
      }

  std::vector<T> Get_solution() const { return solution_; }

private:
  void Check_if_matrix_correct() const
  {
    if (matrix_.Get_matrix_size().first != free_col_.size())
      {
        throw std::runtime_error("Tridiagonal matrix and free colomb sizes are different\n");
      }

    for (size_t i = 0; i < matrix_.Get_matrix_size().first; ++i)
    {
      T sum = std::abs(upper_diag_[i]) + std::abs(lower_diag_[i]);
      if (std::abs(main_diag_[i]) <= sum)
      {
        std::cout << "!!!\nThe condition of strictly diagonal predominance is not fulfilled, the algorithm is unstable\n!!!\n";
        //std::cerr << "!!!\nThe condition of strictly diagonal predominance is not fulfilled, the algorithm is unstable\n!!!\n";
      }
    }
  }

  void Solve()
  {
    size_t size = matrix_.Get_matrix_size().first;
    std::vector<T> p_arr(size);
    std::vector<T> q_arr(size);

    p_arr[0] = - (upper_diag_[0]) / (main_diag_[0]);
    q_arr[0] = (free_col_[0]) / (main_diag_[0]);
    
    for (size_t i = 1; i < size; ++i)
    {
      p_arr[i] = - (upper_diag_[i]) / (main_diag_[i] + lower_diag_[i] * p_arr[i - 1]);
      q_arr[i] = (free_col_[i] - lower_diag_[i] * q_arr[i - 1]) / (main_diag_[i] + lower_diag_[i] * p_arr[i - 1]);
    }

    std::vector<T> solution(size);

    solution[size - 1] = q_arr[size - 1];

    for (int i = size - 2; i >= 0; --i)
    {
      solution[i] = p_arr[i] * solution[i + 1] + q_arr[i];
    }

    solution_ = solution;
  }
  
  Tridiagonal_matrix<T> matrix_;
  std::vector<T> upper_diag_;
  std::vector<T> main_diag_;
  std::vector<T> lower_diag_;
  std::vector<T> free_col_;

  std::vector<T> solution_;
};


#endif