#ifndef QR_H
#define QR_H

#include "./Matrixes/Dense_matrix.h"
#include "./Operations/Some_matrix_operations.h"
#include "./Operations/vector_operations.h"
#include <iostream>
#include <stdexcept>

template <class T>
class QR
{
public:
  QR(const DenseMatrix<T> & matrix)
    : rows_(matrix.Get_matrix_size().first), cols_(matrix.Get_matrix_size().second), R_matrix(matrix)
  {
    if (rows_ < cols_)
    {
      throw std::invalid_argument("QR decomposition is only for tall matrices (rows >= columns)");
    }

    std::vector<std::vector<T>> Q_init(rows_, std::vector<T>(rows_, 0));
    for (size_t i = 0; i < rows_; ++i) { Q_init[i][i] = 1; }
    Q_matrix = DenseMatrix<T>(Q_init, rows_, rows_);

    Decompose();
  }

  std::vector<T> solve_system(const std::vector<T> & b) const
  {
    if (b.size() != rows_)
    {
      throw std::invalid_argument("Vector b must have the same number of rows as the matrix.");
    }
  
    std::vector<T> QTb = Q_matrix.Transpond() * b;
  
    std::vector<T> x(cols_, 0);
    for (int i = cols_ - 1; i >= 0; --i)
    {
      T sum = QTb[i];
      for (size_t j = i + 1; j < cols_; ++j)
      {
        sum -= R_matrix(i, j) * x[j];
      }
      x[i] = sum / R_matrix(i, i);
    }
  
    return x;
  }
  
  

  DenseMatrix<T> GetQ() const { return Q_matrix; }
  DenseMatrix<T> GetR() const { return R_matrix; }

private:
  void Decompose()
  {
    for (size_t col = 0; col < cols_; ++col)
    {
      std::vector<T> e(rows_ - col, 0);
      e[0] = 1;

      auto vector_R = Matrix2vectorCol(R_matrix, col);
      auto xR = Cut(vector_R, col, rows_ - 1);
      auto v = xR - e * VectorNorm(xR);
      auto norm_v = VectorNorm(v);
      if (norm_v == T(0)) continue;

      for (size_t right_matrix_col = col; right_matrix_col < cols_; ++right_matrix_col)
      {
        auto vector_for_change = Cut(Matrix2vectorCol(R_matrix, right_matrix_col), col, rows_ - 1);
        auto R_matrix_col = vector_for_change - v * (2.0 * (v * vector_for_change) / (norm_v * norm_v));

        for (size_t i = col; i < rows_; ++i)
        {
          R_matrix(i, right_matrix_col) = R_matrix_col[i - col];
        }
      }

      for (size_t right_matrix_col = 0; right_matrix_col < rows_; ++right_matrix_col)
      {
        auto vector_for_change = Cut(Matrix2vectorCol(Q_matrix, right_matrix_col), col, rows_ - 1);
        auto Q_matrix_col = vector_for_change - v * (2.0 * (v * vector_for_change) / (norm_v * norm_v));

        for (size_t i = col; i < rows_; ++i)
        {
          Q_matrix(i, right_matrix_col) = Q_matrix_col[i - col];
        }
      }
    }
    Q_matrix.TranspondInPlace();
  }

  DenseMatrix<T> Q_matrix;
  DenseMatrix<T> R_matrix;

  size_t rows_;
  size_t cols_;
};

#endif
