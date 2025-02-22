#ifndef QR_H
#define QR_H

#include "Dense_matrix.h"
#include "Some_matrix_operations.h"
#include "vector_operations.h"

template <class T>
class QR
{
public:
  QR(const DenseMatrix<T> & matrix) 
    : matrix_(matrix), rows_(matrix.Get_matrix_size().first), cols_(matrix.Get_matrix_size().second) 
  {
    if (rows_ < cols_)
    {
      throw std::invalid_argument("QR decomposition is only for tall matrixes (rows >= cols)");
    }
    Decompose(); 
  }

  DenseMatrix<T> GetQ() const { return Q_matrix; }
  DenseMatrix<T> GetR() const { return R_matrix; }

private:
  std::vector<T> Reflection_result(const std::vector<T> & x, const std::vector<T> & e)
  {
    if (VectorNorm(x) == T(0)) throw std::invalid_argument("Can't reflect zero vector");
    std::vector<T> v = x - e * VectorNorm(x);
    return x - v * (2. * (v * x) / (v * v));
  }

  void Decompose() // Householder algorithm
  {
    std::vector<std::vector<T>> Q(rows_, std::vector<T>(rows_, 0));  // Q - rows x rows
    std::vector<std::vector<T>> R(rows_, std::vector<T>(cols_, 0));  // R - (rows x cols)

    // R matrix creation
    size_t row = rows_;

    for (size_t col = 0; col < cols_; ++col)
    {
      std::vector<T> e(row, 0);
      e[0] = 1;

      try
      {
        std::vector<T> matrix_col = Cut(Matrix2vectorCol(matrix_, col), rows_ - row, rows_ - 1);

        R[col] = Reflection_result(matrix_col, e);
      }
      catch(const std::invalid_argument &) // if norm of vector x is zero we do not need to reflect it
      {
        --row;
        continue;
      }

      --row;
    }

    // Q matrix creation

    Eye<T> eye_matrix(rows_);
    row = rows_;

    for (size_t row_ind = 0; row_ind < rows_; ++row_ind)
    {
      std::vector<T> e(row, 0);
      e[0] = 1;

      try
      {
        std::vector<T> eye_row = Cut(Matrix2vectorRow(eye_matrix, row_ind), rows_ - row, rows_ - 1);

        Q[row_ind] = Reflection_result(eye_row, e);
      }
      catch(const std::invalid_argument &) // if norm of vector x is zero we do not need to reflect it
      {
        --row;
        continue;
      }

      --row;
    }

    Q_matrix = DenseMatrix<T>(Q, rows_, rows_).Transpond();
    R_matrix = DenseMatrix<T>(R, rows_, cols_);
  }

  DenseMatrix<T> matrix_;
  DenseMatrix<T> Q_matrix;
  DenseMatrix<T> R_matrix;

  size_t rows_;
  size_t cols_;
};

#endif
