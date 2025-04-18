#ifndef MATRIX_OP
#define MATRIX_OP

#include <limits>
#include <cmath>
#include <iostream>
#include "./Matrixes/Dense_matrix.h"
#include "./Matrixes/CSR_matrix.h"
#include "./Operations/vector_operations.h"


template <class T>
DenseMatrix<T> Vector2matrix(const std::vector<T> & vec)
{
  std::vector<std::vector<T>> matrix(vec.size(), std::vector<T>(1)); 
  for (size_t i = 0; i < vec.size(); ++i)
  {
    matrix[i][0] = vec[i];
  }
  return DenseMatrix<T>(matrix, vec.size(), 1);
}

template <class T>
std::vector<T> Matrix2vectorRow(const IMatrix<T> &  matrix, const size_t i)
{
  auto [rows, cols] = matrix.Get_matrix_size();
  if (i >= rows)
  {
    throw std::out_of_range("Index of row is out of range");
  }
  std::vector<T> vec(cols);
  for (size_t j = 0; j < cols; ++j)
  {
    vec[j] = matrix(i, j);
  }
  return vec;
}

template <class T>
std::vector<T> Matrix2vectorCol(const IMatrix<T> & matrix, const size_t j)
{
  auto [rows, cols] = matrix.Get_matrix_size();
  if (j >= cols)
  {
    throw std::out_of_range("Index of column is out of range");
  }
  std::vector<T> vec(rows);
  for (size_t i = 0; i < rows; ++i)
  {
    vec[i] = matrix(i, j);
  }
  return vec;
}

template <class T>
DenseMatrix<T> operator * (const DenseMatrix<T> & A, const DenseMatrix<T> & B)
{
  auto [A_rows, A_cols] = A.Get_matrix_size();
  auto [B_rows, B_cols] = B.Get_matrix_size();

  if (A_cols != B_rows)
  {
    throw std::invalid_argument("Matrix dimensions do not match for multiplication.");
  }

  std::vector<std::vector<T>> result(A_rows, std::vector<T>(B_cols, T(0)));

  for (size_t i = 0; i < A_rows; ++i)
  {
    for (size_t j = 0; j < B_cols; ++j)
    {
      for (size_t k = 0; k < A_cols; ++k)
      {
        result[i][j] += A(i, k) * B(k, j);
      }
    }
  }

  return DenseMatrix<T>(result, A_rows, B_cols);
}


template <class T>
DenseMatrix<T> operator + (const DenseMatrix<T> & lhs, const DenseMatrix<T> & rhs)
{
  auto [rows1, cols1] = lhs.Get_matrix_size();
  auto [rows2, cols2] = rhs.Get_matrix_size();

  if (rows1 != rows2 || cols1 != cols2)
  {
    throw std::invalid_argument("Matrix sizes are different");
  }

  std::vector<std::vector<T>> result(rows1, std::vector<T>(cols1));

  for (size_t i = 0; i < rows1; ++i)
  {
    for (size_t j = 0; j < cols1; ++j)
    {
      result[i][j] = lhs(i, j) + rhs(i, j);
    }
  }

  return DenseMatrix<T>(result, rows1, cols1);
}

template <class T>
DenseMatrix<T> operator - (const DenseMatrix<T> & lhs, const DenseMatrix<T> & rhs)
{
  auto [rows1, cols1] = lhs.Get_matrix_size();
  auto [rows2, cols2] = rhs.Get_matrix_size();

  if (rows1 != rows2 || cols1 != cols2)
  {
    throw std::invalid_argument("Matrix sizes are different");
  }

  std::vector<std::vector<T>> result(rows1, std::vector<T>(cols1));

  for (size_t i = 0; i < rows1; ++i)
  {
    for (size_t j = 0; j < cols1; ++j)
    {
      result[i][j] = lhs(i, j) - rhs(i, j);
    }
  }

  return DenseMatrix<T>(result, rows1, cols1);
}

template <class T>
DenseMatrix<T> operator + (const Eye<T> & lhs, const Eye<T> & rhs)
{
  auto [size1, _] = lhs.Get_matrix_size();
  auto [size2, __] = rhs.Get_matrix_size();

  if (size1 != size2)
  {
    throw std::invalid_argument("Matrix sizes are different");
  }

  return DenseMatrix<T>(std::vector<std::vector<T>>(size1, std::vector<T>(size1, 2)), size1, size1);
}

template <class T>
DenseMatrix<T> operator - (const Eye<T> & lhs, const Eye<T> & rhs)
{
  auto [size1, _] = lhs.Get_matrix_size();
  auto [size2, __] = rhs.Get_matrix_size();

  if (size1 != size2)
  {
    throw std::invalid_argument("Matrix sizes are different");
  }

  return DenseMatrix<T>(std::vector<std::vector<T>>(size1, std::vector<T>(size1, 0)), size1, size1);
}

template <class T>
DenseMatrix<T> operator + (const DenseMatrix<T> & lhs, const Eye<T> & rhs)
{
  auto [rows, cols] = lhs.Get_matrix_size();
  auto [size, _] = rhs.Get_matrix_size();

  if (rows != size || cols != size)
  {
    throw std::invalid_argument("Matrix sizes are different");
  }

  std::vector<std::vector<T>> result(rows, std::vector<T>(cols));

  for (size_t i = 0; i < rows; ++i)
  {
    for (size_t j = 0; j < cols; ++j)
    {
      result[i][j] = lhs(i, j) + rhs(i, j);
    }
  }

  return DenseMatrix<T>(result, rows, cols);
}

template <class T>
DenseMatrix<T> operator + (const Eye<T> & lhs, const DenseMatrix<T> & rhs)
{
  return rhs + lhs;
}

template <class T>
DenseMatrix<T> operator - (const DenseMatrix<T> & lhs, const Eye<T> & rhs)
{
  auto [rows, cols] = lhs.Get_matrix_size();
  auto [size, _] = rhs.Get_matrix_size();

  if (rows != size || cols != size)
  {
    throw std::invalid_argument("Matrix sizes are different");
  }

  std::vector<std::vector<T>> result(rows, std::vector<T>(cols));

  for (size_t i = 0; i < rows; ++i)
  {
    for (size_t j = 0; j < cols; ++j)
    {
      result[i][j] = lhs(i, j) - rhs(i, j);
    }
  }

  return DenseMatrix<T>(result, rows, cols);
}

template <class T>
DenseMatrix<T> operator - (const Eye<T> & lhs, const DenseMatrix<T> & rhs)
{
  auto [rows, cols] = rhs.Get_matrix_size();
  auto [size, _] = lhs.Get_matrix_size();

  if (rows != size || cols != size)
  {
    throw std::invalid_argument("Matrix sizes are different");
  }

  std::vector<std::vector<T>> result(rows, std::vector<T>(cols));

  for (size_t i = 0; i < rows; ++i)
  {
    for (size_t j = 0; j < cols; ++j)
    {
      result[i][j] = lhs(i, j) - rhs(i, j);
    }
  }

  return DenseMatrix<T>(result, rows, cols);
}


template <class T>
T Estimate_max_eigenvalue(CSR_Matrix<T> & matrix)
{
  std::vector<T> x = std::vector<T>(matrix.Get_matrix_size().second, 1);
  x = x * (1 / VectorNorm(x));

  const size_t max_iteration_number = 10000;
  const T tol = std::numeric_limits<T>::epsilon() * 1e3;
  T lambda = 0;
  T lambda_new = 0;

  for (size_t iteration = 0; iteration < max_iteration_number; ++iteration)
  {
    auto y = matrix * x;
    y = y * (1 / VectorNorm(y));

    auto y_matrix = Vector2matrix(y);
    lambda_new = (y_matrix.Transpond() * Vector2matrix(matrix * y))(0, 0) / (y_matrix.Transpond() * y_matrix)(0, 0);

    if (std::abs(lambda_new - lambda) < tol)
      return lambda_new;

    lambda = lambda_new;
    x = y;
  }

  std::cout << "Iteration number limit had been reached, but the discrepancy is too big. Try to compute eigenvalues yourself.\n";
  return lambda_new;
}


#endif