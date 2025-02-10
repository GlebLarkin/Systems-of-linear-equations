#ifndef DENSEMATRIX_H
#define DENSEMATRIX_H

#include "IMatrix.h"


template <class T>
class DenseMatrix : public IMatrix<T>
{
public:
  DenseMatrix(size_t rows, size_t cols, const std::vector<std::vector<T>> & init_matrix)
    : nx_(cols), ny_(rows), data_(rows * cols)
  {
    if (init_matrix.size() != rows)
    {
      throw std::invalid_argument("Number of rows in init matrix does not match the number rows");
    }
    for (size_t i = 0; i < rows; ++i)
    {
      if (init_matrix[i].size() != cols)
      {
        throw std::invalid_argument("All rows must have the same number of columns");
      }
      for (size_t j = 0; j < cols; ++j)
      {
        data_[i * cols + j] = init_matrix[i][j];
      }
    }
  }

  DenseMatrix(const DenseMatrix<T> & other)
    : nx_(other.nx_), ny_(other.ny_), data_(other.data_) {}

  DenseMatrix(DenseMatrix<T> && other) noexcept
    : nx_(other.nx_), ny_(other.ny_), data_(std::move(other.data_))
  {
    other.nx_ = 0;
    other.ny_ = 0;
  }

  DenseMatrix<T> & operator = (const DenseMatrix<T> & other)
  {
    if (this != &other)
    {
      nx_ = other.nx_;
      ny_ = other.ny_;
      data_ = other.data_;
    }
    return *this;
  }

  DenseMatrix<T> & operator = (DenseMatrix<T> && other) noexcept
  {
    if (this != &other)
    {
      nx_ = other.nx_;
      ny_ = other.ny_;
      data_ = std::move(other.data_);

      other.nx_ = 0;
      other.ny_ = 0;
    }
    return *this;
  }

  T & operator() (size_t i, size_t j)
  {
    if (i >= ny_ || j >= nx_)
    {
      throw std::out_of_range("Trying to access an element outside the dense matrix\n");
    }
    return data_[i * nx_ + j];
  }

  T operator() (size_t i, size_t j) const override
  {
    if (i >= ny_ || j >= nx_)
    {
      throw std::out_of_range("Trying to access an element outside the dense matrix\n");
    }
    return data_[i * nx_ + j];
  }

  std::vector<T> operator * (const std::vector<T> & vec) const
  {
    size_t vec_size = vec.size();

    if (vec_size != nx_)
    {
      throw std::invalid_argument("size of vector and size of matrix do not match");
    }

    std::vector<T> ans;
    ans.reserve(ny_);
    for (size_t i = 0; i < ny_; ++i)
    {
      T ans_line = 0;
      for (size_t j = 0; j < nx_; ++j)
      {
        ans_line += (*this)(i, j) * vec[j];
      }
      ans.push_back(ans_line);
    }
    return ans;
  }
  
  std::pair<size_t, size_t> Get_matrix_size() const override { return std::make_pair(nx_, ny_); }

private:
  size_t nx_; // Number of columns
  size_t ny_; // Number of rows
  std::vector<T> data_;
};

#endif
