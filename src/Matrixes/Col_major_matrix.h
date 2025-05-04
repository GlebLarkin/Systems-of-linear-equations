/*
#ifndef COL_MAJOR_MATRIX_H
#define COL_MAJOR_MATRIX_H

#include <vector>
#include <span>


template<class T>
class Col_Major_Matrix 
{
public:
  Col_Major_Matrix(const std::vector<T> & data, const std::size_t N) : data_(data), N_(N) {}

  Col_Major_Matrix(const std::size_t N) 
  {
    data_.resize(N * N);
    N_ = N;
  }

  std::span<T> col(size_t i) { return std::span<T>(data_.begin() + N_ * i, N_); }

  std::span<const T> col(size_t i) const { return std::span<const T>(data_.begin() + N_ * i, N_); }

private:
  std::size_t N_;
  std::vector<T> data_;
};

#endif
*/