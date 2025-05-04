#ifndef UPPER_TRIANGULAR_MATRIX_H
#define UPPER_TRIANGULAR_MATRIX_H

#include <vector>
#include <stdexcept>

template<class T>
class UpperTriangularMatrix
{
public:
  UpperTriangularMatrix(const std::vector<T> & data, const size_t n)
    : data_(data), n_(n)
  {
    if (data.size() != n * (n + 1) / 2)
    {
      throw std::invalid_argument("Data size does not match upper triangular matrix size.");
    }
  }

  std::pair<size_t, size_t> shape() const { return {n_, n_}; }

  T operator() (size_t i, size_t j) const
  {
    if (i > j) return T{};
    return data_[index(i, j)];
  }

  void set(size_t i, size_t j, T val)
  {
    if (i > j)
    {
      throw std::invalid_argument("Shitty attempt to set value below main diagonal in upper triangular matrix.");
    }
    data_[index(i, j)] = val;
  }

  void addCol(const std::vector<T> & col)
  {=
    for (size_t i = 0; i <= n_; ++i)
    {
      data_.push_back(col[i]);
    }
    ++n_;
  }

private:
  size_t index(size_t i, size_t j) const { return i * n_ - i * (i - 1) / 2 + (j - i); }

  std::vector<T> data_;
  size_t n_;
};

#endif