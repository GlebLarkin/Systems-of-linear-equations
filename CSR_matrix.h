#ifndef CSR_H
#define CSR_H

#include <algorithm>
#include <map>
#include "IMatrix.h"


template <class T>
class CSR_Matrix : public IMatrix<T>
{
public:
  // using DOK initialization
  CSR_Matrix(const std::map<std::pair<size_t, size_t>, T> & init_matrix) 
    : min_sizes(Get_DOK_sizes(init_matrix)) 
    {
      // Reserve memory
      size_t init_matrix_size = init_matrix.size();
      values.reserve(init_matrix_size);
      cols.reserve(init_matrix_size);
      rows.resize(min_sizes.first + 1, 0);

      // start to init
      size_t current_row = 0;

      for (const auto & it : init_matrix)
      {
        for (size_t r = current_row + 1; r <= it.first.first; ++r) 
        {
          rows[r] = values.size();
        }
        current_row = it.first.first;

        values.push_back(it.second);
        cols.push_back(it.first.second);
      }

      rows[min_sizes.first] = values.size();
    }

  T operator() (size_t i, size_t j) const override
  {
    if (i >= min_sizes.first || j >= min_sizes.second) return T(0);

    size_t value_row_begin = rows[i];
    size_t value_row_end = rows[i + 1];

    for (size_t k = value_row_begin; k < value_row_end; ++k) 
    {
      if (cols[k] == j) return values[k];
    }

    return T(0);
  }

  //T& operator() (size_t i, size_t j) = delete; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! спросить семера
  // if it is not deleted, user will be able to change 
  // element with any large indexes 
  // because we can have matrixes with very 
  // large number of zero rows and cols in the end

  std::vector<T> operator* (const std::vector<T> & vec) const
  {
    size_t vec_size = vec.size();
    size_t nx_ = min_sizes.first;
    size_t ny_ = min_sizes.second;

    if (vec_size < nx_ || vec_size < ny_)
    {
      throw std::invalid_argument("size of vector and size of matrix do not match: the dimension of the vector is less than min CSR matrix size");
    }

    std::vector<T> ans(ny_, 0);

    for (size_t i = 0; i < ny_; ++i) 
    {
      size_t row_start = rows[i];
      size_t row_end = rows[i + 1];

      for (size_t j = row_start; j < row_end; ++j) 
      {
        ans[i] += values[j] * vec[cols[j]];
      }
    }

    return ans;
  }

  std::pair<size_t, size_t> Get_matrix_size() const override { return min_sizes; }

private:
  std::pair<size_t, size_t> Get_DOK_sizes(const std::map<std::pair<size_t, size_t>, T> & init_matrix)
  {
    size_t size_x = 0;
    size_t size_y = 0;
    for (const auto & it : init_matrix) 
    {
      size_x = std::max(size_x, it.first.first + 1);
      size_y = std::max(size_y, it.first.second + 1);
    }
    return std::make_pair(size_x, size_y);
  }

  std::vector<T> values;
  std::vector<size_t> cols;
  std::vector<size_t> rows;

  std::pair<size_t, size_t> min_sizes; // real sizes of csr matrix can be bigger than min size because of "0" rows and cols
};

#endif