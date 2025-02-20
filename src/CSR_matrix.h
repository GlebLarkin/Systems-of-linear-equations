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
  CSR_Matrix(const std::map<std::pair<size_t, size_t>, T> & user_matrix, size_t x_size, size_t y_size) 
    : sizes(std::make_pair(x_size, y_size)) 
    {
      auto init_matrix = user_matrix;

      auto last_el = init_matrix.find(std::make_pair(x_size - 1, y_size - 1));
      if (last_el == init_matrix.end())
      {
        init_matrix.insert({std::make_pair(x_size - 1, y_size - 1), T(0)}); // sizes fix
      }

      // Reserve memory
      size_t init_matrix_size = init_matrix.size();
      values.reserve(init_matrix_size);
      cols.reserve(init_matrix_size);
      rows.resize(sizes.first + 1, 0);

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

      rows[sizes.first] = values.size();
    }

  T operator() (size_t i, size_t j) const override
  {
    if (i >= sizes.first || j >= sizes.second) throw std::out_of_range("Trying to get index out of matrix size");

    size_t value_row_begin = rows[i];
    size_t value_row_end = rows[i + 1];

    for (size_t k = value_row_begin; k < value_row_end; ++k) 
    {
      if (cols[k] == j) return values[k];
    }

    return T(0);
  }

  //T& operator() (size_t i, size_t j) = delete;
  // hard realization

  std::vector<T> operator* (const std::vector<T> & vec) const
  {
    size_t vec_size = vec.size();
    size_t nx_ = sizes.first;
    size_t ny_ = sizes.second;

    if (vec_size != ny_)
    {
      throw std::invalid_argument("size of vector and size of matrix do not match");
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

  std::pair<size_t, size_t> Get_matrix_size() const override { return sizes; }

private:
  std::vector<T> values;
  std::vector<size_t> cols;
  std::vector<size_t> rows;

  std::pair<size_t, size_t> sizes;
};

#endif