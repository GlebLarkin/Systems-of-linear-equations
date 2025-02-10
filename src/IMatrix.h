#ifndef IMATRIX_H
#define IMATRIX_H

#include <vector>
#include <utility>


template <class T>
class IMatrix
{
public:
  virtual T operator() (size_t i, size_t j) const = 0; // for reading matrix
  virtual std::pair<size_t, size_t> Get_matrix_size() const = 0;
};


#endif