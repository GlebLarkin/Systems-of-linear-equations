#ifndef GENERATE_ELLIPTIC_MATRIX_H
#define GENERATE_ELLIPTIC_MATRIX_H

#include "./Matrixes/CSR_matrix.h"
#include <map>


template <class T>
CSR_Matrix<T> Generate_elliptic_matrix(size_t N)
{
  std::map<std::pair<size_t, size_t>, T> DOK;

  for (size_t i = 0; i < N * N; ++i)
  {
    DOK[{i, i}] = T(4);
  }

  for (size_t i = 0; i < N * N - (N - 1); i += N)
  {
    for (size_t j = 0; j < N - 1; ++j)
    {
      DOK[{i + j, i + j + 1}] = T(-1);
      DOK[{i + j + 1, i + j}] = T(-1);
    }
  }

  for (size_t i = 0; i < N * N - N; ++i)
  {
    DOK[{i, i + N}] = T(-1);
    DOK[{i + N, i}] = T(-1);
  }

  return CSR_Matrix<T>(DOK, N * N, N * N);
}


#endif
