#ifndef OPERATIONS
#define OPERATIONS

#include <vector>
#include "Dense_matrix.h"
#include "CSR_matrix.h"


template <class T>
std::vector<T> operator + (const std::vector<T> & vec1, const std::vector<T> & vec2)
{
  size_t s1 = vec1.size();
  size_t s2 = vec2.size();
  if (s1 != s2)
  {
    throw std::runtime_error("Can't summarize vectors with different sizes");
  }
  std::vector<T> ans(s1);
  for (size_t i = 0; i < s1; ++i)
  {
    ans[i] = vec1[i] + vec2[i];
  }
  return ans;
}

template <class T>
std::vector<T> operator - (const std::vector<T> & vec1, const std::vector<T> & vec2)
{
  size_t s1 = vec1.size();
  size_t s2 = vec2.size();
  if (s1 != s2)
  {
    throw std::runtime_error("Can't substract vectors with different sizes");
  }
  std::vector<T> ans(s1);
  for (size_t i = 0; i < s1; ++i)
  {
    ans[i] = vec1[i] - vec2[i];
  }
  return ans;
}

template <class T>
T operator * (const std::vector<T> & vec1, const std::vector<T> & vec2) 
// Standard scalar product in an orthonormal basis
{
  size_t s1 = vec1.size();
  size_t s2 = vec2.size();
  if (s1 != s2)
  {
    throw std::runtime_error("Can't scalar product vectors with different sizes");
  }
  T ans = 0;
  for (size_t i = 0; i < s1; ++i)
  {
    ans += vec1[i] * vec2[i];
  }
  return ans;
}

template <class T>
std::vector<T> operator * (const std::vector<T> & vec1, const T scalar)
{
  size_t s1 = vec1.size();\
  std::vector<T> ans(s1);
  for (size_t i = 0; i < s1; ++i)
  {
    ans[i] = vec1[i] * scalar;
  }
  return ans;
}

#endif