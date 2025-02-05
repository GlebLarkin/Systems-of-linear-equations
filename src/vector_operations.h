#ifndef OPERATIONS
#define OPERATIONS

#include <vector>


template <class T>
std::vector<T> operator + (const std::vector<T> & vec1, const std::vector<T> & vec2)
{
  size_t s1 = vec1.size();
  size_t s2 = vec2.size();
  if (s1 != s2)
  {
    throw std::runtime_error("Can't summarize vectors with different sizes: " + std::to_string(s1) + " and " + std::to_string(s2));
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
    throw std::runtime_error("Can't substract vectors with different sizes: " + std::to_string(s1) + " and " + std::to_string(s2));
  }
  std::vector<T> ans(s1);
  for (size_t i = 0; i < s1; ++i)
  {
    ans[i] = vec1[i] - vec2[i];
  }
  return ans;
}


#endif