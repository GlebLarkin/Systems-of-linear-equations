#ifndef TRIDIAGONALMATRIX_H
#define TRIDIAGONALMATRIX_H

#include "IMatrix.h"
#include <iostream>


template <class T>
class Tridiagonal_matrix : public IMatrix<T>
{
public:
    Tridiagonal_matrix(const std::vector<T> & upper_diag, const std::vector<T> & main_diag, 
                       const std::vector<T> & lower_diag, const size_t main_diag_size) :
      main_diag_size_(main_diag_size),
      upper_diag_(upper_diag),
      main_diag_(main_diag),
      lower_diag_(lower_diag)
    {
      //add zeroes to make sizes correct
      lower_diag_.insert(lower_diag_.begin(), 0);
      upper_diag_.push_back(0);

      if (main_diag_.size() != main_diag_size)
      {
        throw std::invalid_argument("Wrong size of main diagonal vector\n");
      }
      if (main_diag_.size() != upper_diag_.size())
      {
        throw std::invalid_argument("Wrong size of upper diagonal vector\n");
      }
      if (main_diag_.size() != lower_diag_.size())
      {
        throw std::invalid_argument("Wrong size of lower diagonal vector\n");
      }
    }

    T& operator() (size_t i, size_t j) //for changing matrix  
    {
      if (i >= main_diag_size_ || j >= main_diag_size_)
      {
        throw std::out_of_range("Trying to access an element outside the tridiagonal matrix\n");
      }
      if (i == j) return main_diag_[i];
      if (i == j - 1) return lower_diag_[i];
      if (i == j + 1) return upper_diag_[j];
      else
      {
        throw std::runtime_error("Trying to change an element out of the diagonals in the tridiagonal matrix\n");
      }
    }

    T operator() (size_t i, size_t j) const override //for reading matrix 
    {
      if (i >= main_diag_size_ || j >= main_diag_size_)
      {
        throw std::out_of_range("Trying to access an element outside the tridiagonal matrix\n");
      }
      if (i == j) return main_diag_[i];
      if (i == j - 1) return lower_diag_[i];
      if (i == j + 1) return upper_diag_[j];
      return 0;
    }

    std::pair<size_t, size_t> Get_matrix_size() const override
    {
      return std::make_pair(main_diag_size_, main_diag_size_);
    }

    const std::vector<T> & Get_lower_diag() const { return lower_diag_; }
    const std::vector<T> & Get_main_diag() const { return main_diag_; }
    const std::vector<T> & Get_upper_diag() const { return upper_diag_; }

    std::vector<T> & Get_lower_diag() { return lower_diag_; }
    std::vector<T> & Get_main_diag() { return main_diag_; }
    std::vector<T> & Get_upper_diag() { return upper_diag_; }

private:
    size_t main_diag_size_;
    std::vector<T> upper_diag_;
    std::vector<T> main_diag_;
    std::vector<T> lower_diag_;
};


#endif