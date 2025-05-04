#ifndef ARNOLDI_H
#define ARNOLDI_H

#include "./Matrixes/Upper_triangular_matrix.h"
#include "./Matrixes/Dense_matrix.h"
#include "./Matrixes/CSR_matrix.h"
#include "./Matrixes/Tridiagonal_matrix.h"
#include "./Operations/Some_matrix_operations.h"
#include "./Operations/vector_operations.h"


template<class T>
struct GivensRotation  
{
  T cos;
  T sin;
  size_t i_;  // every rotation knows what component of vector needs to be rotated
  
  GivensRotation(const std::vector<T> & vec, size_t i) : i_(i) 
  {
    if (vec.size() < 2 || i_ >= vec.size()) 
      throw std::out_of_range("Givens rotation index out of range");

    T a = vec[0];
    T b = vec[i_];
    T r = std::hypot(a, b);

    if (r == T(0)) 
    {
      cos = 1;
      sin = 0;
    } 
    else 
    {
      cos = a / r;
      sin = b / r;
    }
  }

  void rotate(std::vector<T> & vec)
  {
    T v0 = vec[0];
    T vi = vec[i_];

    vec[0] = cos * v0 + sin * vi;
    vec[i_] = -sin * v0 + cos * vi;
  }
};



template<typename T>
struct ArnoldiStepMaker
{
public:
  ArnoldiStepMaker(const UpperTriangularMatrix<T> & R, 
                   const std::vector<GivensRotation<T>> & rotation_vec, 
                   const std::vector<std::vector<T>> & basis)
                   : iteration_(0), R_(R),
                     rotation_vec_(rotation_vec), basis_(basis) {}

    void Apply_all_rotations(std::vector<T> & v)
    {
      for (size_t i = 0; i < rotation_vec_.size(); ++i) { rotation_vec_[i].rotate(v); }
    }

    void Make_step(const CSRMatrix<T> & A)
    {
        std::vector<T> v = A * basis_[iteration_];
        std::vector<T> h(iteration_ + 2);

        for (size_t i = 0; i < iteration_ + 1; ++i)
        {
            h[i] = v * basis_[i];
            v = v - h[i] * basis_[i];
        }

        h[iteration_ + 1] = VectorNorm(v);
        v = v * (1 / h[iteration_ + 1]);
        basis_.emplace_back(v);  // new basis vector is added

        Apply_all_rotations(h);
        GivensRotation<T> new_Givens_rotation(h, h.size() - 1);

        rotation_vec_.emplace_back(new_Givens_rotation);
        new_Givens_rotation.rotate(h);
        R_.addCol(h);
        iteration_++;
    }

    const UpperTriangularMatrix<T> & Get_R() const { return R_; }
    const std::vector<GivensRotation<T>> & Get_rotation_vector() const {return rotation_vec_; }
    const std::vector<std::vector<T>> & Get_basis() const { return basis_; }

private:
  size_t iteration_;
  UpperTriangularMatrix<T> R_;
  std::vector<GivensRotation<T>> rotation_vec_;
  std::vector<std::vector<T>> basis_;
};


#endif