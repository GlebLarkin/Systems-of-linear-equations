#include <gtest/gtest.h>
#include "Thomas_algorithm.h"
#include "CSR_matrix.h"
#include "Dense_matrix.h"
#include "vector_operations.h"
#include <utility>

//=================Tridiagonal matrix tests=================
TEST(TridiagonalMatrixTest, TestSimpleCase) {
  std::vector<double> lower_diag = {1, 1};
  std::vector<double> main_diag = {3, 3, 3};
  std::vector<double> upper_diag = {1, 1};
  std::vector<double> free_col = {1, 1, 1};

  std::vector<double> expected = {2.0 / 7, 1.0 / 7, 2.0 / 7};

  Tridiagonal_matrix<double> matrix(upper_diag, main_diag, lower_diag, 3);
  Thomas_algorithm<double> solver(matrix, free_col);

  std::vector<double> result = solver.Get_solution();

  ASSERT_EQ(result.size(), expected.size());
  for (size_t i = 0; i < result.size(); ++i) {
    ASSERT_NEAR(result[i], expected[i], 1e-9);
  }
}

TEST(TridiagonalMatrixTest, TestAlsoSimpleCase) {
  std::vector<float> lower_diag = {5, 1};
  std::vector<float> main_diag = {6, 8, 4};
  std::vector<float> upper_diag = {2, 2};
  std::vector<float> free_col = {1, 1, 2};

  std::vector<float> expected = {3.0 / 14, -1.0 / 7, 15.0 / 28};

  Tridiagonal_matrix<float> matrix(upper_diag, main_diag, lower_diag, 3);
  Thomas_algorithm<float> solver(matrix, free_col);

  std::vector<float> result = solver.Get_solution();

  ASSERT_EQ(result.size(), expected.size());
  for (size_t i = 0; i < result.size(); ++i) {
    ASSERT_NEAR(result[i], expected[i], 1e-7);
  }
}


//=================Vector operations tests=================
TEST(VectorOperationsTest, TestVectorAddition) {
  std::vector<double> vec1 = {1.0, 2.0, 3.0};
  std::vector<double> vec2 = {4.0, 5.0, 6.0};
  std::vector<double> expected = {5.0, 7.0, 9.0};

  std::vector<double> result = vec1 + vec2;

  ASSERT_EQ(result.size(), expected.size());
  for (size_t i = 0; i < result.size(); ++i) {
    ASSERT_NEAR(result[i], expected[i], 1e-9);
  }
}

TEST(VectorOperationsTest, TestVectorSubtraction1) {
  std::vector<double> vec1 = {10.1, 20.2, 30.3};
  std::vector<double> vec2 = {4.4, 5.5, 6.6};
  std::vector<double> expected = {5.7, 14.7, 23.7};

  std::vector<double> result = vec1 - vec2;

  ASSERT_EQ(result.size(), expected.size());
  for (size_t i = 0; i < result.size(); ++i) {
    ASSERT_NEAR(result[i], expected[i], 1e-9);
  }
}

TEST(VectorOperationsTest, TestVectorSubtraction2) {
  std::vector<double> vec1 = {10.1, 20.2, 30.3};
  std::vector<double> vec2 = {4.4, 5.5};

  ASSERT_THROW(vec1 - vec2, std::runtime_error);
}

TEST(VectorOperationsTest, TestVectorScalarProduct) {
  std::vector<double> vec1 = {1.0, 2.0, 3.0};
  std::vector<double> vec2 = {4.0, 5.0, 6.0};
  double expected = 32.0;

  double result = vec1 * vec2;

  ASSERT_NEAR(result, expected, 1e-9);
}

TEST(VectorOperationsTest, TestVectorScalarMultiplication) {
  std::vector<double> vec = {1.0, 2.0, 3.0};
  double scalar = 2.0;
  std::vector<double> expected = {2.0, 4.0, 6.0};

  std::vector<double> result = vec * scalar;

  ASSERT_EQ(result.size(), expected.size());
  for (size_t i = 0; i < result.size(); ++i) {
    ASSERT_NEAR(result[i], expected[i], 1e-9);
  }
}

//=================Dense matrix tests=================
TEST(DenseMatrixTest, DenseMatrixTest) {
  std::vector<std::vector<double>> init_matrix = {
    {1.0, 2.0, 3.0},
    {4.0, 5.0, 6.0},
    {7.0, 8.0, 9.0}
  };

  DenseMatrix<double> dense_matrix(3, 3, init_matrix);

  ASSERT_NEAR(dense_matrix(0, 0), 1.0, 1e-9);
  ASSERT_NEAR(dense_matrix(1, 2), 6.0, 1e-9);
  ASSERT_NEAR(dense_matrix(2, 1), 8.0, 1e-9);
  ASSERT_NEAR(dense_matrix(0, 2), 3.0, 1e-9);

  std::vector<double> vec = {1.0, 1.0, 1.0};
  std::vector<double> expected_result = {6.0, 15.0, 24.0};
  std::vector<double> result = dense_matrix * vec;

  ASSERT_EQ(result.size(), expected_result.size());
  for (size_t i = 0; i < result.size(); ++i) {
    ASSERT_NEAR(result[i], expected_result[i], 1e-9);
  }
}

//=================CSR matrix tests=================
TEST(CSRMatrixTest, CSRMatrixTest) {
  std::map<std::pair<size_t, size_t>, double> dok_matrix = {
    {{0, 0}, 1.0},
    {{1, 2}, 2.5},
    {{2, 1}, 3.7},
    {{3, 3}, 4.2}
  };

  CSR_Matrix<double> csr_matrix(dok_matrix, 5, 4);

  ASSERT_NEAR(csr_matrix(0, 0), 1.0, 1e-9);
  ASSERT_NEAR(csr_matrix(1, 2), 2.5, 1e-9);
  ASSERT_NEAR(csr_matrix(2, 1), 3.7, 1e-9);
  ASSERT_NEAR(csr_matrix(3, 3), 4.2, 1e-9);
  ASSERT_NEAR(csr_matrix(1, 1), 0.0, 1e-9);
  ASSERT_NEAR(csr_matrix(4, 3), 0.0, 1e-9);

  std::vector<double> vec = {1.0, 1.0, 1.0, 1.0};
  std::vector<double> expected_result = {1.0, 2.5, 3.7, 4.2};
  std::vector<double> result = csr_matrix * vec;

  ASSERT_EQ(result.size(), expected_result.size());
  for (size_t i = 0; i < result.size(); ++i) {
    ASSERT_NEAR(result[i], expected_result[i], 1e-9);
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
