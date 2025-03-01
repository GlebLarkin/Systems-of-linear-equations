#include "includes.h"

#include <gtest/gtest.h>
#include <utility>


//=================Tridiagonal matrix tests=================
TEST(TridiagonalMatrixTest, TestSimpleCase) 
{
  std::vector<double> lower_diag = {1, 1};
  std::vector<double> main_diag = {3, 3, 3};
  std::vector<double> upper_diag = {1, 1};
  std::vector<double> free_col = {1, 1, 1};

  std::vector<double> expected = {2.0 / 7, 1.0 / 7, 2.0 / 7};

  Tridiagonal_matrix<double> matrix(upper_diag, main_diag, lower_diag, 3);
  Thomas_algorithm<double> solver(matrix, free_col);

  std::vector<double> result = solver.Get_solution();

  ASSERT_EQ(result.size(), expected.size());
  for (size_t i = 0; i < result.size(); ++i) 
  {
    ASSERT_NEAR(result[i], expected[i], 1e-9);
  }
}

TEST(TridiagonalMatrixTest, TestAlsoSimpleCase) 
{
  std::vector<float> lower_diag = {5, 1};
  std::vector<float> main_diag = {6, 8, 4};
  std::vector<float> upper_diag = {2, 2};
  std::vector<float> free_col = {1, 1, 2};

  std::vector<float> expected = {3.0 / 14, -1.0 / 7, 15.0 / 28};

  Tridiagonal_matrix<float> matrix(upper_diag, main_diag, lower_diag, 3);
  Thomas_algorithm<float> solver(matrix, free_col);

  std::vector<float> result = solver.Get_solution();

  ASSERT_EQ(result.size(), expected.size());
  for (size_t i = 0; i < result.size(); ++i) 
  {
    ASSERT_NEAR(result[i], expected[i], 1e-7);
  }
}


//=================Vector operations tests=================
TEST(VectorOperationsTest, TestVectorAddition) 
{
  std::vector<double> vec1 = {1.0, 2.0, 3.0};
  std::vector<double> vec2 = {4.0, 5.0, 6.0};
  std::vector<double> expected = {5.0, 7.0, 9.0};

  std::vector<double> result = vec1 + vec2;

  ASSERT_EQ(result.size(), expected.size());
  for (size_t i = 0; i < result.size(); ++i) 
  {
    ASSERT_NEAR(result[i], expected[i], 1e-9);
  }
}

TEST(VectorOperationsTest, TestVectorSubtraction1) 
{
  std::vector<double> vec1 = {10.1, 20.2, 30.3};
  std::vector<double> vec2 = {4.4, 5.5, 6.6};
  std::vector<double> expected = {5.7, 14.7, 23.7};

  std::vector<double> result = vec1 - vec2;

  ASSERT_EQ(result.size(), expected.size());
  for (size_t i = 0; i < result.size(); ++i) 
  {
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

TEST(VectorOperationsTest, TestVectorScalarMultiplication) 
{
  std::vector<double> vec = {1.0, 2.0, 3.0};
  double scalar = 2.0;
  std::vector<double> expected = {2.0, 4.0, 6.0};

  std::vector<double> result = vec * scalar;

  ASSERT_EQ(result.size(), expected.size());
  for (size_t i = 0; i < result.size(); ++i) 
  {
    ASSERT_NEAR(result[i], expected[i], 1e-9);
  }
}


//=================Dense matrix tests=================
TEST(DenseMatrixTest, DenseMatrixTest) 
{
  std::vector<std::vector<double>> init_matrix = 
  {
    {1.0, 2.0, 3.0},
    {4.0, 5.0, 6.0},
    {7.0, 8.0, 9.0}
  };

  DenseMatrix<double> dense_matrix(init_matrix, 3, 3);

  ASSERT_NEAR(dense_matrix(0, 0), 1.0, 1e-9);
  ASSERT_NEAR(dense_matrix(1, 2), 6.0, 1e-9);
  ASSERT_NEAR(dense_matrix(2, 1), 8.0, 1e-9);
  ASSERT_NEAR(dense_matrix(0, 2), 3.0, 1e-9);

  std::vector<double> vec = {1.0, 1.0, 1.0};
  std::vector<double> expected_result = {6.0, 15.0, 24.0};
  std::vector<double> result = dense_matrix * vec;

  ASSERT_EQ(result.size(), expected_result.size());
  for (size_t i = 0; i < result.size(); ++i) 
  {
    ASSERT_NEAR(result[i], expected_result[i], 1e-9);
  }
}

TEST(DenseMatrixTest, TranspondTest) 
{
  std::vector<std::vector<double>> init_matrix = 
  {
    {1.0, 2.0, 3.0},
    {4.0, 5.0, 6.0},
    {7.0, 8.0, 9.0}
  };

  DenseMatrix<double> dense_matrix(init_matrix, 3, 3);

  auto transponded_matrix = dense_matrix.Transpond();

  ASSERT_NEAR(transponded_matrix(0, 0), 1.0, 1e-9);
  ASSERT_NEAR(transponded_matrix(2, 1), 6.0, 1e-9);
  ASSERT_NEAR(transponded_matrix(1, 2), 8.0, 1e-9);
  ASSERT_NEAR(transponded_matrix(2, 0), 3.0, 1e-9);

}


//=================CSR matrix tests=================
TEST(CSRMatrixTest, CSRMatrixTest) 
{
	std::map<std::pair<size_t, size_t>, double> dok_matrix = 
	{
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
	std::vector<double> expected_result = {1.0, 2.5, 3.7, 4.2, 0.0};
	std::vector<double> result = csr_matrix * vec;

	ASSERT_EQ(result.size(), expected_result.size());
	for (size_t i = 0; i < result.size(); ++i) 
	{
		ASSERT_NEAR(result[i], expected_result[i], 1e-9);
	}
}


//=================QR decomposition tests=================
TEST(QRTest, QRDecomposition) {
  DenseMatrix<double> matrix({
      {12, -51, 4},
      {6, 167, -68},
      {-4, 24, -41}
  }, 3, 3);

  QR<double> qr(matrix);
  DenseMatrix<double> Q = qr.GetQ();
  DenseMatrix<double> R = qr.GetR();

  DenseMatrix<double> QTQ = Q.Transpond() * Q;
  for (size_t i = 0; i < QTQ.Get_matrix_size().first; ++i) {
      for (size_t j = 0; j < QTQ.Get_matrix_size().second; ++j) {
          if (i == j) {
              EXPECT_NEAR(QTQ(i, j), 1.0, 1e-6);
          } else {
              EXPECT_NEAR(QTQ(i, j), 0.0, 1e-6);
          }
      }
  }

  DenseMatrix<double> QR = Q * R;
  for (size_t i = 0; i < matrix.Get_matrix_size().first; ++i) {
      for (size_t j = 0; j < matrix.Get_matrix_size().second; ++j) {
          EXPECT_NEAR(QR(i, j), matrix(i, j), 1e-6);
      }
  }
}

TEST(QRTest, SolveSystem) {
  DenseMatrix<double> matrix({
      {12, -51, 4},
      {6, 167, -68},
      {-4, 24, -41}
  }, 3, 3);

  std::vector<double> b = {1, 0, 0};
  QR<double> qr(matrix);
  std::vector<double> x = qr.solve_system(b);

  std::vector<double> Ax = matrix * x;
  for (size_t i = 0; i < b.size(); ++i) {
      EXPECT_NEAR(Ax[i], b[i], 1e-6);
  }
}


//=================Vector to matrix conversation tests=================
TEST(Vector2MatrixTest, ConvertsVectorCorrectly)
{
  std::vector<double> vec = {1.0, 2.0, 3.0};
  DenseMatrix<double> mat = Vector2matrix(vec);

  EXPECT_EQ(mat.Get_matrix_size().first, 3); 
  EXPECT_EQ(mat.Get_matrix_size().second, 1);

  EXPECT_DOUBLE_EQ(mat(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(mat(1, 0), 2.0);
  EXPECT_DOUBLE_EQ(mat(2, 0), 3.0);
}


//=================Iterative methods tests=================
TEST(IterativeMethodsSolverTest, JacobiMethod)
{
  std::map<std::pair<size_t, size_t>, double> dok_matrix =
  {
    {{0, 0}, 3.0},
    {{0, 1}, 1.0},
    {{0, 2}, 1.0},
    {{1, 0}, 1.0},
    {{1, 1}, 3.0},
    {{1, 2}, 1.0},
    {{2, 0}, 1.0},
    {{2, 1}, 1.0},
    {{2, 2}, 3.0}
  };
  CSR_Matrix<double> A(dok_matrix, 3, 3);

  std::vector<double> test_b = {5, 5, 5};
  std::vector<double> expected_solution = {1, 1, 1};

  IterativeMethodsSolver<double> solver(A, test_b);

  std::vector<double> result = solver.Jacoby_method();

  for (size_t i = 0; i < result.size(); ++i)
  {
    EXPECT_NEAR(result[i], expected_solution[i], 1e-5);
  }
}

TEST(IterativeMethodsSolverTest, GaussSeidelMethod)
{
  std::map<std::pair<size_t, size_t>, double> dok_matrix =
  {
    {{0, 0}, 4.0},
    {{0, 1}, 1.0},
    {{0, 2}, 2.0},
    {{1, 0}, 1.0},
    {{1, 1}, 3.0},
    {{1, 2}, 1.0},
    {{2, 0}, 1.0},
    {{2, 1}, 1.0},
    {{2, 2}, 6.0}
  };
  CSR_Matrix<double> A(dok_matrix, 3, 3);

  std::vector<double> test_b = {7, 5, 8};
  std::vector<double> expected_solution = {1, 1, 1};

  IterativeMethodsSolver<double> solver(A, test_b);

  std::vector<double> result = solver.Gauss_Seidel_method();

  for (size_t i = 0; i < result.size(); ++i)
  {
    EXPECT_NEAR(result[i], expected_solution[i], 1e-5);
  }
}

TEST(IterativeMethodsSolverTest, SimpleIterationMethod)
{
  std::map<std::pair<size_t, size_t>, double> dok_matrix =
  {
    {{0, 0}, 3.0},
    {{0, 1}, 1.0},
    {{0, 2}, 1.0},
    {{1, 0}, 1.0},
    {{1, 1}, 3.0},
    {{1, 2}, 1.0},
    {{2, 0}, 1.0},
    {{2, 1}, 1.0},
    {{2, 2}, 3.0}
  };
  CSR_Matrix<double> A(dok_matrix, 3, 3);

  std::vector<double> test_b = {5, 5, 5};
  std::vector<double> expected_solution = {1, 1, 1};

  IterativeMethodsSolver<double> solver(A, test_b);

  std::vector<double> result = solver.Simple_iteration_method(0.1);

  for (size_t i = 0; i < result.size(); ++i)
  {
    EXPECT_NEAR(result[i], expected_solution[i], 1e-5);
  }
}



int main(int argc, char **argv) 
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
