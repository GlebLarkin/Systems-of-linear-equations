#include "Solver.h"

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
















// Дополнительная функция для сравнения двух матриц с заданной точностью
template <typename T>
bool MatricesApproxEqual(const DenseMatrix<T>& A, const DenseMatrix<T>& B, T tol)
{
  auto [rowsA, colsA] = A.Get_matrix_size();
  auto [rowsB, colsB] = B.Get_matrix_size();
  if (rowsA != rowsB || colsA != colsB)
  {
    return false;
  }
  for (size_t i = 0; i < rowsA; i++)
  {
    for (size_t j = 0; j < colsA; j++)
    {
      if (std::abs(A(i, j) - B(i, j)) > tol)
      {
        return false;
      }
    }
  }
  return true;
}

// Функция для проверки ортогональности матрицы Q: Q^T * Q == I
template <typename T>
bool IsOrthogonal(const DenseMatrix<T>& Q, T tol)
{
  DenseMatrix<T> QT = Q.Transpond();
  DenseMatrix<T> I_calc = QT * Q; // Должно получиться I, размерности Q.Get_matrix_size().first x Q.Get_matrix_size().first
  auto [m, n] = Q.Get_matrix_size();
  std::vector<std::vector<T>> ident(m, std::vector<T>(m, T(0)));
  for (size_t i = 0; i < m; i++)
  {
    ident[i][i] = T(1);
  }
  DenseMatrix<T> I(ident, m, m);
  return MatricesApproxEqual(I_calc, I, tol);
}

// Тест 1: Разложение для высокой матрицы 3x2
TEST(QRDecompositionTest, TallMatrixDecomposition)
{
  // Исходная матрица A (3x2)
  std::vector<std::vector<double>> A_data = { {12, -51},
                                              { 6, 167},
                                              {-4, 24} };
  DenseMatrix<double> A(A_data, 3, 2);

  // Выполняем QR-разложение
  QR<double> qr(A);
  DenseMatrix<double> Q = qr.GetQ();
  DenseMatrix<double> R = qr.GetR();

  // Проверяем ортогональность Q: Q^T * Q = I (размер 3x3)
  EXPECT_TRUE(IsOrthogonal(Q, 1e-6));

  // Проверяем, что Q * R приблизительно равна A
  DenseMatrix<double> A_reconstructed = Q * R;
  EXPECT_TRUE(MatricesApproxEqual(A, A_reconstructed, 1e-6));
}

// Тест 2: Разложение для квадратной матрицы 2x2
TEST(QRDecompositionTest, SquareMatrixDecomposition)
{
  std::vector<std::vector<double>> A_data = { {1, 2},
                                              {3, 4} };
  DenseMatrix<double> A(A_data, 2, 2);

  QR<double> qr(A);
  DenseMatrix<double> Q = qr.GetQ();
  DenseMatrix<double> R = qr.GetR();

  EXPECT_TRUE(IsOrthogonal(Q, 1e-6));
  DenseMatrix<double> A_reconstructed = Q * R;
  EXPECT_TRUE(MatricesApproxEqual(A, A_reconstructed, 1e-6));
}

// Тест 3: Проверка, что для широкой матрицы (m < n) выбрасывается исключение
TEST(QRDecompositionTest, WideMatrixThrows)
{
  std::vector<std::vector<double>> A_data = { {1, 2, 3},
                                              {4, 5, 6} }; // 2x3 матрица (широкая)
  DenseMatrix<double> A(A_data, 2, 3);
  EXPECT_THROW(QR<double> qr(A), std::invalid_argument);
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



int main(int argc, char **argv) 
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
