#ifndef CSR_H
#define CSR_H

#include <algorithm>
#include <map>
#include <vector>
#include <stdexcept>
#include "IMatrix.h"

template <class T>
class CSR_Matrix : public IMatrix<T>
{
public:
	CSR_Matrix(const std::map<std::pair<size_t, size_t>, T> & init_matrix, size_t y_size, size_t x_size)
		: sizes(std::make_pair(y_size, x_size))
	{
		size_t init_matrix_size = init_matrix.size();
		values.reserve(init_matrix_size);
		cols.reserve(init_matrix_size);
		rows.resize(y_size + 1, 0);

		for (const auto & element : init_matrix)
		{
			values.push_back(element.second);
			cols.push_back(element.first.second);
			rows[element.first.first + 1]++;
		}

		for (size_t i = 1; i <= y_size; ++i)
		{
			rows[i] += rows[i - 1];
		}
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

	std::vector<T> operator* (const std::vector<T> & vec) const
	{
		size_t vec_size = vec.size();
		size_t y_size = sizes.first; 
		size_t x_size = sizes.second;

		if (vec_size != x_size)
		{
			throw std::invalid_argument("Size of vector and size of matrix do not match");
		}

		std::vector<T> ans(y_size, 0);

		for (size_t i = 0; i < y_size; ++i)
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
