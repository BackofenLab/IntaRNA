#ifndef INTARNA_MATRIX_H_
#define INTARNA_MATRIX_H_

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <limits>
#include <mdspan>
#include <stdexcept>
#include <utility>
#include <vector>

namespace IntaRNA {

namespace matrix_detail {
inline std::size_t product(std::size_t rows, std::size_t columns) {
	if (columns != 0 && rows > std::numeric_limits<std::size_t>::max() / columns)
		throw std::length_error("matrix dimensions overflow");
	return rows * columns;
}
inline std::size_t triangle(std::size_t n) {
	if (n == std::numeric_limits<std::size_t>::max())
		throw std::length_error("matrix dimensions overflow");
	return n % 2 == 0 ? product(n / 2, n + 1) : product(n, (n + 1) / 2);
}
}

/** Owning row-major matrix. Views are constructed on access, so copies, moves
 * and resizes cannot leave a cached mdspan pointing into another allocation.
 * resize preserves the overlapping rectangle unless preserve=false is passed.
 * Preserving resize value-initializes new cells. With preserve=false, callers
 * must overwrite all cells before reading. clear resets cells without changing
 * shape.
 */
template<class T>
class Matrix {
	std::vector<T> values;
	std::size_t rows = 0, columns = 0;
	using Extents = std::dextents<std::size_t, 2>;
public:
	using value_type = T;
	Matrix() = default;
	Matrix(std::size_t rows, std::size_t columns, const T &value = T{})
		: values(matrix_detail::product(rows, columns), value), rows(rows), columns(columns) {}
	Matrix(const Matrix &) = default;
	Matrix &operator=(const Matrix &) = default;
	Matrix(Matrix &&other) noexcept : values(std::move(other.values)),
		rows(std::exchange(other.rows, 0)), columns(std::exchange(other.columns, 0)) {}
	Matrix &operator=(Matrix &&other) noexcept {
		if (this != &other) {
			Matrix moved(std::move(other));
			swap(moved);
		}
		return *this;
	}
	std::size_t size1() const noexcept { return rows; }
	std::size_t size2() const noexcept { return columns; }
	std::size_t storageSize() const noexcept { return values.size(); }
	T &operator()(std::size_t i, std::size_t j) {
		assert(i < rows && j < columns);
		return std::mdspan<T, Extents>(values.data(), rows, columns)[i, j];
	}
	const T &operator()(std::size_t i, std::size_t j) const {
		assert(i < rows && j < columns);
		return std::mdspan<const T, Extents>(values.data(), rows, columns)[i, j];
	}
	void clear() { std::fill(values.begin(), values.end(), T{}); }
	void swap(Matrix &other) noexcept {
		values.swap(other.values);
		std::swap(rows, other.rows);
		std::swap(columns, other.columns);
	}
	void resize(std::size_t newRows, std::size_t newColumns, bool preserve = true) {
		if (rows == newRows && columns == newColumns) return;
		if (!preserve) {
			values.resize(matrix_detail::product(newRows, newColumns));
			rows = newRows;
			columns = newColumns;
			return;
		}
		Matrix next(newRows, newColumns);
		for (std::size_t i = 0; i < std::min(rows, newRows); ++i)
			for (std::size_t j = 0; j < std::min(columns, newColumns); ++j)
				next(i, j) = (*this)(i, j);
		swap(next);
	}
};

/** Square upper-triangular matrix with exactly n*(n+1)/2 stored cells.
 * Rows are packed consecutively. A one-dimensional mdspan addresses the
 * packed storage because triangular row lengths are not an affine 2D layout.
 * Const access below the diagonal returns zero; writes require i <= j.
 */
template<class T>
class UpperTriangularMatrix {
	std::vector<T> values;
	std::size_t n = 0;
	using Extents = std::dextents<std::size_t, 1>;
	std::size_t offset(std::size_t i, std::size_t j) const noexcept {
		const auto remaining = n - i;
		// Constructor has already checked that every triangular size fits.
		const auto tail = remaining % 2 == 0
			? (remaining / 2) * (remaining + 1) : remaining * ((remaining + 1) / 2);
		return values.size() - tail + j - i;
	}
	static std::size_t count(std::size_t rows, std::size_t columns) {
		if (rows != columns) throw std::invalid_argument("upper-triangular matrix must be square");
		return matrix_detail::triangle(rows);
	}
public:
	using value_type = T;
	UpperTriangularMatrix() = default;
	UpperTriangularMatrix(std::size_t rows, std::size_t columns)
		: values(count(rows, columns)), n(rows) {}
	UpperTriangularMatrix(const UpperTriangularMatrix &) = default;
	UpperTriangularMatrix &operator=(const UpperTriangularMatrix &) = default;
	UpperTriangularMatrix(UpperTriangularMatrix &&other) noexcept
		: values(std::move(other.values)), n(std::exchange(other.n, 0)) {}
	UpperTriangularMatrix &operator=(UpperTriangularMatrix &&other) noexcept {
		if (this != &other) {
			UpperTriangularMatrix moved(std::move(other));
			swap(moved);
		}
		return *this;
	}
	std::size_t size1() const noexcept { return n; }
	std::size_t size2() const noexcept { return n; }
	std::size_t storageSize() const noexcept { return values.size(); }
	T &operator()(std::size_t i, std::size_t j) {
		assert(i <= j && j < n);
		return std::mdspan<T, Extents>(values.data(), values.size())[offset(i, j)];
	}
	const T &operator()(std::size_t i, std::size_t j) const {
		assert(i < n && j < n);
		static const T zero{};
		return i > j ? zero : std::mdspan<const T, Extents>(values.data(), values.size())[offset(i, j)];
	}
	void clear() { std::fill(values.begin(), values.end(), T{}); }
	void swap(UpperTriangularMatrix &other) noexcept {
		values.swap(other.values);
		std::swap(n, other.n);
	}
	void resize(std::size_t rows, std::size_t columns, bool preserve = true) {
		const auto cells = count(rows, columns);
		if (rows == n) return;
		if (!preserve) {
			values.resize(cells);
			n = rows;
			return;
		}
		UpperTriangularMatrix next(rows, columns);
		for (std::size_t i = 0; i < std::min(n, rows); ++i)
			for (std::size_t j = i; j < std::min(n, rows); ++j)
				next(i, j) = (*this)(i, j);
		swap(next);
	}
};

/** Upper band including the diagonal and 'upper' superdiagonals.
 * Stores rows * min(columns, upper+1) cells, never a dense square for a
 * narrow band. Logical (i,j) maps to mdspan[i,j-i] in the owning Matrix.
 * Const access outside the band returns zero; writes must be in the band.
 */
template<class T>
class UpperBandedMatrix {
	Matrix<T> band;
	std::size_t columns = 0;
	static std::size_t width(std::size_t columns, std::size_t lower, std::size_t upper) {
		if (lower != 0) throw std::invalid_argument("upper-banded matrix requires lower=0");
		return upper >= columns ? columns : upper + 1;
	}
public:
	using value_type = T;
	UpperBandedMatrix() = default;
	UpperBandedMatrix(std::size_t rows, std::size_t columns, std::size_t lower, std::size_t upper)
		: band(rows, width(columns, lower, upper)), columns(columns) {}
	UpperBandedMatrix(const UpperBandedMatrix &) = default;
	UpperBandedMatrix &operator=(const UpperBandedMatrix &) = default;
	UpperBandedMatrix(UpperBandedMatrix &&other) noexcept
		: band(std::move(other.band)), columns(std::exchange(other.columns, 0)) {}
	UpperBandedMatrix &operator=(UpperBandedMatrix &&other) noexcept {
		if (this != &other) {
			UpperBandedMatrix moved(std::move(other));
			swap(moved);
		}
		return *this;
	}
	std::size_t size1() const noexcept { return band.size1(); }
	std::size_t size2() const noexcept { return columns; }
	std::size_t storageSize() const noexcept { return band.storageSize(); }
	T &operator()(std::size_t i, std::size_t j) {
		assert(i < size1() && j < columns && i <= j && j - i < band.size2());
		return band(i, j - i);
	}
	const T &operator()(std::size_t i, std::size_t j) const {
		assert(i < size1() && j < columns);
		static const T zero{};
		return i > j || j - i >= band.size2() ? zero : band(i, j - i);
	}
	void clear() { band.clear(); }
	void swap(UpperBandedMatrix &other) noexcept {
		band.swap(other.band);
		std::swap(columns, other.columns);
	}
	void resize(std::size_t rows, std::size_t newColumns, std::size_t lower,
			std::size_t upper, bool preserve = true) {
		UpperBandedMatrix next(rows, newColumns, lower, upper);
		if (preserve) {
			for (std::size_t i = 0; i < std::min(size1(), rows); ++i)
				for (std::size_t d = 0; d < std::min(band.size2(), next.band.size2())
						&& i < std::min(columns, newColumns)
						&& d < std::min(columns, newColumns) - i; ++d)
					next.band(i, d) = band(i, d);
		}
		swap(next);
	}
};

} // namespace IntaRNA
#endif
