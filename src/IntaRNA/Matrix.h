#ifndef INTARNA_MATRIX_H_
#define INTARNA_MATRIX_H_

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "IntaRNA/intarna_config.h"
#if INTARNA_USE_STD_MDSPAN
#include <mdspan>
#else
#include "mdspan/mdspan.hpp"
#endif

namespace IntaRNA {

namespace matrix_detail {
#if INTARNA_USE_STD_MDSPAN
namespace md = std;
#else
namespace md = MDSPAN_IMPL_STANDARD_NAMESPACE;
#endif

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
	//! Owned contiguous storage.
	std::vector<T> values;
	//! Logical dimensions.
	std::size_t rows = 0, columns = 0;
	using Extents = matrix_detail::md::dextents<std::size_t, 2>;
public:
	//! Type of one stored cell.
	using value_type = T;
	/** Construct an empty matrix. */
	Matrix() = default;
	/**
	 * Allocate a dense matrix initialized to one value.
	 * @param rows logical row count
	 * @param columns logical column count
	 * @param value initial value for each cell
	 * @throws std::length_error if the storage size overflows
	 */
	Matrix(std::size_t rows, std::size_t columns, const T &value = T{});
	/** Copy the shape and values into independent storage. */
	Matrix(const Matrix &) = default;
	/** Copy the shape and values into independent storage.
	 * @return this matrix
	 */
	Matrix &operator=(const Matrix &) = default;
	/**
	 * Take ownership of another matrix and leave it empty.
	 * @param other matrix whose storage is moved
	 */
	Matrix(Matrix &&other) noexcept;
	/**
	 * Take ownership of another matrix; self-move leaves this matrix unchanged.
	 * @param other matrix whose storage is moved; emptied on a non-self move
	 * @return this matrix
	 */
	Matrix &operator=(Matrix &&other) noexcept;
	/**
	 * Number of logical rows.
	 * @return logical row count
	 */
	std::size_t size1() const noexcept;
	/**
	 * Number of logical columns.
	 * @return logical column count
	 */
	std::size_t size2() const noexcept;
	/**
	 * Number of stored elements, including any band padding.
	 * @return stored element count
	 */
	std::size_t storageSize() const noexcept;
	/**
	 * Access a cell within the logical shape.
	 * @param i zero-based row within the logical shape
	 * @param j zero-based column within the logical shape
	 * @return mutable reference to the cell
	 */
	T &operator()(std::size_t i, std::size_t j);
	/**
	 * Access a cell within the logical shape.
	 * @param i zero-based row within the logical shape
	 * @param j zero-based column within the logical shape
	 * @return read-only reference to the cell
	 */
	const T &operator()(std::size_t i, std::size_t j) const;
	/**
	 * Reset stored cells to default values without changing the shape.
	 */
	void clear();
	/**
	 * Exchange shape and owned storage.
	 * @param other matrix to exchange with
	 */
	void swap(Matrix &other) noexcept;
	/**
	 * Resize, preserving overlapping logical cells when requested.
	 * @param newRows new row count
	 * @param newColumns new column count
	 * @param preserve retain overlapping cells and default-initialize new cells;
	 * otherwise callers must overwrite cells before reading
	 * @throws std::length_error if the storage size overflows
	 */
	void resize(std::size_t newRows, std::size_t newColumns, bool preserve = true);
};

template<class T>
inline
Matrix<T>::Matrix(std::size_t rows, std::size_t columns, const T &value)
	: values(matrix_detail::product(rows, columns), value), rows(rows), columns(columns)
{}

template<class T>
inline
Matrix<T>::Matrix(Matrix &&other) noexcept : values(std::move(other.values)),
	rows(std::exchange(other.rows, 0)), columns(std::exchange(other.columns, 0))
{}

template<class T>
inline
Matrix<T> &Matrix<T>::operator=(Matrix &&other) noexcept
{
	if (this != &other) {
		Matrix moved(std::move(other));
		swap(moved);
	}
	return *this;
}

template<class T>
inline
std::size_t Matrix<T>::size1() const noexcept
{ return rows; }

template<class T>
inline
std::size_t Matrix<T>::size2() const noexcept
{ return columns; }

template<class T>
inline
std::size_t Matrix<T>::storageSize() const noexcept
{ return values.size(); }

template<class T>
inline
T &Matrix<T>::operator()(std::size_t i, std::size_t j)
{
	assert(i < rows && j < columns);
	return matrix_detail::md::mdspan<T, Extents>(values.data(), rows, columns)[i, j];
}

template<class T>
inline
const T &Matrix<T>::operator()(std::size_t i, std::size_t j) const
{
	assert(i < rows && j < columns);
	return matrix_detail::md::mdspan<const T, Extents>(values.data(), rows, columns)[i, j];
}

template<class T>
inline
void Matrix<T>::clear()
{ std::fill(values.begin(), values.end(), T{}); }

template<class T>
inline
void Matrix<T>::swap(Matrix &other) noexcept
{
	values.swap(other.values);
	std::swap(rows, other.rows);
	std::swap(columns, other.columns);
}

template<class T>
inline
void Matrix<T>::resize(std::size_t newRows, std::size_t newColumns, bool preserve)
{
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


/** Square upper-triangular matrix with exactly n*(n+1)/2 stored cells.
 * Rows are packed consecutively. A one-dimensional mdspan addresses the
 * packed storage because triangular row lengths are not an affine 2D layout.
 * Const access below the diagonal returns zero; writes require i <= j.
 */
template<class T>
class UpperTriangularMatrix {
	//! Owned contiguous storage.
	std::vector<T> values;
	//! Logical row and column count.
	std::size_t n = 0;
	using Extents = matrix_detail::md::dextents<std::size_t, 1>;
	/**
	 * Offset of a stored upper-triangular cell in the packed allocation.
	 */
	std::size_t offset(std::size_t i, std::size_t j) const noexcept;
	/**
	 * Validate a square shape and return its packed storage size.
	 */
	static std::size_t count(std::size_t rows, std::size_t columns);
public:
	//! Type of one stored cell.
	using value_type = T;
	/** Construct an empty matrix. */
	UpperTriangularMatrix() = default;
	/**
	 * Allocate a default-initialized square upper triangle.
	 * @param rows logical row count
	 * @param columns logical column count; must equal rows
	 * @throws std::invalid_argument if the shape is not square
	 * @throws std::length_error if the storage size overflows
	 */
	UpperTriangularMatrix(std::size_t rows, std::size_t columns);
	/** Copy the shape and values into independent storage. */
	UpperTriangularMatrix(const UpperTriangularMatrix &) = default;
	/** Copy the shape and values into independent storage.
	 * @return this matrix
	 */
	UpperTriangularMatrix &operator=(const UpperTriangularMatrix &) = default;
	/**
	 * Take ownership of another matrix and leave it empty.
	 * @param other matrix whose storage is moved
	 */
	UpperTriangularMatrix(UpperTriangularMatrix &&other) noexcept;
	/**
	 * Take ownership of another matrix; self-move leaves this matrix unchanged.
	 * @param other matrix whose storage is moved; emptied on a non-self move
	 * @return this matrix
	 */
	UpperTriangularMatrix &operator=(UpperTriangularMatrix &&other) noexcept;
	/**
	 * Number of logical rows.
	 * @return logical row count
	 */
	std::size_t size1() const noexcept;
	/**
	 * Number of logical columns.
	 * @return logical column count
	 */
	std::size_t size2() const noexcept;
	/**
	 * Number of stored elements, including any band padding.
	 * @return stored element count
	 */
	std::size_t storageSize() const noexcept;
	/**
	 * Access a cell with i <= j < size2().
	 * @param i zero-based row within the logical shape
	 * @param j zero-based column within the logical shape
	 * @return mutable reference to the cell
	 */
	T &operator()(std::size_t i, std::size_t j);
	/**
	 * Read a logical cell; structural zeros are returned outside stored cells.
	 * @param i zero-based row within the logical shape
	 * @param j zero-based column within the logical shape
	 * @return read-only reference to the cell
	 */
	const T &operator()(std::size_t i, std::size_t j) const;
	/**
	 * Reset stored cells to default values without changing the shape.
	 */
	void clear();
	/**
	 * Exchange shape and owned storage.
	 * @param other matrix to exchange with
	 */
	void swap(UpperTriangularMatrix &other) noexcept;
	/**
	 * Resize, preserving overlapping logical cells when requested.
	 * @param rows new row count
	 * @param columns new column count; must equal rows
	 * @param preserve retain overlapping cells and default-initialize new cells;
	 * otherwise callers must overwrite cells before reading
	 * @throws std::length_error if the storage size overflows
	 * @throws std::invalid_argument if the requested shape or band is unsupported
	 */
	void resize(std::size_t rows, std::size_t columns, bool preserve = true);
};

template<class T>
inline
std::size_t UpperTriangularMatrix<T>::offset(std::size_t i, std::size_t j) const noexcept
{
	const auto remaining = n - i;
	// Constructor has already checked that every triangular size fits.
	const auto tail = remaining % 2 == 0
		? (remaining / 2) * (remaining + 1) : remaining * ((remaining + 1) / 2);
	return values.size() - tail + j - i;
}

template<class T>
inline
std::size_t UpperTriangularMatrix<T>::count(std::size_t rows, std::size_t columns)
{
	if (rows != columns) throw std::invalid_argument("upper-triangular matrix must be square");
	return matrix_detail::triangle(rows);
}

template<class T>
inline
UpperTriangularMatrix<T>::UpperTriangularMatrix(std::size_t rows, std::size_t columns)
	: values(count(rows, columns)), n(rows)
{}

template<class T>
inline
UpperTriangularMatrix<T>::UpperTriangularMatrix(UpperTriangularMatrix &&other) noexcept
	: values(std::move(other.values)), n(std::exchange(other.n, 0))
{}

template<class T>
inline
UpperTriangularMatrix<T> &UpperTriangularMatrix<T>::operator=(UpperTriangularMatrix &&other) noexcept
{
	if (this != &other) {
		UpperTriangularMatrix moved(std::move(other));
		swap(moved);
	}
	return *this;
}

template<class T>
inline
std::size_t UpperTriangularMatrix<T>::size1() const noexcept
{ return n; }

template<class T>
inline
std::size_t UpperTriangularMatrix<T>::size2() const noexcept
{ return n; }

template<class T>
inline
std::size_t UpperTriangularMatrix<T>::storageSize() const noexcept
{ return values.size(); }

template<class T>
inline
T &UpperTriangularMatrix<T>::operator()(std::size_t i, std::size_t j)
{
	assert(i <= j && j < n);
	return matrix_detail::md::mdspan<T, Extents>(values.data(), values.size())[offset(i, j)];
}

template<class T>
inline
const T &UpperTriangularMatrix<T>::operator()(std::size_t i, std::size_t j) const
{
	assert(i < n && j < n);
	static const T zero{};
	return i > j ? zero : matrix_detail::md::mdspan<const T, Extents>(values.data(), values.size())[offset(i, j)];
}

template<class T>
inline
void UpperTriangularMatrix<T>::clear()
{ std::fill(values.begin(), values.end(), T{}); }

template<class T>
inline
void UpperTriangularMatrix<T>::swap(UpperTriangularMatrix &other) noexcept
{
	values.swap(other.values);
	std::swap(n, other.n);
}

template<class T>
inline
void UpperTriangularMatrix<T>::resize(std::size_t rows, std::size_t columns, bool preserve)
{
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


/** Upper band including the diagonal and 'upper' superdiagonals.
 * Stores rows * min(columns, upper+1) cells, never a dense square for a
 * narrow band. Logical (i,j) maps to mdspan[i,j-i] in the owning Matrix.
 * Const access outside the band returns zero; writes must be in the band.
 */
template<class T>
class UpperBandedMatrix {
	//! Physical rows of the upper band, including row-end padding.
	Matrix<T> band;
	//! Logical column count, independent of stored band width.
	std::size_t columns = 0;
	/**
	 * Validate the lower bandwidth and clamp the stored width to the columns.
	 */
	static std::size_t width(std::size_t columns, std::size_t lower, std::size_t upper);
public:
	//! Type of one stored cell.
	using value_type = T;
	/** Construct an empty matrix. */
	UpperBandedMatrix() = default;
	/**
	 * Allocate a default-initialized upper band.
	 * @param rows logical row count
	 * @param columns logical column count
	 * @param lower lower bandwidth; must be zero
	 * @param upper number of stored superdiagonals
	 * @throws std::invalid_argument if lower is not zero
	 * @throws std::length_error if the storage size overflows
	 */
	UpperBandedMatrix(std::size_t rows, std::size_t columns, std::size_t lower, std::size_t upper);
	/** Copy the shape and values into independent storage. */
	UpperBandedMatrix(const UpperBandedMatrix &) = default;
	/** Copy the shape and values into independent storage.
	 * @return this matrix
	 */
	UpperBandedMatrix &operator=(const UpperBandedMatrix &) = default;
	/**
	 * Take ownership of another matrix and leave it empty.
	 * @param other matrix whose storage is moved
	 */
	UpperBandedMatrix(UpperBandedMatrix &&other) noexcept;
	/**
	 * Take ownership of another matrix; self-move leaves this matrix unchanged.
	 * @param other matrix whose storage is moved; emptied on a non-self move
	 * @return this matrix
	 */
	UpperBandedMatrix &operator=(UpperBandedMatrix &&other) noexcept;
	/**
	 * Number of logical rows.
	 * @return logical row count
	 */
	std::size_t size1() const noexcept;
	/**
	 * Number of logical columns.
	 * @return logical column count
	 */
	std::size_t size2() const noexcept;
	/**
	 * Number of stored elements, including any band padding.
	 * @return stored element count
	 */
	std::size_t storageSize() const noexcept;
	/**
	 * Access a cell within the stored upper band.
	 * @param i zero-based row within the logical shape
	 * @param j zero-based column within the logical shape
	 * @return mutable reference to the cell
	 */
	T &operator()(std::size_t i, std::size_t j);
	/**
	 * Read a logical cell; structural zeros are returned outside stored cells.
	 * @param i zero-based row within the logical shape
	 * @param j zero-based column within the logical shape
	 * @return read-only reference to the cell
	 */
	const T &operator()(std::size_t i, std::size_t j) const;
	/**
	 * View the contiguous stored cells (i,i), (i,i+1), ... in one row.
	 * Excludes structural zeros and row-end padding. The view aliases this
	 * matrix and must not outlive its storage or be retained across resize,
	 * move, swap or assignment.
	 * @param i zero-based row, less than size1()
	 * @return mutable stored cells; empty when i >= size2() or the band is empty
	 */
	std::span<T> row(std::size_t i);
	/**
	 * Read-only view of one stored row, with the same bounds and lifetime as row().
	 * @param i zero-based row, less than size1()
	 * @return read-only stored cells, excluding padding and structural zeros
	 */
	std::span<const T> row(std::size_t i) const;
	/**
	 * Reset stored cells to default values without changing the shape.
	 */
	void clear();
	/**
	 * Exchange shape and owned storage.
	 * @param other matrix to exchange with
	 */
	void swap(UpperBandedMatrix &other) noexcept;
	/**
	 * Resize, preserving overlapping logical cells when requested.
	 * @param rows new row count
	 * @param newColumns new column count
	 * @param lower lower bandwidth; must be zero
	 * @param upper number of stored superdiagonals
	 * @param preserve retain overlapping cells and default-initialize new cells;
	 * otherwise callers must overwrite cells before reading
	 * @throws std::length_error if the storage size overflows
	 * @throws std::invalid_argument if the requested shape or band is unsupported
	 */
	void resize(std::size_t rows, std::size_t newColumns, std::size_t lower,
			std::size_t upper, bool preserve = true);
};

template<class T>
inline
std::size_t UpperBandedMatrix<T>::width(std::size_t columns, std::size_t lower, std::size_t upper)
{
	if (lower != 0) throw std::invalid_argument("upper-banded matrix requires lower=0");
	return upper >= columns ? columns : upper + 1;
}

template<class T>
inline
UpperBandedMatrix<T>::UpperBandedMatrix(std::size_t rows, std::size_t columns, std::size_t lower, std::size_t upper)
	: band(rows, width(columns, lower, upper)), columns(columns)
{}

template<class T>
inline
UpperBandedMatrix<T>::UpperBandedMatrix(UpperBandedMatrix &&other) noexcept
	: band(std::move(other.band)), columns(std::exchange(other.columns, 0))
{}

template<class T>
inline
UpperBandedMatrix<T> &UpperBandedMatrix<T>::operator=(UpperBandedMatrix &&other) noexcept
{
	if (this != &other) {
		UpperBandedMatrix moved(std::move(other));
		swap(moved);
	}
	return *this;
}

template<class T>
inline
std::size_t UpperBandedMatrix<T>::size1() const noexcept
{ return band.size1(); }

template<class T>
inline
std::size_t UpperBandedMatrix<T>::size2() const noexcept
{ return columns; }

template<class T>
inline
std::size_t UpperBandedMatrix<T>::storageSize() const noexcept
{ return band.storageSize(); }

template<class T>
inline
T &UpperBandedMatrix<T>::operator()(std::size_t i, std::size_t j)
{
	assert(i < size1() && j < columns && i <= j && j - i < band.size2());
	return band(i, j - i);
}

template<class T>
inline
const T &UpperBandedMatrix<T>::operator()(std::size_t i, std::size_t j) const
{
	assert(i < size1() && j < columns);
	static const T zero{};
	return i > j || j - i >= band.size2() ? zero : band(i, j - i);
}

template<class T>
inline
std::span<T> UpperBandedMatrix<T>::row(std::size_t i)
{
	assert(i < size1());
	const auto count = i < columns ? std::min(band.size2(), columns-i) : 0;
	return count == 0 ? std::span<T>{} : std::span<T>{&band(i, 0), count};
}

template<class T>
inline
std::span<const T> UpperBandedMatrix<T>::row(std::size_t i) const
{
	assert(i < size1());
	const auto count = i < columns ? std::min(band.size2(), columns-i) : 0;
	return count == 0 ? std::span<const T>{} : std::span<const T>{&band(i, 0), count};
}

template<class T>
inline
void UpperBandedMatrix<T>::clear()
{ band.clear(); }

template<class T>
inline
void UpperBandedMatrix<T>::swap(UpperBandedMatrix &other) noexcept
{
	band.swap(other.band);
	std::swap(columns, other.columns);
}

template<class T>
inline
void UpperBandedMatrix<T>::resize(std::size_t rows, std::size_t newColumns, std::size_t lower,
		std::size_t upper, bool preserve)
{
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


} // namespace IntaRNA
#endif
