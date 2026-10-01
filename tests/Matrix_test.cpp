#include "catch.hpp"
#include "IntaRNA/Matrix.h"

#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <string>

using namespace IntaRNA;
namespace ublas = boost::numeric::ublas;

namespace {
template<class Actual, class Reference>
void compare(const Actual &actual, const Reference &reference) {
	REQUIRE(actual.size1() == reference.size1());
	REQUIRE(actual.size2() == reference.size2());
	for (std::size_t i = 0; i < actual.size1(); ++i)
		for (std::size_t j = 0; j < actual.size2(); ++j)
			REQUIRE(actual(i, j) == reference(i, j));
}

// Build a reference explicitly: uBLAS's preserving band resize itself can
// throw on rectangular/empty shapes (e.g. 1x2 -> 3x3 -> 0x1 in Boost 1.85).
void resizeReferenceBand(ublas::banded_matrix<int> &matrix, std::size_t rows,
		std::size_t columns, std::size_t upper) {
	ublas::banded_matrix<int> next(rows, columns, 0, upper);
	next.clear();
	for (std::size_t i = 0; i < std::min(rows, matrix.size1()); ++i)
		for (std::size_t j = i; j < std::min(columns, matrix.size2())
				&& j-i <= std::min(upper, matrix.upper()); ++j)
			next(i, j) = matrix(i, j);
	matrix.swap(next);
}

template<class M>
void checkOwnership(M original) {
	original(0, 0) = 7;
	M copy(original);
	copy(0, 0) = 9;
	REQUIRE(original(0, 0) == 7);
	M assigned;
	assigned = copy;
	copy.clear();
	REQUIRE(assigned(0, 0) == 9);
	M moved(std::move(assigned));
	REQUIRE(moved(0, 0) == 9);
	REQUIRE(assigned.size1() == 0);
	REQUIRE(assigned.size2() == 0);
	assigned = std::move(moved);
	REQUIRE(assigned(0, 0) == 9);
	REQUIRE(moved.size1() == 0);
	assigned.swap(original);
	REQUIRE(assigned(0, 0) == 7);
	REQUIRE(original(0, 0) == 9);
}
}

TEST_CASE("Dense matrix preserves logical cells across reshape", "[Matrix]") {
	Matrix<int> actual(3, 5);
	ublas::matrix<int> reference(3, 5);
	for (std::size_t i = 0; i < 3; ++i)
		for (std::size_t j = 0; j < 5; ++j)
			actual(i, j) = reference(i, j) = 10 * i + j;
	compare(actual, reference);
	for (auto shape : {std::pair{5u, 3u}, {2u, 7u}, {2u, 7u}, {0u, 4u}, {4u, 0u}, {1u, 1u}}) {
		const auto oldRows = reference.size1(), oldColumns = reference.size2();
		actual.resize(shape.first, shape.second);
		reference.resize(shape.first, shape.second);
		// uBLAS leaves new primitive cells uninitialized: only retained cells
		// have a defined reference value. Match vector's zero initialization.
		for (std::size_t i = 0; i < reference.size1(); ++i)
			for (std::size_t j = 0; j < reference.size2(); ++j)
				if (i >= oldRows || j >= oldColumns) reference(i, j) = 0;
		compare(actual, reference);
		REQUIRE(actual.storageSize() == shape.first * shape.second);
	}
	actual.resize(2, 3, false);
	actual(1, 2) = 42;
	REQUIRE(std::as_const(actual)(1, 2) == 42);
	actual.clear();
	REQUIRE(actual.size1() == 2);
	REQUIRE(actual(1, 2) == 0);
	checkOwnership(Matrix<int>(2, 3));

	Matrix<std::pair<int, std::string>> objects(2, 2, {42, "seed"});
	auto copy = objects;
	objects.resize(3, 1);
	REQUIRE(objects(1, 0).second == "seed");
	copy(1, 1).second = "helix";
	REQUIRE(objects(1, 0).second == "seed");
}

TEST_CASE("Upper triangular matrix retains packed indexing", "[Matrix]") {
	for (std::size_t n = 0; n <= 12; ++n) {
		UpperTriangularMatrix<int> actual(n, n);
		ublas::triangular_matrix<int, ublas::upper> reference(n, n);
		REQUIRE(actual.storageSize() == n * (n + 1) / 2);
		for (std::size_t i = 0; i < n; ++i)
			for (std::size_t j = i; j < n; ++j)
				actual(i, j) = reference(i, j) = 100 * i + j + 1;
		compare(actual, reference);
		for (auto size : {n + 3, n / 2, std::size_t(0)}) {
			const auto oldSize = reference.size1();
			actual.resize(size, size);
			reference.resize(size, size);
			for (std::size_t i = 0; i < size; ++i)
				for (std::size_t j = std::max(i, oldSize); j < size; ++j)
					reference(i, j) = 0;
			compare(actual, reference);
		}
		actual.resize(n, n, false);
		actual.clear();
		if (n) REQUIRE(actual(n-1, n-1) == 0);
	}
	checkOwnership(UpperTriangularMatrix<int>(3, 3));
	REQUIRE_THROWS_AS(UpperTriangularMatrix<int>(2, 3), std::invalid_argument);
}

TEST_CASE("Upper band retains zeros and the last superdiagonal", "[Matrix]") {
	for (std::size_t rows = 0; rows <= 7; ++rows) {
		for (std::size_t columns = 0; columns <= 7; ++columns) {
			for (std::size_t upper = 0; upper <= 8; ++upper) {
				CAPTURE(rows, columns, upper);
				UpperBandedMatrix<int> actual(rows, columns, 0, upper);
				ublas::banded_matrix<int> reference(rows, columns, 0, upper);
				REQUIRE(actual.storageSize() == rows * std::min(columns, upper + 1));
				for (std::size_t i = 0; i < rows; ++i)
					for (std::size_t j = i; j < columns && j-i <= upper; ++j)
						actual(i, j) = reference(i, j) = 100 * i + j + 1;
				compare(actual, reference);
				actual.resize(rows + 2, columns + 1, 0, upper + 1);
				resizeReferenceBand(reference, rows + 2, columns + 1, upper + 1);
				compare(actual, reference);
				actual.resize(rows / 2, columns / 2, 0, upper / 2);
				resizeReferenceBand(reference, rows / 2, columns / 2, upper / 2);
				compare(actual, reference);
			}
		}
	}
	checkOwnership(UpperBandedMatrix<int>(3, 3, 0, 1));
	UpperBandedMatrix<int> narrow(10000, 10000, 0, 10);
	REQUIRE(narrow.storageSize() == 110000);
	narrow.resize(2, 2, 0, 0, false);
	narrow(1, 1) = 42;
	REQUIRE(std::as_const(narrow)(1, 1) == 42);
	REQUIRE(std::as_const(narrow)(0, 1) == 0);
	REQUIRE_THROWS_AS(UpperBandedMatrix<int>(2, 2, 1, 1), std::invalid_argument);
}

TEST_CASE("Oversized matrix dimensions fail before allocation", "[Matrix]") {
	const auto max = std::numeric_limits<std::size_t>::max();
	REQUIRE_THROWS_AS(Matrix<int>(max, 2), std::length_error);
	REQUIRE_THROWS_AS(UpperTriangularMatrix<int>(max, max), std::length_error);
	REQUIRE_THROWS_AS(UpperTriangularMatrix<int>(max / 2, max / 2), std::length_error);
	REQUIRE_THROWS_AS(UpperBandedMatrix<int>(max, 3, 0, 2), std::length_error);
	Matrix<int> matrix(1, 1, 42);
	REQUIRE_THROWS_AS(matrix.resize(max, 2, false), std::length_error);
	REQUIRE(matrix.size1() == 1);
	REQUIRE(matrix(0, 0) == 42);
}

TEST_CASE("Upper band row views alias valid cells without padding", "[Matrix]") {
	for (std::size_t rows : {0u, 2u, 5u}) {
		for (std::size_t columns : {0u, 2u, 5u}) {
			for (std::size_t upper : {0u, 1u, 5u}) {
				UpperBandedMatrix<int> matrix(rows, columns, 0, upper);
				for (std::size_t i = 0; i < rows; ++i) {
					auto row = matrix.row(i);
					const auto constRow = std::as_const(matrix).row(i);
					REQUIRE(row.size() == (i < columns ? std::min(upper+1, columns-i) : 0));
					REQUIRE(constRow.size() == row.size());
					for (std::size_t k = 0; k < row.size(); ++k) {
						row[k] = static_cast<int>(100*i+k);
						REQUIRE(&constRow[k] == &matrix(i, i+k));
						REQUIRE(constRow[k] == static_cast<int>(100*i+k));
					}
				}
			}
		}
	}
}
