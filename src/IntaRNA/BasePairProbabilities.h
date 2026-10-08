#ifndef INTARNA_BASEPAIRPROBABILITIES_H_
#define INTARNA_BASEPAIRPROBABILITIES_H_
#include "IntaRNA/general.h"
#include "IntaRNA/IndexRange.h"
#include "IntaRNA/Matrix.h"
#include <utility>
#include <vector>

namespace IntaRNA {
/** Success-only accumulator for one target/query pair. Raw numerators and Z are
 * merged over disjoint searched rectangles, then normalized explicitly. All
 * stored coordinates follow the original sequences' 5'-to-3' order. This is
 * conditional on an allowed interaction, without an unbound-state contribution.
 * A destructor never publishes data. Not thread safe: use one owner per pair.
 */
class BasePairProbabilities {
public:
	/** Pending accumulation, successful nonempty/empty completion, or failure. */
	enum class Status { pending, nonempty, empty, failed };
	/** Allocate raw masses for the full original target/query sequences.
	 * @param targetLength number of target positions
	 * @param queryLength number of query positions
	 */
	BasePairProbabilities(size_t targetLength,size_t queryLength);
	/** Commit one completed region, after validating ALL sums and disjointness.
	 * @param target inclusive original target region
	 * @param query inclusive original query region
	 * @param denominator raw region Z
	 * @param masses raw region numerators in original orientation (local indices)
	 * @throws on invalid/overlapping regions or numerical failure, marking failed
	 */
	void addRegion(const IndexRange & target,const IndexRange & query,
			Z_type denominator,const Matrix<Z_type> & masses);
	/** Finalize only after all requested regions succeeded; validate normalized
	 * masses, row/column exclusivity and numerical range before publication.
	 * @throws on failed computation or inconsistent masses
	 */
	void finalize();
	/** Mark this pair failed after any error or cancellation. */
	void fail() noexcept;
	/** @return explicit lifecycle status */
	Status status() const;
	/** @return the raw merged denominator (not usable as a result after failure) */
	Z_type getZ() const;
	/** @return raw merged pair numerators (original coordinates) */
	const Matrix<Z_type> & rawMasses() const;
	/** @return validated matrix; only available for finalized nonempty results.
	 * Successful empty output is represented by status()==empty, not zeroes.
	 */
	Matrix<Z_type> probabilities() const;
private:
	Status state=Status::pending;
	Z_type z=0;
	Matrix<Z_type> mass;
	std::vector<std::pair<IndexRange,IndexRange>> regions;
};
inline BasePairProbabilities::Status BasePairProbabilities::status() const { return state; }
inline Z_type BasePairProbabilities::getZ() const { return z; }
inline const Matrix<Z_type> & BasePairProbabilities::rawMasses() const { return mass; }
inline void BasePairProbabilities::fail() noexcept { state=Status::failed; }
}
#endif
