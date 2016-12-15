
#ifndef SEEDCONSTRAINT_H_
#define SEEDCONSTRAINT_H_


#include "general.h"
#include "IndexRangeList.h"

#include <cstddef>
#include <iostream>

/**
 * Encodes seed constraints to be used for interaction prediction.
 *
 * @author Martin Mann
 *
 */
class SeedConstraint {

public:

	/**
	 * Construction
	 *
	 * @param number of base pairs a seed has to have (>= 2)
	 * @param maxUnpairedOverall_ the maximal (summed) number of unpaired bases
	 *        within both seq1 and seq2 allowed within a seed
	 * @param maxUnpaired1 the maximal number of unpaired bases within seq1
	 *        allowed within a seed
	 * @param maxUnpaired2 the maximal number of unpaired bases within seq2
	 *        allowed within a seed
	 * @param maxE maximal energy a seed is allowed to have
	 * @param ranges1 the index ranges of seq1 to be searched for seeds
	 * @param ranges2reversed the index reversed ranges of seq2 to be searched for seeds
	 */
	SeedConstraint(  const size_t bp
				, const size_t maxUnpairedOverall
				, const size_t maxUnpaired1
				, const size_t maxUnpaired2
				, const E_type maxE
				, const IndexRangeList & ranges1
				, const IndexRangeList & ranges2reversed
				);

	virtual ~SeedConstraint();


	/**
	 * Provides the number of base pairs to be present within a seed
	 *
	 * @return the number of base pairs a seed has to have
	 */
	size_t
	getBasePairs() const;

	/**
	 * Provides the overall maximal number of unpaired bases within a seed
	 *
	 * @return the overall maximal number of unpaired bases within
	 *         a seed is allowed to have
	 */
	size_t
	getMaxUnpairedOverall() const;

	/**
	 * Provides the maximal number of unpaired bases within the first sequence
	 * within a seed
	 *
	 * @return the maximal number of unpaired bases within the first sequence
	 *         a seed is allowed to have
	 */
	size_t
	getMaxUnpaired1() const;

	/**
	 * Provides the maximal number of unpaired bases within the second sequence
	 * within a seed
	 *
	 * @return the maximal number of unpaired bases within the second sequence
	 *         a seed is allowed to have
	 */
	size_t
	getMaxUnpaired2() const;

	/**
	 * Provides the maximally allowed energy for seeds to be considered
	 *
	 * @return the maximally allowed energy for a seed
	 */
	E_type
	getMaxE() const;

	/**
	 * Provides the maximal length of the seed in seq1
	 * @return the maximal length of the seed in seq1
	 */
	size_t
	getMaxLength1() const;

	/**
	 * Provides the maximal length of the seed in seq2
	 * @return the maximal length of the seed in seq2
	 */
	size_t
	getMaxLength2() const;

	/**
	 * Index ranges in seq1 to be searched for seeds or empty if all indices
	 * are to be considered.
	 * @return index ranges to be searched or empty list if all indices relevant
	 */
	const IndexRangeList &
	getRanges1() const;

	/**
	 * Index ranges in seq2 to be searched for seeds or empty if all indices
	 * are to be considered.
	 * @return index ranges to be searched or empty list if all indices relevant
	 */
	const IndexRangeList &
	getRanges2() const;

	/**
	 * Index ranges in seq1 to be searched for seeds or empty if all indices
	 * are to be considered.
	 * @return index ranges to be searched or empty list if all indices relevant
	 */
	IndexRangeList &
	getRanges1();

	/**
	 * Index ranges in seq2 to be searched for seeds or empty if all indices
	 * are to be considered.
	 * @return index ranges to be searched or empty list if all indices relevant
	 */
	IndexRangeList &
	getRanges2();

	/**
	 * Prints the seed constraint details to stream
	 * @param out the ostream to write to
	 * @param c the object to add
	 * @return the altered stream out
	 */
	friend std::ostream& operator<<(std::ostream& out, const SeedConstraint& c);

protected:

	//! the number of base pairs to be present in a seed
	size_t bp;

	//! the overall summed maximally allowed number of unpaired bases in a seed
	size_t maxUnpairedOverall;

	//! the maximally allowed number of unpaired bases in seed seq1
	size_t maxUnpaired1;

	//! the maximally allowed number of unpaired bases in seed seq2
	size_t maxUnpaired2;

	//! the maximal energy allowed for a seed
	E_type maxE;

	//! the index ranges of seq1 to be searched for seeds
	IndexRangeList ranges1;

	//! the index ranges of seq2 to be searched for seeds
	IndexRangeList ranges2;

};

#endif /* SEEDCONSTRAINT_H_ */
