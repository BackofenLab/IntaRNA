
#ifndef SEEDHANDLERIDXOFFSET_H_
#define SEEDHANDLERIDXOFFSET_H_

#include "SeedHandler.h"

#include <boost/multi_array.hpp>

#include <boost/numeric/ublas/matrix.hpp>


/**
 * Wrapper around a SeedHandler object where indices are shifted by
 * a given positive offset (shifted towards infinity).
 * This is useful for local interaction computations.
 *
 * @author Martin Mann
 *
 */
class SeedHandlerIdxOffset
{

public:

	/**
	 * Construction
	 * @param energy the energy function to be used for seed prediction (already offset)
	 * @param seedConstraint the seed constraint to be applied
	 */
	SeedHandlerIdxOffset(
			const InteractionEnergy & energy
			, const SeedConstraint & seedConstraint
			);

	/**
	 * destruction
	 */
	virtual ~SeedHandlerIdxOffset();

	/**
	 * Access to the currently used index offset for sequence 1
	 * @return the index offset for sequence 1 used
	 */
	size_t getOffset1() const;

	/**
	 * Sets the index offset to be used for sequence 1
	 * @param offset1 the index offset for sequence 1 to be used
	 */
	void setOffset1(size_t offset1);

	/**
	 * Access to the currently used index offset for sequence 2
	 * @return the index offset for sequence 2 used
	 */
	size_t getOffset2() const;

	/**
	 * Sets the index offset to be used for sequence 2
	 * @param offset1 the index offset for sequence 2 to be used
	 */
	void setOffset2(size_t offset2);


	///////////////  OVERWRITTEN MEMBERS USING OFFSET  /////////////////


	/**
	 * Access to the underlying seed constraint
	 * @return the used seed constraint
	 */
	virtual
	const SeedConstraint&
	getConstraint() const;


	/**
	 * Computes the seed matrix for the given interval boundaries
	 *
	 * Note, the indices are shifted by an offset for computation.
	 *
	 * @param i1 the first index of seq1 that might interact
	 * @param j1 the last index of seq1 that might interact
	 * @param i2 the first index of seq2 that might interact
	 * @param j2 the last index of seq2 that might interact
	 */
	virtual
	void
	fillSeed(const size_t i1, const size_t j1, const size_t i2, const size_t j2);

	/**
	 * Identifies the base pairs of the mfe seed interaction starting at i1,i2
	 * and writes them to the provided container
	 *
	 * NOTE: the right most base pair is excluded!
	 *
	 * Note, the indices are shifted by an offset for computation.
	 *
	 * @param interaction the container to add the base pairs too
	 * @param i1 the start of the seed in seq1
	 * @param i2 the start of the seed in seq2
	 */
	virtual
	void
	traceBackSeed( Interaction & interaction, const size_t i1, const size_t i2);


	/**
	 * Access to the mfe of any seed with left-most base pair (i1,i2)
	 *
	 * Note, the indices are shifted by an offset for computation.
	 *
	 * @param i1 the left most interacting base of seq1
	 * @param i2 the left most interacting base of seq2
	 * @return the mfe of any seed starting at (i1,i2) or E_INF if none possible
	 */
	virtual
	E_type
	getSeedE( const size_t i1, const size_t i2 ) const;

	/**
	 * Access to the length in seq1 of the mfe seed with left-most base pair (i1,i2)
	 *
	 * Note, the indices are shifted by an offset for computation.
	 *
	 * @param i1 the left most interacting base of seq1
	 * @param i2 the left most interacting base of seq2
	 * @return the length in seq1 of the mfe seed starting at (i1,i2) or 0 if none possible
	 */
	virtual
	size_t
	getSeedLength1( const size_t i1, const size_t i2 ) const;

	/**
	 * Access to the length in seq2 of the mfe seed with left-most base pair (i1,i2)
	 *
	 * Note, the indices are shifted by an offset for computation.
	 *
	 * @param i1 the left most interacting base of seq1
	 * @param i2 the left most interacting base of seq2
	 * @return the length in seq2 of the mfe seed starting at (i1,i2) or 0 if none possible
	 */
	virtual
	size_t
	getSeedLength2( const size_t i1, const size_t i2 ) const;



protected:

	//! the index shifted seed constraint
	SeedHandler seedHandlerOriginal;

	//! the index shifted seed constraint
	SeedConstraint seedConstraintOffset;

	//! offset for indices in seq1
	size_t idxOffset1;

	//! offset for indices in seq2
	size_t idxOffset2;

};

#endif /* SEEDHANDLERIDXOFFSET_H_ */
