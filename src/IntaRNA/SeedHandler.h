
#ifndef INTARNA_SEEDHANDLER_H_
#define INTARNA_SEEDHANDLER_H_

#include "IntaRNA/InteractionEnergy.h"
#include "IntaRNA/SeedConstraint.h"

#include <boost/multi_array.hpp>

#include <boost/numeric/ublas/matrix.hpp>

namespace IntaRNA {

/**
 * Abstract interface to provide seed interaction information
 */
class SeedHandler
{

public:

	/**
	 * Construction
	 * @param energy the energy function to be used for seed prediction
	 * @param seedConstraint the seed constraint to be applied
	 */
	SeedHandler(
			const InteractionEnergy & energy
			, const SeedConstraint & seedConstraint
			);

	/**
	 * destruction
	 */
	virtual ~SeedHandler();

	/**
	 * Access to the underlying seed constraint
	 * @return the used seed constraint
	 */
	virtual
	const SeedConstraint&
	getConstraint() const;

	/**
	 * Access to the underlying interaction energy function
	 * @return the used energy function
	 */
	virtual
	const InteractionEnergy&
	getInteractionEnergy() const;


	/**
	 * Computes the seed information for the given interval boundaries
	 * @param i1 the first index of seq1 that might interact
	 * @param j1 the last index of seq1 that might interact
	 * @param i2 the first index of seq2 that might interact
	 * @param j2 the last index of seq2 that might interact
	 * @return the number of potential seed interactions
	 */
	virtual
	size_t
	fillSeed(const size_t i1, const size_t j1, const size_t i2, const size_t j2) = 0;

	/**
	 * Identifies the base pairs of the mfe seed interaction starting at i1,i2
	 * and writes them to the provided container
	 *
	 * NOTE: the left- and right-most base pairs are excluded!
	 *
	 * @param interaction the container to add the base pairs too
	 * @param i1 the start of the seed in seq1
	 * @param i2 the start of the seed in seq2
	 */
	virtual
	void
	traceBackSeed( Interaction & interaction, const size_t i1, const size_t i2) const = 0;


	/**
	 * Access to the mfe of any seed with left-most base pair (i1,i2)
	 * @param i1 the left most interacting base of seq1
	 * @param i2 the left most interacting base of seq2
	 * @return the mfe of any seed starting at (i1,i2) or E_INF if none possible
	 */
	virtual
	E_type
	getSeedE( const size_t i1, const size_t i2 ) const = 0;

	/**
	 * Checks whether or not a given base pair is the left-most base pair of
	 * any seed
	 * @param i1 the interacting base of seq1
	 * @param i2 the interacting base of seq2
	 * @return true if (i1,i2) is the left most base pair of some seed; false
	 *         otherwise
	 */
	virtual
	bool
	isSeedBound( const size_t i1, const size_t i2 ) const = 0;

	/**
	 * Access to the length in seq1 of the mfe seed with left-most base pair (i1,i2)
	 * @param i1 the left most interacting base of seq1
	 * @param i2 the left most interacting base of seq2
	 * @return the length in seq1 of the mfe seed starting at (i1,i2) or 0 if none possible
	 */
	virtual
	size_t
	getSeedLength1( const size_t i1, const size_t i2 ) const = 0;

	/**
	 * Access to the length in seq2 of the mfe seed with left-most base pair (i1,i2)
	 * @param i1 the left most interacting base of seq1
	 * @param i2 the left most interacting base of seq2
	 * @return the length in seq2 of the mfe seed starting at (i1,i2) or 0 if none possible
	 */
	virtual
	size_t
	getSeedLength2( const size_t i1, const size_t i2 ) const = 0;


	/**
	 * Replace the input variables i1 and i2 to values to within the given range
	 * that correspond to
	 *
	 * - the first seed (if the given index pair is no valid seed start or one
	 *   of the indices is out of range bounds)
	 * - the next seed according to some seed order
	 *
	 * The indices are not updated if the last seed within the range is given
	 * or no seed within the range could be found.
	 * It returns whether or not the input variables have been updated.
	 *
	 * Note, if changed, only the seed left-most base pair is within the range
	 * but the full seed indices might exceed i1max or i2max.
	 *
	 * @param i1 seq1 seed index to be changed
	 * @param i2 seq2 seed index to be changed
	 * @param i1min first position within seq1 (inclusive)
	 * @param i1max last position within seq1 (inclusive)
	 * @param i2min first position within seq2 (inclusive)
	 * @param i2max last position within seq2 (inclusive)
	 * @return true if the input variables have been changed; false otherwise
	 */
	virtual
	bool
	updateToNextSeed( size_t & i1, size_t & i2
			, const size_t i1min = 0, const size_t i1max = RnaSequence::lastPos
			, const size_t i2min = 0, const size_t i2max = RnaSequence::lastPos
			) const;

	/**
	 * Adds all seeds to a given interaction that are completely covered by its
	 * base pairs.
	 *
	 * @param i the interaction to update with indexing in global in-sequence
	 *          order
	 */
	virtual
	void
	addSeeds( Interaction & i ) const;

	/**
	 * Checks whether or not two seeds are loop overlapping, ie. given i1 < k1, the
	 * last x (>1) base pairs of seed(i1,i2) are equal to the first x base pairs of
	 * seed(k1,k2), or vice versa if k1 < i1.
	 *
	 * @param i1 the left most interacting base of seq1 of seedA
	 * @param i2 the left most interacting base of seq2 of seedA
	 * @param k1 the left most interacting base of seq1 of seedB
	 * @param k2 the left most interacting base of seq2 of seedB
	 * @return true if seedA and seedB exist and they are loop overlapping;
	 *         false otherwise
	 */
	virtual
	bool
	areLoopOverlapping( const size_t i1, const size_t i2
					, const size_t k1, const size_t k2 ) const;

	/**
	 * Checks whether or not a given index pair is a valid seed base pair
	 *
	 * @param i1 the interacting base of seq1
	 * @param i2 the interacting base of seq2
	 * @param atEndOfSeed whether or not the (i1,i2) form the end of a seed
	 * @return true if (i1,i2) are complementary, ED(i,i) < maxED, and both i1
	 *              and i2 are in allowed regions; false otherwise
	 */
	virtual
	bool
	isFeasibleSeedBasePair( const size_t i1
						, const size_t i2
						, const bool atEndOfSeed = false ) const;

protected:

	//! the used energy function
	const InteractionEnergy& energy;

	//! the seed constraint to be applied
	const SeedConstraint & seedConstraint;



};


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


inline
SeedHandler::SeedHandler(
		const InteractionEnergy & energy
		, const SeedConstraint & seedConstraint
		)
	:
		energy(energy)
		, seedConstraint(seedConstraint)
{
}

////////////////////////////////////////////////////////////////////////////

inline
SeedHandler::~SeedHandler()
{
}

////////////////////////////////////////////////////////////////////////////

inline
const SeedConstraint&
SeedHandler::
getConstraint() const
{
	return seedConstraint;
}

//////////////////////////////////////////////////////////////////////////

inline
const InteractionEnergy&
SeedHandler::
getInteractionEnergy() const
{
	return energy;
}

//////////////////////////////////////////////////////////////////////////

} // namespace

#endif /* SEEDHANDLER_H_ */
