
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
	traceBackSeed( Interaction & interaction, const size_t i1, const size_t i2) = 0;


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

protected:

	//! the used energy function
	const InteractionEnergy& energy;

	//! the seed constraint to be applied
	const SeedConstraint & seedConstraint;


protected:

	/**
	 * Checks whether or not a given index pair is a valid seed base pair
	 *
	 * @param i1 the interacting base of seq1
	 * @param i2 the interacting base of seq2
	 * @return true if (i1,i2) are complementary, ED(i,i) < maxED, and both i1
	 *              and i2 are in allowed regions; false otherwise
	 */
	bool
	isFeasibleSeedBasePair( const size_t i1, const size_t i2 ) const;


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

inline
bool
SeedHandler::
isFeasibleSeedBasePair( const size_t i1, const size_t i2 ) const
{
#if INTARNA_IN_DEBUG_MODE
	if ( i1 >= energy.size1() ) throw std::runtime_error("SeedHandler::isFeasibleSeedBasePair: i1("+toString(i1)+") >= energy.size1("+toString(energy.size1())+")");
	if ( i2 >= energy.size2() ) throw std::runtime_error("SeedHandler::isFeasibleSeedBasePair: i2("+toString(i2)+") >= energy.size2("+toString(energy.size2())+")");
#endif

	return		i1 < energy.size1() && i2 < energy.size2()
			&&	energy.isAccessible1(i1)
			&&	energy.isAccessible2(i2)
			&&	energy.areComplementary(i1,i2)
			&&	seedConstraint.getMaxED() >= energy.getED1( i1,i1 )
			&&	seedConstraint.getMaxED() >= energy.getED2( i2,i2 )
			&&  (seedConstraint.isGUallowed() || (!energy.isGU( i1, i2 ))) // check for GU bp constraint
			&&	(seedConstraint.getRanges1().empty() || seedConstraint.getRanges1().covers(i1))
			&&	(seedConstraint.getRanges2().empty() || seedConstraint.getRanges2().covers(i2))
			;
}

//////////////////////////////////////////////////////////////////////////

inline
bool
SeedHandler::
updateToNextSeed( size_t & i1_out, size_t & i2_out
		, const size_t i1min, const size_t i1max
		, const size_t i2min, const size_t i2max
		) const
{
	size_t i1=i1_out, i2=i2_out;
	// find true max value
	const size_t i1maxVal = std::min(energy.size1()-1,i1max)
				, i2maxVal = std::min(energy.size2()-1,i2max);
	// find first seed if out of bound or no valid seed boundary
	if ( i1 < i1min || i1 > i1maxVal
		|| i2 < i2min || i2 > i2maxVal
		|| E_isINF(getSeedE(i1,i2)) )
	{
		// first potential seed position
		i1 = i1min;
		i2 = i2min;
	} else {
		// update to next potential seed position
		if (++i1 >= i1maxVal) {
			i1 = i1min;
			i2++;
		}
	}

	// find next valid seed start within range
	while( i2 < i2maxVal && E_isINF(getSeedE(i1,i2))) {
		// update seed position within range
		if (++i1 == i1maxVal) {
			i1 = i1min;
			i2++;
		}
	}

	// check if we found a valid seed in the range
	if (i1 < i1maxVal && i2< i2maxVal) {
		i1_out = i1;
		i2_out = i2;
		return true;
	}
	// no valid next seed found
	return false;
}

//////////////////////////////////////////////////////////////////////////

} // namespace

#endif /* SEEDHANDLER_H_ */
