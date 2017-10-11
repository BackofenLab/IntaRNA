
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
	 * NOTE: the right most base pair is excluded!
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
