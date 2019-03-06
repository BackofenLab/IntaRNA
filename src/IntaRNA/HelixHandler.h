
#ifndef INTARNA_HELIXHANDLER_H_
#define INTARNA_HELIXHANDLER_H_

#include "IntaRNA/InteractionEnergy.h"
#include "IntaRNA/HelixConstraint.h"
#include "IntaRNA/SeedHandler.h"

namespace IntaRNA {

 /**
  * Abstract interface to provide helix interaction information
  * If SeedHandler is provided, also helices containing seeds will be computed
  */
class HelixHandler
{

public:

	/**
	 * Construction
	 * @param energy the energy function to be used for helix prediction
	 * @param helixConstraint the helix constraint to be applied
	 */
	HelixHandler();

	/**
	 * Provides a newly allocated helixHandler according to the user defined
	 * parameters
	 * @param energy the interaction energy handler to be used
	 * @param helixConstraint the helixConstraint to be used
	 * @return the newly allocated HelixHandler object
	 */
	static
	HelixHandler* getHelixHandler( const InteractionEnergy & energy
			, const HelixConstraint & helixConstraint
			, SeedHandler * const seedHandler = NULL
	);

	/**
	 * destruction
	 */
	virtual ~HelixHandler();

	/**
	 * Access to the underlying helix constraint
	 * @return the used helix constraint
	 */
	virtual
	const HelixConstraint&
	getConstraint() const = 0;

	/**
	 * Access to the underlying interaction energy function
	 * @return the used energy function
	 */
	virtual
	const InteractionEnergy&
	getInteractionEnergy() const = 0;


	/**
	 * Computes the helix information for the given interval boundaries
	 * @param i1 the first index of seq1 that might interact
	 * @param j1 the last index of seq1 that might interact
	 * @param i2 the first index of seq2 that might interact
	 * @param j2 the last index of seq2 that might interact
	 * @return the number of potential helix interactions
	 */
	virtual
	size_t
	fillHelix(const size_t i1, const size_t j1, const size_t i2, const size_t j2) = 0;

	/**
	 * Computes the helix information for the given interval boundaries for helices containing a seed
	 * @param i1 the first index of seq1 that might interact
	 * @param j1 the last index of seq1 that might interact
	 * @param i2 the first index of seq2 that might interact
	 * @param j2 the last index of seq2 that might interact
	 * @return the number of potential helix interactions
	 */
	virtual
	size_t
	fillHelixSeed(const size_t i1, const size_t j1, const size_t i2, const size_t j2) = 0;

	/**
	 * Identifies the base pairs of the mfe helix interaction starting at i1,i2
	 * and writes them to the provided container
	 *
	 * NOTE: the left- and right-most base pair is excluded!
	 *
	 * @param interaction the container to add the base pairs too
	 * @param i1 the start of the helix in seq1
	 * @param i2 the start of the helix in seq2
	 */
	virtual
	void
	traceBackHelix( Interaction & interaction, const size_t i1, const size_t i2) = 0;

	/**
	 * Identifies the base pairs of the mfe helix interaction,
	 * containing a seed, starting at i1,i2
	 * and writes them to the provided container
	 *
	 * NOTE: the left- and right-most base pair is excluded!
	 *
	 * @param interaction the container to add the base pairs too
	 * @param i1 the start of the helix in seq1
	 * @param i2 the start of the helix in seq2
	 */
	virtual
	void
	traceBackHelixSeed( Interaction & interaction, const size_t i1, const size_t i2) = 0;

	/**
	 * Access to the mfe of any helix with left-most base pair (i1,i2)
	 * @param i1 the left most interacting base of seq1
	 * @param i2 the left most interacting base of seq2
	 * @return the mfe of any helix starting at (i1,i2) or E_INF if none possible
	 */
	virtual
	E_type
	getHelixE( const size_t i1, const size_t i2 ) const = 0;

	/**
	 * Access to the mfe of any helix with left-most base pair (i1,i2)
	 * @param i1 the left most interacting base of seq1
	 * @param i2 the left most interaction base of seq2
	 * @return the mfe of any helix starting at (i1,i2) or E_INF if none possible
	 */
	virtual
	E_type
	getHelixSeedE( const size_t i1, const size_t i2 ) const = 0;

	/**
	 * Access to the length in seq1 of the mfe helix with left-most base pair (i1,i2)
	 * @param i1 the left most interacting base of seq1
	 * @param i2 the left most interacting base of seq2
	 * @return the length in seq1 of the mfe helix starting at (i1,i2) or 0 if none possible
	 */
	virtual
	size_t
	getHelixLength1( const size_t i1, const size_t i2 ) const = 0;

	/**
	 * Access to the length in seq2 of the mfe helix with left-most base pair (i1,i2)
	 * @param i1 the left most interacting base of seq1
	 * @param i2 the left most interacting base of seq2
	 * @return the length in seq2 of the mfe helix starting at (i1,i2) or 0 if none possible
	 */
	virtual
	size_t
	getHelixLength2( const size_t i1, const size_t i2 ) const = 0;

	/**
	 * Access to the length in seq1 of the mfe helix with left-most base pair (i1,i2) containing a seed
	 * @param i1 the left most interacting base of seq1
	 * @param i2 the left most interacting base of seq2
	 * @return the length in seq1 of the mfe helix starting at (i1,i2) or 0 if none possible
	 */
	virtual
	size_t
	getHelixSeedLength1( const size_t i1, const size_t i2 ) const = 0;

	/**
	 * Access to the length in seq2 of the mfe helix with left-most base pair (i1,i2) containing a seed
	 * @param i1 the left most interacting base of seq1
	 * @param i2 the left most interacting base of seq2
	 * @return the length in seq2 of the mfe helix starting at (i1,i2) or 0 if none possible
	 */
	virtual
	size_t
	getHelixSeedLength2( const size_t i1, const size_t i2 ) const = 0;

	/**
	 * Set the seedHandler used in the helixHandlerSeed computation
 	 */
	virtual
	void
	setSeedHandler( SeedHandler & seedHandler ) = 0;

};


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


inline
HelixHandler::HelixHandler()
{
}

////////////////////////////////////////////////////////////////////////////

inline
HelixHandler::~HelixHandler()
{
}

//////////////////////////////////////////////////////////////////////////

} // namespace
#endif //INTARNA_HELIXHANDLER_H_
