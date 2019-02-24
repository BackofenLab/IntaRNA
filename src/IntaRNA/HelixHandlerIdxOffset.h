
#ifndef INTARNA_HELIXHANDLERIDXOFFSET_H_
#define INTARNA_HELIXHANDLERIDXOFFSET_H_

#include "IntaRNA/HelixHandler.h"
#include "IntaRNA/SeedHandlerIdxOffset.h"

#include <boost/multi_array.hpp>

#include <boost/numeric/ublas/matrix.hpp>

namespace IntaRNA {

/**
 * Wrapper around a HelixHandler object where indices are shifted by
 * a given positive offset (shifted towards infinity).
 * This is useful for local interaction computations.
 *
 * @author Rick Gelhausen
 */
class HelixHandlerIdxOffset : public HelixHandler
{

public:

	/**
	 * Construction
	 *
	 */
	HelixHandlerIdxOffset( HelixHandler * helixHandler );

	/**
	 * destruction
	 */
	virtual ~HelixHandlerIdxOffset();

	/**
	 * Access to the underlying interaction energy function
	 * @return the used energy function
	 */
	virtual
	const InteractionEnergy&
	getInteractionEnergy() const;

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
	 * @param offset2 the index offset for sequence 2 to be used
	 */
	void setOffset2(size_t offset2);

	///////////////  OVERWRITTEN MEMBERS USING OFFSET  /////////////////

	/**
	 * Access to the underlying helix constraint
	 * @return the used helix constraint
	 */
	virtual
	const HelixConstraint&
	getConstraint() const;

	/**
	 * Computes the helix matrix for the given interval boundaries
	 *
	 * Note, the indices are shifted by an offset for computation.
	 *
	 * @param i1 the first index of seq1 that might interact
	 * @param j1 the last index of seq1 that might interact
	 * @param i2 the first index of seq2 that might interact
	 * @param j2 the last index of seq2 that might interact
	 *
	 * @return number of possible helix interactions
	 */
	virtual
	size_t
	fillHelix(const size_t i1, const size_t j1, const size_t i2, const size_t j2);

	/**
	 * Computes the helix matrix for the given interval boundaries for helices containing a seed
	 *
	 * Note, the indices are shifted by an offset for computation.
	 *
	 * @param i1 the first index of seq1 that might interact
	 * @param j1 the last index of seq1 that might interact
	 * @param i2 the first index of seq2 that might interact
	 * @param j2 the last index of seq2 that might interact
	 * @return
	 */
	virtual
	size_t
	fillHelixSeed(const size_t i1, const size_t j1, const size_t i2, const size_t j2);

	/**
	 * Identifies the base pairs of the mfe helix interaction starting at i1,i2
	 * and writes them to the provided container
	 *
	 * NOTE: the right most base pair is excluded!
	 *
	 * Note, the indices are shifted by an offset for computation.
	 *
	 * @param interaction the container to add the base pairs too
	 * @param i1 the start of the helix in seq1
	 * @param i2 the start of the helix in seq2
	 */
	virtual
	void
	traceBackHelix( Interaction & interaction, const size_t i1, const size_t i2);

	/**
	 * Identifies the base pairs of the mfe helix interaction, containing a seed, starting at i1,i2
	 * and writes them to the provided container
	 *
	 * Note, the indices are shifted by an offset for computation.
	 *
	 * @param interaction the container to add the base pairs too
	 * @param i1 the start of the helix in seq1
	 * @param i2 the start of the helix in seq2
	 */
	virtual
	void
	traceBackHelixSeed( Interaction & interaction, const size_t i1, const size_t i2);

	/**
	 * Access to the mfe of any helix with left-most base pair (i1,i2)
	 *
	 * Note, the indices are shifted by an offset for computation
	 *
	 * @param i1 the left most interacting base of seq1
	 * @param i2 the left most interacting base of seq2
	 * @return the mfe of any helix starting at (i1,i2) or E_INF if none possible
	 */
	virtual
	E_type
	getHelixE( const size_t i1, const size_t i2 ) const;

	/**
	 * Access to the mfe of any helix with left-most base pair (i1,i2) and containing a seed
	 *
	 * Note, the indices are shifted by an offset for computation
	 *
	 * @param i1 the left most interacting base of seq1
	 * @param i2 the left most interacting base of seq2
	 * @return the mfe of any helix starting at (i1,i2) or E_INF if none possible
	 */
	virtual
	E_type
	getHelixSeedE( const size_t i1, const size_t i2 ) const;

	/**
	 * Access to the length in seq1 of the mfe helix with left-most base pair (i1,i2)
	 *
	 * Note, the indices are shifted by an offset for computation.
	 *
	 * @param i1 the left most interaction base of seq1
	 * @param i2 the left most interaction base of seq2
	 * @return the length in seq1 of the mfe helix starting at (i1,i2) or 0 if none possible
	 */
	virtual
	size_t
	getHelixLength1( const size_t i1, const size_t i2 ) const;


	/**
	 * Access to the length in seq2 of the mfe helix with left-most base pair (i1,i2)
	 *
	 * Note, the indices are shifted by an offset for computation.
	 *
	 * @param i1 the left most interaction base of seq1
	 * @param i2 the left most interaction base of seq2
	 * @return the length in seq2 of the mfe helix starting at (i1,i2) or 0 if none possible
	 */
	virtual
	size_t
	getHelixLength2( const size_t i1, const size_t i2 ) const;

	/**
	 * Access to the length in seq1 of the mfe helix with left-most base pair (i1,i2), containing a seed
	 *
	 * Note, the indices are shifted by an offset for computation.
	 *
	 * @param i1 the left most interaction base of seq1
	 * @param i2 the left most interaction base of seq2
	 * @return the length in seq1 of the mfe helix starting at (i1,i2) or 0 if none possible
	 */
	virtual
	size_t
	getHelixSeedLength1( const size_t i1, const size_t i2 ) const;

	/**
	 * Access to the length in seq2 of the mfe helix with left-most base pair (i1,i2), containing a seed
	 *
	 * Note, the indices are shifted by an offset for computation.
	 *
	 * @param i1 the left most interaction base of seq1
	 * @param i2 the left most interaction base of seq2
	 * @return the length in seq2 of the mfe helix starting at (i1,i2) or 0 if none possible
	 */
	virtual
	size_t
	getHelixSeedLength2( const size_t i1, const size_t i2 ) const;

	/**
	 * Check whether the given seedHandler is of type SeedHandlerIdxOffset
	 * Set the seedHandler of helixHandlerOriginal
	 *
	 * @param seedHandler the seedHandler with offset
	 */
	virtual
	void
	setSeedHandler( SeedHandler & seedHandler);

protected:

	//! the index shifted helixHandler
	HelixHandler * helixHandlerOriginal;

	//! the index shifted helix constraint
	HelixConstraint helixConstraintOffset;

	//! the seedHandler with offset
	//! Note: (might be needed later)
	SeedHandlerIdxOffset * seedHandlerIdxOffset;

	//! offset for indices in seq1
	size_t idxOffset1;

	//! offset for indices in seq2
	size_t idxOffset2;

};


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


inline
HelixHandlerIdxOffset::
HelixHandlerIdxOffset( HelixHandler * helixHandlerInstance )
	:
		helixHandlerOriginal( helixHandlerInstance )
		, helixConstraintOffset( helixHandlerOriginal->getConstraint() )
		, idxOffset1(0)
		, idxOffset2(0)
		, seedHandlerIdxOffset(NULL)
{

}

////////////////////////////////////////////////////////////////////////////

inline
HelixHandlerIdxOffset::~HelixHandlerIdxOffset()
{
	if (helixHandlerOriginal != NULL) {
		delete helixHandlerOriginal;
		helixHandlerOriginal = NULL;
	}
}

////////////////////////////////////////////////////////////////////////////

inline
const InteractionEnergy&
HelixHandlerIdxOffset::
getInteractionEnergy() const
{
	return helixHandlerOriginal->getInteractionEnergy();
}

////////////////////////////////////////////////////////////////////////////

inline
const HelixConstraint&
HelixHandlerIdxOffset::
getConstraint() const
{
	return helixConstraintOffset;
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerIdxOffset::
fillHelix(const size_t i1min, const size_t i1max, const size_t i2min, const size_t i2max)
{
	return helixHandlerOriginal->fillHelix( i1min+idxOffset1, i1max+idxOffset1, i2min+idxOffset2, i2max+idxOffset2 );
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerIdxOffset::
fillHelixSeed(const size_t i1min, const size_t i1max, const size_t i2min, const size_t i2max)
{
	return helixHandlerOriginal->fillHelixSeed( i1min+idxOffset1, i1max+idxOffset1, i2min+idxOffset2, i2max+idxOffset2 );
}

////////////////////////////////////////////////////////////////////////////

inline
void
HelixHandlerIdxOffset::
traceBackHelix(Interaction &interaction
		, const size_t i1
		, const size_t i2
		)
{
	helixHandlerOriginal->traceBackHelix( interaction, i1+idxOffset1, i2+idxOffset2 );
}

////////////////////////////////////////////////////////////////////////////

inline
void
HelixHandlerIdxOffset::
traceBackHelixSeed(Interaction &interaction
		, const size_t i1
		, const size_t i2
)
{
	helixHandlerOriginal->traceBackHelixSeed( interaction, i1+idxOffset1, i2+idxOffset2 );
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
HelixHandlerIdxOffset::
getHelixE(const size_t i1, const size_t i2) const
{
	return helixHandlerOriginal->getHelixE( i1+idxOffset1, i2+idxOffset2 );
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
HelixHandlerIdxOffset::
getHelixSeedE(const size_t i1, const size_t i2) const
{
	return helixHandlerOriginal->getHelixSeedE( i1+idxOffset1, i2+idxOffset2 );
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerIdxOffset::
getHelixLength1(const size_t i1, const size_t i2) const
{
	return helixHandlerOriginal->getHelixLength1( i1+idxOffset1, i2+idxOffset2 );
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerIdxOffset::
getHelixLength2(const size_t i1, const size_t i2) const
{
	return helixHandlerOriginal->getHelixLength2( i1+idxOffset1, i2+idxOffset2 );
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerIdxOffset::
getHelixSeedLength1(const size_t i1, const size_t i2) const
{
	return helixHandlerOriginal->getHelixSeedLength1( i1+idxOffset1, i2+idxOffset2 );
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerIdxOffset::
getHelixSeedLength2(const size_t i1, const size_t i2) const
{
	return helixHandlerOriginal->getHelixSeedLength2( i1+idxOffset1, i2+idxOffset2 );
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerIdxOffset::
getOffset1() const
{
	return idxOffset1;
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerIdxOffset::
getOffset2() const
{
	return idxOffset2;
}

////////////////////////////////////////////////////////////////////////////

inline
void
HelixHandlerIdxOffset::
setOffset1( const size_t offset )
{
#if INTARNA_IN_DEBUG_MODE
	if (offset >= helixHandlerOriginal->getInteractionEnergy().size1()) {
		throw std::runtime_error("HelixHandlerIdxOffset.setOffset1("+toString(offset)
		+") offset > seq1.length "+toString(helixHandlerOriginal->getInteractionEnergy().size1()));
	}
#endif
	// set idx offset
	this->idxOffset1 = offset;
}

////////////////////////////////////////////////////////////////////////////

inline
void
HelixHandlerIdxOffset::
setOffset2( const size_t offset )
{
#if INTARNA_IN_DEBUG_MODE
	if (offset >= helixHandlerOriginal->getInteractionEnergy().size2()) {
		throw std::runtime_error("HelixHandlerIdxOffset.setOffset1("+toString(offset)
								 +") offset > seq1.length "+toString(helixHandlerOriginal->getInteractionEnergy().size2()));
	}
#endif
	// set idx offset
	this->idxOffset2 = offset;
}

////////////////////////////////////////////////////////////////////////////

inline
void
HelixHandlerIdxOffset::
setSeedHandler(SeedHandler & seedHandler)
{
	// cast input
	SeedHandlerIdxOffset * shOffset = dynamic_cast<SeedHandlerIdxOffset*>(&seedHandler);

#if INTARNA_IN_DEBUG_MODE
	// sanity check
	if (shOffset == NULL) {
		throw std::runtime_error("HelixHandlerIdxOffset.setSeedHandler(). Given seedHandler is not of type SeedHandlerIdxOffset.");
	}
#endif

	// store locally (might be useful later)
	seedHandlerIdxOffset = shOffset;

	// forward seedHandler to helixHandler
	helixHandlerOriginal->setSeedHandler( shOffset->getOriginalSeedHandler());
}

} // namespace
#endif /* HELIXHANDLERIDXOFFSET_H_ */
