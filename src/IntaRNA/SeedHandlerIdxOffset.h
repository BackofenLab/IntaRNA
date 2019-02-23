
#ifndef INTARNA_SEEDHANDLERIDXOFFSET_H_
#define INTARNA_SEEDHANDLERIDXOFFSET_H_

#include "IntaRNA/SeedHandler.h"
#include "IntaRNA/InteractionEnergyIdxOffset.h"

#include <boost/multi_array.hpp>

#include <boost/numeric/ublas/matrix.hpp>


namespace IntaRNA {

/**
 * Wrapper around a SeedHandler object where indices are shifted by
 * a given positive offset (shifted towards infinity).
 * This is useful for local interaction computations.
 *
 * @author Martin Mann
 *
 */
class SeedHandlerIdxOffset : public SeedHandler
{

public:

	/**
	 * Construction
	 * @param energy the energy function to be used for seed prediction (already offset)
	 * @param seedHandler the seed handler to be wrapped
	 */
	SeedHandlerIdxOffset( SeedHandler * seedHandler );

	/**
	 * destruction of this object and the wrapped seed handler
	 */
	virtual ~SeedHandlerIdxOffset();

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
	 *
	 * @return number of possible seed interactions
	 */
	virtual
	size_t
	fillSeed(const size_t i1, const size_t j1, const size_t i2, const size_t j2);

	/**
	 * Identifies the base pairs of the mfe seed interaction starting at i1,i2
	 * and writes them to the provided container
	 *
	 * NOTE: the left- and right-most base pairs are excluded!
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


	/**
	 * Access to the wrapped SeedHandler instance.
	 *
	 * @return the internally used SeedHandler without offset
	 */
	virtual
	SeedHandler&
	getOriginalSeedHandler();


protected:

	//! the index shifted seed handler
	SeedHandler * seedHandlerOriginal;

	//! the index shifted seed constraint
	SeedConstraint seedConstraintOffset;

	//! offset for indices in seq1
	size_t idxOffset1;

	//! offset for indices in seq2
	size_t idxOffset2;

	//! dedicated energy object to avoid
	InteractionEnergyIdxOffset energyIdxOffset;

};




////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


inline
SeedHandlerIdxOffset::
SeedHandlerIdxOffset( SeedHandler * seedHandlerInstance )
	:
		SeedHandler(seedHandlerInstance->getInteractionEnergy(), seedHandlerInstance->getConstraint() )
		, seedHandlerOriginal( seedHandlerInstance )
		, seedConstraintOffset( seedHandlerOriginal->getConstraint() )
		, idxOffset1(0)
		, idxOffset2(0)
		, energyIdxOffset(seedHandlerInstance->getInteractionEnergy(),idxOffset1,idxOffset2)
{
}

////////////////////////////////////////////////////////////////////////////

inline
SeedHandlerIdxOffset::~SeedHandlerIdxOffset()
{
	if (seedHandlerOriginal != NULL) {
		delete seedHandlerOriginal;
		seedHandlerOriginal = NULL;
	}
}

////////////////////////////////////////////////////////////////////////////

inline
const SeedConstraint&
SeedHandlerIdxOffset::
getConstraint() const
{
	return seedConstraintOffset;
}

////////////////////////////////////////////////////////////////////////////

inline
const InteractionEnergy&
SeedHandlerIdxOffset::
getInteractionEnergy() const
{
	return this->energyIdxOffset;
}

//////////////////////////////////////////////////////////////////////////

inline
size_t
SeedHandlerIdxOffset::
fillSeed( const size_t i1min, const size_t i1max, const size_t i2min, const size_t i2max)
{
	return seedHandlerOriginal->fillSeed( i1min+idxOffset1, i1max+idxOffset1, i2min+idxOffset2, i2max+idxOffset2 );
}

//////////////////////////////////////////////////////////////////////////

inline
void
SeedHandlerIdxOffset::
traceBackSeed( Interaction & interaction
		, const size_t i1
		, const size_t i2
		)
{
	seedHandlerOriginal->traceBackSeed( interaction, i1+idxOffset1, i2+idxOffset2 );
}

//////////////////////////////////////////////////////////////////////////

inline
E_type
SeedHandlerIdxOffset::
getSeedE( const size_t i1, const size_t i2 ) const
{
	return seedHandlerOriginal->getSeedE( i1+idxOffset1, i2+idxOffset2 );
}

//////////////////////////////////////////////////////////////////////////

inline
size_t
SeedHandlerIdxOffset::
getSeedLength1( const size_t i1, const size_t i2 ) const
{
	return seedHandlerOriginal->getSeedLength1( i1+idxOffset1, i2+idxOffset2 );
}

//////////////////////////////////////////////////////////////////////////

inline
size_t
SeedHandlerIdxOffset::
getSeedLength2( const size_t i1, const size_t i2 ) const
{
	return seedHandlerOriginal->getSeedLength2( i1+idxOffset1, i2+idxOffset2 );
}

//////////////////////////////////////////////////////////////////////////

inline
size_t
SeedHandlerIdxOffset::
getOffset1() const
{
	return idxOffset1;
}

//////////////////////////////////////////////////////////////////////////

inline
size_t
SeedHandlerIdxOffset::
getOffset2() const
{
	return idxOffset2;
}

//////////////////////////////////////////////////////////////////////////

inline
void
SeedHandlerIdxOffset::
setOffset1( const size_t offset )
{
#if INTARNA_IN_DEBUG_MODE
	if (offset >= seedHandlerOriginal->getInteractionEnergy().size1()) {
		throw std::runtime_error("SeedHandlerIdxOffset.setOffset1("+toString(offset)
				+") offset > seq1.length "+toString(seedHandlerOriginal->getInteractionEnergy().size1()));
	}
#endif
	// set idx offset
	this->idxOffset1 = offset;
	this->energyIdxOffset.setOffset1(idxOffset1);
	// update ranges of seed constraint
	seedConstraintOffset.getRanges1() = seedHandlerOriginal->getConstraint().getRanges1().shift( -(int)offset, seedHandlerOriginal->getInteractionEnergy().size1()-1-offset );
}

//////////////////////////////////////////////////////////////////////////

inline
void
SeedHandlerIdxOffset::
setOffset2( const size_t offset )
{
#if INTARNA_IN_DEBUG_MODE
	if (offset >= seedHandlerOriginal->getInteractionEnergy().size2()) {
		throw std::runtime_error("SeedHandlerIdxOffset.setOffset2("+toString(offset)
				+") offset > seq2.length "+toString(seedHandlerOriginal->getInteractionEnergy().size2()));
	}
#endif
	// set idx offset
	this->idxOffset2 = offset;
	this->energyIdxOffset.setOffset2(idxOffset2);
	// update ranges of seed constraint
	seedConstraintOffset.getRanges2() = seedHandlerOriginal->getConstraint().getRanges2().shift( -(int)offset, seedHandlerOriginal->getInteractionEnergy().size2()-1-offset );
}

////////////////////////////////////////////////////////////////////////////

inline
SeedHandler&
SeedHandlerIdxOffset::
getOriginalSeedHandler()
{
	return *seedHandlerOriginal;
}

////////////////////////////////////////////////////////////////////////////

} // namespace

#endif /* SEEDHANDLERIDXOFFSET_H_ */
