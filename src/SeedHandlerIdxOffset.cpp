
#include "SeedHandlerIdxOffset.h"

////////////////////////////////////////////////////////////////////////////

SeedHandlerIdxOffset::SeedHandlerIdxOffset(
		const InteractionEnergy & energy
		, const SeedConstraint & seedConstraint
		)
	:
		seedHandlerOriginal( energy, seedConstraint )
		, seedConstraintOffset(seedConstraint)
		, idxOffset1(0)
		, idxOffset2(0)
{
}

////////////////////////////////////////////////////////////////////////////

SeedHandlerIdxOffset::~SeedHandlerIdxOffset()
{
}

////////////////////////////////////////////////////////////////////////////

const SeedConstraint&
SeedHandlerIdxOffset::
getConstraint() const
{
	return seedConstraintOffset;
}

//////////////////////////////////////////////////////////////////////////

void
SeedHandlerIdxOffset::
fillSeed( const size_t i1min, const size_t i1max, const size_t i2min, const size_t i2max)
{
	seedHandlerOriginal.fillSeed( i1min+idxOffset1, i1max+idxOffset1, i2min+idxOffset2, i2max+idxOffset2 );
}

//////////////////////////////////////////////////////////////////////////

void
SeedHandlerIdxOffset::
traceBackSeed( Interaction & interaction
		, const size_t i1
		, const size_t i2
		)
{
	seedHandlerOriginal.traceBackSeed( interaction, i1+idxOffset1, i2+idxOffset2 );
}

//////////////////////////////////////////////////////////////////////////

E_type
SeedHandlerIdxOffset::
getSeedE( const size_t i1, const size_t i2 ) const
{
	return seedHandlerOriginal.getSeedE( i1+idxOffset1, i2+idxOffset2 );
}

//////////////////////////////////////////////////////////////////////////

size_t
SeedHandlerIdxOffset::
getSeedLength1( const size_t i1, const size_t i2 ) const
{
	return seedHandlerOriginal.getSeedLength1( i1+idxOffset1, i2+idxOffset2 );
}

//////////////////////////////////////////////////////////////////////////

size_t
SeedHandlerIdxOffset::
getSeedLength2( const size_t i1, const size_t i2 ) const
{
	return seedHandlerOriginal.getSeedLength2( i1+idxOffset1, i2+idxOffset2 );
}

//////////////////////////////////////////////////////////////////////////

size_t
SeedHandlerIdxOffset::
getOffset1() const
{
	return idxOffset1;
}

//////////////////////////////////////////////////////////////////////////

size_t
SeedHandlerIdxOffset::
getOffset2() const
{
	return idxOffset2;
}

//////////////////////////////////////////////////////////////////////////

void
SeedHandlerIdxOffset::
setOffset1( const size_t offset )
{
#if IN_DEBUG_MODE
	if (offset >= seedHandlerOriginal.getInteractionEnergy().size1()) {
		throw std::runtime_error("SeedHandlerIdxOffset.setOffset1("+toString(offset)
				+") offset > seq1.length "+toString(seedHandlerOriginal.getInteractionEnergy().size1()));
	}
#endif
	// set idx offset
	this->idxOffset1 = offset;
	// update ranges of seed constraint
	seedConstraintOffset.getRanges1() = seedHandlerOriginal.getConstraint().getRanges1().shift( -(int)offset, seedHandlerOriginal.getInteractionEnergy().size1()-1-offset );
}

//////////////////////////////////////////////////////////////////////////

void
SeedHandlerIdxOffset::
setOffset2( const size_t offset )
{
#if IN_DEBUG_MODE
	if (offset >= seedHandlerOriginal.getInteractionEnergy().size2()) {
		throw std::runtime_error("SeedHandlerIdxOffset.setOffset2("+toString(offset)
				+") offset > seq2.length "+toString(seedHandlerOriginal.getInteractionEnergy().size2()));
	}
#endif
	// set idx offset
	this->idxOffset2 = offset;
	// update ranges of seed constraint
	seedConstraintOffset.getRanges2() = seedHandlerOriginal.getConstraint().getRanges2().shift( -(int)offset, seedHandlerOriginal.getInteractionEnergy().size2()-1-offset );
}

////////////////////////////////////////////////////////////////////////////
