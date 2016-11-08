
#include "InteractionEnergyIdxOffset.h"

#include <boost/foreach.hpp>

//////////////////////////////////////////////////////////////////////////

InteractionEnergyIdxOffset::
InteractionEnergyIdxOffset( const InteractionEnergy & energyOriginal
			, const size_t offset1
			, const size_t offset2 )

 :
	InteractionEnergy( energyOriginal )
	, energyOriginal(energyOriginal)
	, offset1(offset1)
	, offset2(offset2)
{
#if IN_DEBUG_MODE
	// input sanity checks
	setOffset1(offset1);
	setOffset2(offset2);
#endif
}

//////////////////////////////////////////////////////////////////////////

InteractionEnergyIdxOffset::~InteractionEnergyIdxOffset()
{
}

//////////////////////////////////////////////////////////////////////////

size_t
InteractionEnergyIdxOffset::
getOffset1() const
{
	return offset1;
}

//////////////////////////////////////////////////////////////////////////

void
InteractionEnergyIdxOffset::
setOffset1(size_t offset1)
{
#if IN_DEBUG_MODE
	if (offset1 >= energyOriginal.size1()) {
		throw std::runtime_error("InteractionEnergyIdxOffset : offset1 "+toString(offset1)
				+" > seq1.length "+toString(energyOriginal.size1()));
	}
#endif
	this->offset1 = offset1;
}

//////////////////////////////////////////////////////////////////////////

size_t
InteractionEnergyIdxOffset::
getOffset2() const
{
	return offset2;
}

//////////////////////////////////////////////////////////////////////////

void
InteractionEnergyIdxOffset::
setOffset2(size_t offset2)
{
#if IN_DEBUG_MODE
	if (offset2 >= energyOriginal.size2()) {
		throw std::runtime_error("InteractionEnergyIdxOffset : offset2 "+toString(offset2)
				+" > seq2.length "+toString(energyOriginal.size2()));
	}
#endif
	this->offset2 = offset2;
}

//////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergyIdxOffset::
getED1( const size_t i1, const size_t j1 ) const
{
	return energyOriginal.getED1( i1+offset1, j1+offset1 );
}

//////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergyIdxOffset::
getED2( const size_t i2, const size_t j2 ) const
{
	return energyOriginal.getED2( i2+offset2, j2+offset2 );
}

////////////////////////////////////////////////////////////////////////////

bool
InteractionEnergyIdxOffset::
isAccessible1( const size_t i ) const
{
	return energyOriginal.isAccessible1(i+offset1);
}

////////////////////////////////////////////////////////////////////////////

bool
InteractionEnergyIdxOffset::
isAccessible2( const size_t i ) const
{
	return energyOriginal.isAccessible2(i+offset2);
}

//////////////////////////////////////////////////////////////////////////

bool
InteractionEnergyIdxOffset::
areComplementary( const size_t i1, const size_t i2 ) const
{
	return energyOriginal.areComplementary( i1+offset1, i2+offset2 );
}

//////////////////////////////////////////////////////////////////////////

size_t
InteractionEnergyIdxOffset::
size1() const
{
	return energyOriginal.size1()-offset1;
}

////////////////////////////////////////////////////////////////////////////

size_t
InteractionEnergyIdxOffset::
size2() const
{
	return energyOriginal.size2()-offset2;
}


E_type
InteractionEnergyIdxOffset::
getE_init( ) const
{
	return energyOriginal.getE_init();
}

//////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergyIdxOffset::
getE_interLeft( const size_t i1, const size_t j1, const size_t i2, const size_t j2 ) const
{
	return energyOriginal.getE_interLeft(i1+offset1, j1+offset1, i2+offset2, j2+offset2);
}


//////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergyIdxOffset::
getE_danglingLeft( const size_t i1, const size_t i2 ) const
{
	return energyOriginal.getE_danglingLeft( i1+offset1, i2+offset2 );
}

//////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergyIdxOffset::
getE_danglingRight( const size_t j1, const size_t j2 ) const
{
	return energyOriginal.getE_danglingRight( j1+offset1, j2+offset2 );
}

//////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergyIdxOffset::
getE_endLeft( const size_t i1, const size_t i2 ) const
{
	return energyOriginal.getE_endLeft( i1+offset1, i2+offset2 );
}

//////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergyIdxOffset::
getE_endRight( const size_t j1, const size_t j2 ) const
{
	return energyOriginal.getE_endRight( j1+offset1, j2+offset2 );
}

//////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergyIdxOffset::
getPr_danglingLeft( const size_t i1, const size_t j1, const size_t i2, const size_t j2 ) const
{
	return energyOriginal.getPr_danglingLeft(i1+offset1, j1+offset1, i2+offset2, j2+offset2);
}

//////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergyIdxOffset::
getPr_danglingRight( const size_t i1, const size_t j1, const size_t i2, const size_t j2 ) const
{
	return energyOriginal.getPr_danglingRight(i1+offset1, j1+offset1, i2+offset2, j2+offset2);
}

//////////////////////////////////////////////////////////////////////////

Interaction::BasePair
InteractionEnergyIdxOffset::
getBasePair( const size_t i1, const size_t i2 ) const
{
	return energyOriginal.getBasePair( i1+offset1, i2+offset2 );
}

//////////////////////////////////////////////////////////////////////////

size_t
InteractionEnergyIdxOffset::
getIndex1( const Interaction::BasePair & bp ) const
{
#if IN_DEBUG_MODE
	if (energyOriginal.getIndex1(bp)<offset1) throw std::runtime_error("InteractionEnergyIdxOffset::getIndex1("+toString(energyOriginal.getIndex1(bp))+") < offset1 = "+toString(offset1));
#endif
	return energyOriginal.getIndex1(bp)-offset1;
}

////////////////////////////////////////////////////////////////////////////

size_t
InteractionEnergyIdxOffset::
getIndex2( const Interaction::BasePair & bp ) const
{
#if IN_DEBUG_MODE
	if (energyOriginal.getIndex2(bp)<offset2) throw std::runtime_error("InteractionEnergyIdxOffset::getIndex2("+toString(energyOriginal.getIndex2(bp))+") < offset2 = "+toString(offset2));
#endif
	return energyOriginal.getIndex2(bp)-offset2;
}

////////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergyIdxOffset::
getRT() const
{
	return energyOriginal.getRT();
}
////////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergyIdxOffset::
getBestE_interLoop() const
{
	return energyOriginal.getBestE_interLoop();
}

////////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergyIdxOffset::
getBestE_dangling() const
{
	return energyOriginal.getBestE_dangling();
}


////////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergyIdxOffset::
getBestE_end() const
{
	return energyOriginal.getBestE_end();
}

