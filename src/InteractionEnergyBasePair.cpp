
#include "InteractionEnergyBasePair.h"


////////////////////////////////////////////////////////////////////////////

InteractionEnergyBasePair::InteractionEnergyBasePair(
		const Accessibility & accS1
		, const ReverseAccessibility & accS2
		, const size_t maxInternalLoopSize1
		, const size_t maxInternalLoopSize2
	)
 :
	InteractionEnergy(accS1, accS2, maxInternalLoopSize1, maxInternalLoopSize2)
{
}

////////////////////////////////////////////////////////////////////////////

InteractionEnergyBasePair::~InteractionEnergyBasePair()
{
}


////////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergyBasePair::
getInterLoopE( const size_t i1, const size_t j1, const size_t i2, const size_t j2 ) const
{
	// if valid internal loop
	if ( isValidInternalLoop(i1,j1,i2,j2) ) {
		// return negated number of gained base pairs by closing this loop = -1
		return getBestStackingEnergy();
	} else {
		return E_INF;
	}
}

////////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergyBasePair::
getDanglingLeft( const size_t i1, const size_t i2 ) const
{
	// no dangling end contribution
	return (E_type)0;
}

////////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergyBasePair::
getDanglingRight( const size_t j1, const size_t j2 ) const
{
	// no dangling end contribution
	return (E_type)0;
}

////////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergyBasePair::
getBestStackingEnergy() const
{
	return getBestInitEnergy();
}

////////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergyBasePair::
getBestInitEnergy() const
{
	return (E_type)-1.0;
}

////////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergyBasePair::
getBestDangleEnergy() const
{
	return (E_type)0.0;
}

////////////////////////////////////////////////////////////////////////////

