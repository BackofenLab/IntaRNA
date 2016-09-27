
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
getInterLoopE( const size_t i1, const size_t j1, const size_t i2, const size_t j2 )
{
	// if region indices and loop sizes are ok
	if (InteractionEnergy::isAllowedLoopRegion(accS1.getSequence(), i1, j1, maxInternalLoopSize1)
		&& InteractionEnergy::isAllowedLoopRegion(accS2.getSequence(), i2, j2, maxInternalLoopSize2)
		// TODO check if interaction base pair is possible
	) {
		// return negated number of gained base pairs by closing this loop = -1
		return (E_type)-1.0;
	} else {
		return E_MAX;
	}
}

////////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergyBasePair::
getDanglingLeft( const size_t i1, const size_t i2 )
{
	// no dangling end contribution
	return (E_type)0;
}

////////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergyBasePair::
getDanglingRight( const size_t j1, const size_t j2 )
{
	// no dangling end contribution
	return (E_type)0;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

