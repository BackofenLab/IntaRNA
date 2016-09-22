/*
 * EnergyBasePair.cpp
 *
 *  Created on: 27.06.2014
 *      Author: Mmann
 */

#include "EnergyBasePair.h"

extern "C" {
	#include <ViennaRNA/utils.h>
	#include <ViennaRNA/fold_vars.h>
	#include <ViennaRNA/loop_energies.h>
}


////////////////////////////////////////////////////////////////////////////

EnergyBasePair::EnergyBasePair(
		const Accessibility & accS1
		, const ReverseAccessibility & accS2
		, const size_t maxInternalLoopSize1
		, const size_t maxInternalLoopSize2
	)
 :
	Energy(accS1, accS2, maxInternalLoopSize1, maxInternalLoopSize2)
{
}

////////////////////////////////////////////////////////////////////////////

EnergyBasePair::~EnergyBasePair()
{
}


////////////////////////////////////////////////////////////////////////////

E_type
EnergyBasePair::
getInterLoopE( const size_t i1, const size_t j1, const size_t i2, const size_t j2 )
{
	// if region indices and loop sizes are ok
	if (Energy::isAllowedLoopRegion(accS1.getSequence(), i1, j1, maxInternalLoopSize1)
		&& Energy::isAllowedLoopRegion(accS2.getSequence(), i2, j2, maxInternalLoopSize2)
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
EnergyBasePair::
getDanglingLeft( const size_t i1, const size_t i2 )
{
	// no dangling end contribution
	return (E_type)0;
}

////////////////////////////////////////////////////////////////////////////

E_type
EnergyBasePair::
getDanglingRight( const size_t j1, const size_t j2 )
{
	// no dangling end contribution
	return (E_type)0;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

