/*
 * Interaction.cpp
 *
 *  Created on: 14.07.2014
 *      Author: Mmann
 */

#include "Interaction.h"
#include "general.h"

#include <algorithm>
#include <stdexcept>


////////////////////////////////////////////////////////////////////////////

bool
Interaction::
isValid() const
{

	// no or single interaction
	if (basePairs.size() < 2) {
		// check if at least one interacting base pair
		return !isEmpty();
	}

	// multiple interacting base pairs
	PairingVec::const_iterator i = basePairs.begin(), j = basePairs.begin();
	// index order and duplicate check
	bool isValid = true;
	for (++j; isValid && j!=basePairs.end(); ++i,++j ) {
		isValid = (i->first < j->first) && (i->second > j->second);
	}

	return isValid;

}

////////////////////////////////////////////////////////////////////////////

void
Interaction::
setSeedRange( const BasePair bp_left, const BasePair bp_right, const E_type energy )
{
	// check if container available
	if (seed == NULL) {
		// create new seed information
		seed = new Seed();
	}
	// set seed data
	seed->bp_i = bp_left;
	seed->bp_j = bp_right;
	seed->energy = energy;

}

////////////////////////////////////////////////////////////////////////////

Interaction &
Interaction::
operator= ( const Interaction & toCopy )
{
#if IN_DEBUG_MODE
	if (!toCopy.isValid())
		throw std::runtime_error("Interaction::=("+toString(toCopy)+") not valid!");
#endif
	// clear current interactions
	basePairs.clear();


	// copy sequence handles
	s1 = toCopy.s1;
	s2 = toCopy.s2;

	// copy base pair
	basePairs = toCopy.basePairs;

	// copy energy value
	energy = toCopy.energy;

	// copy seed data
	if (toCopy.seed != NULL) {
		// create seed info if not existing
		if (seed == NULL) { seed = new Seed(); }
		// copy data
		*seed = *(toCopy.seed);
	} else {
		// remove seed information if present
		CLEANUP(seed);
	}

	return *this;
}

////////////////////////////////////////////////////////////////////////////

Interaction &
Interaction::
operator= ( const InteractionRange & range )
{
#if IN_DEBUG_MODE
	if (!range.isSane())
		throw std::runtime_error("Interaction::=("+toString(range)+") not sane!");
#endif
	// clear current interactions
	basePairs.clear();

	// undo seed information
	CLEANUP(seed);

	// copy sequence handles
	s1 = range.s1;
	s2 = range.s2;

	// add left boundary base pair
	basePairs.push_back( BasePair(range.r1.from, range.r2.from) );
	// add right boundary base pair if both not singleton ranges
	if ( range.r1.from != range.r1.to || range.r2.from != range.r2.to ) {
		basePairs.push_back( BasePair(range.r1.to, range.r2.to) );
	}

	// copy energy value
	energy = range.energy;

	return *this;
}

////////////////////////////////////////////////////////////////////////////

std::string
Interaction::
dotBar( const Interaction & i )
{
#if IN_DEBUG_MODE
	if (!i.isValid())
		throw std::runtime_error("Interaction::dotBar("+toString(i)+") not valid!");
#endif
	// compile overall dot-bracket representation
	return	toString(i.basePairs.begin()->first +1)
			+ dotSomething(i.basePairs.begin(), i.basePairs.end(), true, '|')
			+"&"
			+toString(i.basePairs.rbegin()->second +1)
			+ dotSomething(i.basePairs.rbegin(), i.basePairs.rend(), false, '|')
			;
}

////////////////////////////////////////////////////////////////////////////

std::string
Interaction::
dotBracket( const Interaction & i, const char symOpen, const char symClose )
{
#if IN_DEBUG_MODE
	if (!i.isValid())
		throw std::runtime_error("Interaction::dotBracket("+toString(i)+") not valid!");
#endif
	// compile overall dot-bracket representation
	return	dotSomething(i.basePairs.begin(), i.basePairs.end(), true, symOpen)
			+"&"
			+ dotSomething(i.basePairs.rbegin(), i.basePairs.rend(), false, symClose)
			;
}

////////////////////////////////////////////////////////////////////////////
