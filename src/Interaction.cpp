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
setSeedRange( const BasePair ij1, const BasePair ij2, const E_type energy )
{
	// store seed information
	if (seedRange == NULL) {
		// create new seed information
		seedRange = new InteractionRange(
				// sequences
				*(s1), *(s2),
				// seed ranges
				IndexRange(ij1.first,ij2.first), IndexRange(ij1.second,ij2.second),
				// hybridization loop energies only
				energy
				);
	} else {
		// overwrite
		assert(s1 == seedRange->s1);
		assert(s2 == seedRange->s2);
		// seed ranges
		seedRange->r1 = IndexRange(ij1.first,ij2.first);
		seedRange->r2 = IndexRange(ij1.second,ij2.second),
		// hybridization loop energies only
		seedRange->energy = energy;
	}

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

	// undo seed information
	CLEANUP(seedRange);

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
