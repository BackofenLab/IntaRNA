/*
 * Interaction.cpp
 *
 *  Created on: 14.07.2014
 *      Author: Mmann
 */

#include "InteractionRange.h"
#include "Interaction.h"
#include "general.h"

#include <algorithm>
#include <stdexcept>


////////////////////////////////////////////////////////////////////////////

Interaction::Interaction( const RnaSequence & s1, const RnaSequence & s2 )
:
	s1(&s1)
	, s2(&s2)
	, basePairs()
	, energy( std::numeric_limits<E_type>::signaling_NaN() )
	, seedRange( NULL )
{
}

////////////////////////////////////////////////////////////////////////////

Interaction::Interaction( const InteractionRange & range )
:
	s1(NULL)
	, s2(NULL)
	, basePairs()
	, energy( std::numeric_limits<E_type>::signaling_NaN() )
	, seedRange( NULL )
{
	// init data
	this->operator =( range );
}

////////////////////////////////////////////////////////////////////////////

Interaction::~Interaction()
{
	if (seedRange != NULL) delete seedRange;
}

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

bool
Interaction::
isEmpty() const
{
	return basePairs.size()==0;
}

////////////////////////////////////////////////////////////////////////////

void
Interaction::
sort()
{
	// sort based on STL functionalities for vector and pair
	std::sort(basePairs.begin(), basePairs.end());
}

////////////////////////////////////////////////////////////////////////////

void
Interaction::
clear()
{
	// clear interaction base pairing information
	basePairs.clear();
	// clear energy
	energy = std::numeric_limits<E_type>::signaling_NaN();
	// undo seed information
	if (seedRange != NULL) delete seedRange;

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
	if (seedRange != NULL) delete seedRange;

	return *this;
}



////////////////////////////////////////////////////////////////////////////

bool
Interaction::
compareEnergy( const E_type & energy, const Interaction & hasLargerE )
{
	return energy < hasLargerE.energy && !(E_equal(energy,hasLargerE.energy));
}

////////////////////////////////////////////////////////////////////////////

std::ostream&
operator<<(std::ostream& out, const Interaction::BasePair& bp)
{
	out <<"("<<bp.first<<"-"<<bp.second<<")";
	return out;
}

////////////////////////////////////////////////////////////////////////////

std::ostream&
operator<<(std::ostream& out, const Interaction& i)
{
	for (int p=0; p<i.basePairs.size(); p++) {
		out <<(p==0?"":",") <<i.basePairs.at(p);
	}
	return out;
}

////////////////////////////////////////////////////////////////////////////
