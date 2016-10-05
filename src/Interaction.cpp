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
{
}

////////////////////////////////////////////////////////////////////////////

Interaction::Interaction( const InteractionRange & range )
:
	s1(NULL)
	, s2(NULL)
	, basePairs()
	, energy( std::numeric_limits<E_type>::signaling_NaN() )
{
	// init data
	this->operator =( range );
}

////////////////////////////////////////////////////////////////////////////

Interaction::~Interaction()
{
}

////////////////////////////////////////////////////////////////////////////

void
Interaction::
addInteraction( const size_t i1, const size_t i2 )
{
#ifdef NDEBUG           /* required by ANSI standard */
	// no check
#else
	// check if sane indices
	if (i1 >= s1->size()) {
		throw std::runtime_error("Interaction::addInteraction: index i1="+toString(i1)+" exceed first sequence's length ("+toString(s1->size())+")");
	}
	if (i2 >= s2->size()) {
		throw std::runtime_error("Interaction::addInteraction: index i2="+toString(i2)+" exceed second sequence's length ("+toString(s2->size())+")");
	}
	if ( ! RnaSequence::areComplementary( *s1, *s2, i1, i2 ) ) {
		throw std::runtime_error("Interaction::addInteraction: positions "+toString(i1)+" and "+toString(i2)+" are not complementary");
	}
#endif

	// push to according interaction container
	basePairs.push_back( BasePair(i1,i2) );
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
}


////////////////////////////////////////////////////////////////////////////

Interaction &
Interaction::
operator= ( const InteractionRange & range )
{
#ifdef NDEBUG           /* required by ANSI standard */
	// no check
#else
	if (!range.isSane())
		throw std::runtime_error("Interaction::=("+toString(range)+") not sane!");
#endif
	// clear current interactions
	basePairs.clear();

	// copy sequence handles
	s1 = range.s1;
	s2 = range.s2;

	// add left boundary base pair
	addInteraction( range.r1.from, range.r2.from );
	// add right boundary base pair if both not singleton ranges
	if ( range.r1.from != range.r1.to || range.r2.from != range.r2.to ) {
		addInteraction( range.r1.to, range.r2.to );
	}

	// copy energy value
	energy = range.energy;

	return *this;
}


////////////////////////////////////////////////////////////////////////////
