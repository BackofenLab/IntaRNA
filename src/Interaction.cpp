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

Interaction::Interaction( const RnaSequence & s1, const RnaSequence & s2 )
:
	s1(s1)
	, s2(s2)
	, interaction()
	, energy( std::numeric_limits<E_type>::signaling_NaN() )
{
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
	if (i1 >= s1.size()) {
		throw std::runtime_error("Interaction::addInteraction: index i1="+toString(i1)+" exceed first sequence's length ("+toString(s1.size())+")");
	}
	if (i2 >= s2.size()) {
		throw std::runtime_error("Interaction::addInteraction: index i2="+toString(i2)+" exceed second sequence's length ("+toString(s2.size())+")");
	}
	if ( ! RnaSequence::areComplementary( s1, s2, i1, i2 ) ) {
		throw std::runtime_error("Interaction::addInteraction: positions "+toString(i1)+" and "+toString(i2)+" are not complementary");
	}
#endif

	// push to according interaction container
	interaction.push_back( BasePair(i1,i2) );
}

////////////////////////////////////////////////////////////////////////////

bool
Interaction::
isValid() const
{
	bool isValid = true;

	// no or single interaction
	if (interaction.size() < 2) {
		// check if at least one interacting base pair
		return interaction.size()>0;
	}

	// multiple interactions
	PairingVec::const_iterator i = interaction.begin(), j = interaction.begin();
	// index order and duplicate check
	for (++j; isValid && j!=interaction.end(); ++i,++j ) {
		isValid = (i->first < j->first) && (i->second > j->second);
	}

	return isValid;

}

////////////////////////////////////////////////////////////////////////////

void
Interaction::
sort()
{
	// sort based on STL functionalities for vector and pair
	std::sort(interaction.begin(), interaction.end());
}

////////////////////////////////////////////////////////////////////////////

void
Interaction::
clear()
{
	// clear interaction base pairing information
	interaction.clear();
}

////////////////////////////////////////////////////////////////////////////

const Interaction::PairingVec &
Interaction::
getBasePairs() const
{
	// access to interacting base pairs
	return interaction;
}


////////////////////////////////////////////////////////////////////////////

Interaction::PairingVec &
Interaction::
getBasePairs()
{
	// access to interacting base pairs
	return interaction;
}


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
