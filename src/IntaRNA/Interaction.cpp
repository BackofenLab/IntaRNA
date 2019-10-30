/*
 * Interaction.cpp
 *
 *  Created on: 14.07.2014
 *      Author: Mmann
 */

#include "IntaRNA/Interaction.h"
#include "IntaRNA/general.h"

#include <algorithm>
#include <stdexcept>




namespace IntaRNA {

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

Interaction::Seed::
Seed()
	: bp_i(Interaction::BasePair(RnaSequence::lastPos, RnaSequence::lastPos))
	, bp_j(Interaction::BasePair(RnaSequence::lastPos, RnaSequence::lastPos))
	, energy(E_INF)
{

}

////////////////////////////////////////////////////////////////////////////

Interaction::Seed::
Seed( const Interaction::BasePair & bp_i
		, const Interaction::BasePair & bp_j
		, const E_type energy )
	: bp_i(bp_i)
	, bp_j(bp_j)
	, energy(energy)
{

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

void
Interaction::
setSeedRange( const BasePair bp_left, const BasePair bp_right, const E_type energy )
{
	// check if container available
	if (seed == NULL) {
		// create new seed information
		seed = new SeedSet();
	}
	// set seed data
	seed->insert( Seed(bp_left, bp_right, energy) );
}

////////////////////////////////////////////////////////////////////////////

IndexRangeList
Interaction::
getSeedRanges1() const
{
	IndexRangeList ranges( true );

	// copy all seed ranges for seq1
	if (seed != NULL) {
		for ( const auto & s : (*seed) ) {
			ranges.insert( IndexRange(s.bp_i.first,s.bp_j.first) );
		}
	}

	return ranges;
}

////////////////////////////////////////////////////////////////////////////

IndexRangeList
Interaction::
getSeedRanges2() const
{
	IndexRangeList ranges( true );

	// copy all reversed seed ranges for seq2
	if (seed != NULL) {
		for ( const auto & s : (*seed) ) {
			ranges.insert( IndexRange(s.bp_j.second,s.bp_i.second) );
		}
	}

	return ranges;
}

////////////////////////////////////////////////////////////////////////////

Interaction &
Interaction::
operator= ( const Interaction & toCopy )
{
#if INTARNA_IN_DEBUG_MODE
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
		if (seed == NULL) { seed = new SeedSet(); }
		else { seed->clear(); }
		// copy data
		seed->insert(toCopy.seed->begin(), toCopy.seed->end());
	} else {
		// remove seed information if present
		INTARNA_CLEANUP(seed);
	}

	return *this;
}

////////////////////////////////////////////////////////////////////////////

Interaction &
Interaction::
operator= ( const InteractionRange & range )
{
#if INTARNA_IN_DEBUG_MODE
	if (!range.isSane())
		throw std::runtime_error("Interaction::=("+toString(range)+") not sane!");
#endif
	// clear current interactions
	basePairs.clear();

	// undo seed information
	INTARNA_CLEANUP(seed);

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

bool
Interaction::
operator == ( const Interaction &i ) const
{
	return 	   s1 == i.s1
			&& s2 == i.s2
			&& E_equal( energy, i.energy )
			&& basePairs == i.basePairs
			&& (seed == i.seed || *seed == *(i.seed))
			;
}

////////////////////////////////////////////////////////////////////////////

std::string
Interaction::
dotBar( const Interaction & i, const bool fullLength )
{
#if INTARNA_IN_DEBUG_MODE
	if (!i.isValid())
		throw std::runtime_error("Interaction::dotBar("+toString(i)+") not valid!");
#endif
	// compile dot-bar representation
	return toString(i.s1->getInOutIndex(fullLength?0:i.basePairs.begin()->first))
			+ (fullLength?std::string(i.basePairs.begin()->first, '.' ):"") // leading unpaired s1
			+ dotSomething(i.basePairs.begin(), i.basePairs.end(), true, '|') // s1 structure
			+ (fullLength?std::string(i.s1->size() - i.basePairs.rbegin()->first -1, '.' ):"") // trailing unpaired s1
			+ "&"
			+ toString(i.s2->getInOutIndex(fullLength?0:i.basePairs.rbegin()->second))
			+ (fullLength?std::string(i.basePairs.rbegin()->second, '.' ):"") // trailing unpaired s2
			+ dotSomething(i.basePairs.rbegin(), i.basePairs.rend(), false, '|') // s2 structure
			+ (fullLength?std::string( i.s2->size() - i.basePairs.begin()->second -1, '.' ):"") // leading unpaired s2
			;
}

////////////////////////////////////////////////////////////////////////////

std::string
Interaction::
dotBracket( const Interaction & i, const char symOpen, const char symClose, const bool fullLength )
{
#if INTARNA_IN_DEBUG_MODE
	if (!i.isValid())
		throw std::runtime_error("Interaction::dotBracket("+toString(i)+") not valid!");
#endif

		// compile dot-bracket representation
	return	(fullLength?std::string(i.basePairs.begin()->first, '.' ):"") // leading unpaired s1
				+ dotSomething(i.basePairs.begin(), i.basePairs.end(), true, symOpen) // s1 structure
			+ (fullLength?std::string(i.s1->size() - i.basePairs.rbegin()->first -1, '.' ):"") // trailing unpaired s1
				+"&"
			+ (fullLength?std::string(i.basePairs.rbegin()->second, '.' ):"") // leading unpaired s2
				+ dotSomething(i.basePairs.rbegin(), i.basePairs.rend(), false, symClose) // s2 structure
			+ (fullLength?std::string(i.s2->size() - i.basePairs.begin()->second -1, '.' ):"") // trailing unpaired s2
			;
}

////////////////////////////////////////////////////////////////////////////

} // namespace
