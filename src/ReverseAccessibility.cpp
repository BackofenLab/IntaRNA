/*
 * ReverseAccessibility.cpp
 *
 *  Created on: 30.06.2014
 *      Author: Mmann
 */

#include "ReverseAccessibility.h"


////////////////////////////////////////////////////////////////////////////

ReverseAccessibility::ReverseAccessibility( Accessibility & origAcc )
 :
	Accessibility( origAcc.getSequence(), origAcc.getMaxLength(), &origAcc.getAccConstraint() )
	, origAcc(origAcc)
	, seqReversed( seq.getId(), getReversedString(seq) )
	, accConstrReversed( origAcc.getAccConstraint(), true )
{
}

////////////////////////////////////////////////////////////////////////////

ReverseAccessibility::~ReverseAccessibility() {
}

////////////////////////////////////////////////////////////////////////////

std::string
ReverseAccessibility::
getReversedString( const RnaSequence & seq )
{
	const RnaSequence::String_type& seqAsString = seq.asString();
	// create string container
	std::string revStr(seqAsString.size(),'_');
	// copy in reversed order
	std::string::reverse_iterator revStrIt = revStr.rbegin();
	for (RnaSequence::String_type::const_iterator nuc = seqAsString.begin();
			nuc != seqAsString.end(); ++nuc,++revStrIt)
	{
		*revStrIt = *nuc;
	}
	// return reversed string
	return revStr;
}

////////////////////////////////////////////////////////////////////////////

std::string
ReverseAccessibility::
getReversedString( const std::string & seqAsString )
{
	// create string container
	std::string revStr(seqAsString.size(),'_');
	// copy in reversed order
	std::string::reverse_iterator revStrIt = revStr.rbegin();
	for (RnaSequence::String_type::const_iterator nuc = seqAsString.begin();
			nuc != seqAsString.end(); ++nuc,++revStrIt)
	{
		*revStrIt = *nuc;
	}
	// return reversed string
	return revStr;
}

////////////////////////////////////////////////////////////////////////////

const RnaSequence &
ReverseAccessibility::
getSequence() const
{
	// access reversed sequence
	return seqReversed;
}


////////////////////////////////////////////////////////////////////////////

const AccessibilityConstraint&
ReverseAccessibility::
getAccConstraint() const
{
	// access reversed accessibility constraint
	return accConstrReversed;
}


////////////////////////////////////////////////////////////////////////////


E_type
ReverseAccessibility::
getED( const size_t from, const size_t to ) const
{
	// check indices
	checkIndices(from,to);
	// reversed ED access
	return origAcc.getED( seq.size()-to-1, seq.size()-from-1 );
}


////////////////////////////////////////////////////////////////////////////

const Accessibility &
ReverseAccessibility::
getAccessibilityOrigin() const
{
	return origAcc;
}

////////////////////////////////////////////////////////////////////////////

size_t
ReverseAccessibility::
getReversedIndex( const size_t i ) const
{
	// check indices
	checkIndices(i,i);
	// compute reverse index
	return this->seq.size() -i -1;
}

////////////////////////////////////////////////////////////////////////////
