/*
 * ReverseAccessibility.cpp
 *
 *  Created on: 30.06.2014
 *      Author: Mmann
 */

#include "IntaRNA/ReverseAccessibility.h"

namespace IntaRNA {

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

} // namespace
