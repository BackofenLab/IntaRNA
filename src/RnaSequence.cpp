/*
 * Sequence.cpp
 *
 *  Created on: 24.06.2014
 *      Author: Mmann
 */

#include "RnaSequence.h"
#include <stdexcept>


extern "C" {
	#include "ViennaRNA/pair_mat.h"
}



// setup static members
std::locale RnaSequence::codeLocale = std::locale();

// setup allowed alphabet in lower and upper case
const std::string RnaSequence::SequenceAlphabet	= "aAuUcCgGnN";


////////////////////////////////////////////////////////////////////////////

const size_t RnaSequence::lastPos = std::string::npos;

/////////////////////////////////////////////////////////////////////////////

RnaSequence::
Code_type
RnaSequence::
getCodeForChar( const char nucleotide )
{
	// check if nucleotide character is NOT supported
	if (SequenceAlphabet.find(nucleotide) == std::string::npos)
	{
		// raise exception
		throw std::runtime_error("RnaSequence::getCodeForChar() : unsupported nucleotide character '"+toString(nucleotide)+"' in sequence");
	}
	// otherwise get encoding:
	// use nucleotide character encoding from Vienna RNA package
	return (Code_type)encode_char(nucleotide);
}

/////////////////////////////////////////////////////////////////////////////

bool
RnaSequence::
areComplementary( const RnaSequence & s1, const RnaSequence & s2,
					const size_t p1, const size_t p2 )
{
	// check if valid positions
	if (p1<s1.size() && p2<s2.size()) {
		// check via VRNA util
		return BP_pair[s1.seqCode.at(p1)][s2.seqCode.at(p2)] > 0;
	} else {
		throw std::runtime_error("RnaSequence::areComplementary : index positions p1/p2 ("
				+ toString(p1)+"/"+toString(p2)
				+ ") are out of bounds s1/s2 ("
				+ toString(s1.size())+"/"+toString(s2.size())
				+")"
				);
	}
}

/////////////////////////////////////////////////////////////////////////////

