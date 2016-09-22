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


/////////////////////////////////////////////////////////////////////////////

RnaSequence::RnaSequence(
		const std::string & id
		, const std::string & seqString )
 :
	id(id)
	, seqString(getUpperCase(seqString))
	, seqCode(getCodeForString(this->seqString))
	, ambiguous(seqString.find_first_of("nN")!=std::string::npos)
{
}

/////////////////////////////////////////////////////////////////////////////

RnaSequence::~RnaSequence()
{
}

/////////////////////////////////////////////////////////////////////////////

const std::string&
RnaSequence::
getId() const
{
	return id;
}

/////////////////////////////////////////////////////////////////////////////

size_t
RnaSequence::
size() const
{
	return seqString.size();
}

/////////////////////////////////////////////////////////////////////////////

const
RnaSequence::
String_type&
RnaSequence::
asString() const
{
	return seqString;
}

/////////////////////////////////////////////////////////////////////////////

const
RnaSequence::
CodeSeq_type&
RnaSequence::
asCodes() const
{
	return seqCode;
}

/////////////////////////////////////////////////////////////////////////////

bool
RnaSequence::
isAmbiguous() const
{
	return ambiguous;
}

/////////////////////////////////////////////////////////////////////////////


RnaSequence::
String_type
RnaSequence::
getUpperCase( const std::string & seqString )
{
	// create container to fill
	String_type seqRet(seqString.size(),'_');

	for (size_t i=0; i<seqString.size(); ++i)
	{
		// get upper case characters
		seqRet[i] = std::toupper(seqString.at(i),codeLocale);
	}

	return seqRet;
}


/////////////////////////////////////////////////////////////////////////////


RnaSequence::
CodeSeq_type
RnaSequence::
getCodeForString( const String_type& seqString )
{
	// create container to fill and init with 'X'
	CodeSeq_type seqCode(seqString.size());

	for (size_t i=0; i<seqString.size(); ++i)
	{
		seqCode[i] = getCodeForChar( seqString.at(i) );
	}

	return seqCode;
}


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
isValidSequence( const std::string & sequence )
{
	// check whether or not the string contains unsupported characters
	return sequence.find_first_not_of(SequenceAlphabet) == std::string::npos;
}

/////////////////////////////////////////////////////////////////////////////

