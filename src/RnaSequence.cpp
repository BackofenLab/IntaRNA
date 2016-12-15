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

RnaSequence::RnaSequence(
		const std::string & id
		, const std::string & seqString )
 :
	id(id)
	, seqString(getUpperCase(seqString))
	, seqCode(getCodeForString(this->seqString))
	, ambiguous(seqString.find_first_of("nN")!=std::string::npos)
{
#if IN_DEBUG_MODE
	if (id.size() == 0) {
		throw std::runtime_error("RnaSequence::RnaSequence : id empty");
	}
	if (seqString.size() == 0) {
		throw std::runtime_error("RnaSequence::RnaSequence : seqString empty");
	}
#endif

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

bool
RnaSequence::
isAmbiguous( const size_t i ) const
{
	// TODO ensure sequence is upper case and reduce the test accordingly
	return this->seqString.at(i) == 'N' || this->seqString.at(i) == 'n';
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

std::ostream&
operator<<(std::ostream& out, const RnaSequence& rna)
{
	// add ID(SEQUENCE) to stream
	out <<rna.id <<'(' <<rna.asString() <<')';
	return out;
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

