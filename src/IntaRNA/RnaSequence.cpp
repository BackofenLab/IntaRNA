/*
 * Sequence.cpp
 *
 *  Created on: 24.06.2014
 *      Author: Mmann
 */

#include "IntaRNA/RnaSequence.h"
#include <stdexcept>

namespace IntaRNA {


// setup static members
std::locale RnaSequence::codeLocale = std::locale();

// setup allowed alphabet in lower and upper case
const std::string RnaSequence::SequenceAlphabet	= "ACGUN";
const std::string RnaSequence::SequenceAlphabetIUPAC	= "aAuUcCgGtTrRyYsSwWkKmMbBdDhHvVnN";


////////////////////////////////////////////////////////////////////////////

const size_t RnaSequence::lastPos = std::string::npos;

/////////////////////////////////////////////////////////////////////////////


} // namespace
