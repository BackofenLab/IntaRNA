
#include "IntaRNA/RnaSequence.h"
#include <stdexcept>

namespace IntaRNA {


// setup static members
std::locale RnaSequence::codeLocale = std::locale();

// setup allowed alphabet in lower and upper case
const std::string RnaSequence::SequenceAlphabet	= "ACGUN";
const std::string RnaSequence::SequenceAlphabetIUPAC	= "aAuUcCgGtTrRyYsSwWkKmMbBdDhHvVnN";

// setup GU base pair code information
int RnaSequence::bpGUcodes[] = { BP_pair[RnaSequence::getCodeForChar('G')][RnaSequence::getCodeForChar('U')]
				,BP_pair[RnaSequence::getCodeForChar('U')][RnaSequence::getCodeForChar('G')]};


////////////////////////////////////////////////////////////////////////////

const size_t RnaSequence::lastPos = std::string::npos;

/////////////////////////////////////////////////////////////////////////////


} // namespace
