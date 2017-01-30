
#include "AccessibilityFromStream.h"

#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>

/////////////////////////////////////////////////////////////////////////

AccessibilityFromStream::
AccessibilityFromStream(
		const RnaSequence& sequence
		, const size_t maxLength
		, const AccessibilityConstraint * const accConstraint
		, std::istream & inStream
		, const InStreamType inStreamType
		, const E_type RT
		)
 :	Accessibility( sequence, maxLength, accConstraint )
	, edValues( getSequence().size(), getSequence().size(), 0, getMaxLength() )
{
	switch( inStreamType ) {

	case Pu_RNAplfold_Text :
		parsePu_RNAplfold_Text( inStream, RT );
		break;

	}
}

/////////////////////////////////////////////////////////////////////////

AccessibilityFromStream::
~AccessibilityFromStream()
{
}

/////////////////////////////////////////////////////////////////////////


void
AccessibilityFromStream::
parsePu_RNAplfold_Text( std::istream & inStream, const E_type RT )
{
	// assume VRNA v2* style = matrix with maxLength rows

	// skip leading white spaces
	inStream >>std::skipws;

	// parse first comment line
	std::string line;
	if ( !std::getline( inStream, line ) ) {
		throw std::runtime_error("AccessibilityFromStream::parsePu_RNAplfold_Text() : nothing readable");
	}
	if ( ! boost::regex_match(line,boost::regex("^\\s*#unpaired probabilities\\s*$"), boost::match_perl) ) {
		throw std::runtime_error("AccessibilityFromStream::parsePu_RNAplfold_Text() : first line != expected '#unpaired probabilities' header");
	}
	// parse second line = available lengths
	if ( !std::getline( inStream, line ) ) {
		throw std::runtime_error("AccessibilityFromStream::parsePu_RNAplfold_Text() : length header (2nd line) not found");
	}
	if ( ! boost::regex_match(line,boost::regex("^\\s*#i.\\s+l=1(\\s+\\d+)*\\s*$"), boost::match_perl) ) {
		throw std::runtime_error("AccessibilityFromStream::parsePu_RNAplfold_Text() : second line is no proper lengths header");
	}
	// check if maxLength <= max available length
	size_t cutEnd   = line.find_last_of("1234567890");
	size_t cutStart = line.find_last_not_of("1234567890", cutEnd );
	size_t maxAvailLength = boost::lexical_cast<size_t>( line.substr(cutStart+1,cutEnd-cutStart));
	if (maxAvailLength < getMaxLength()) {
		throw std::runtime_error("AccessibilityFromStream::parsePu_RNAplfold_Text() : available maximal window length "
				+toString(maxAvailLength)+" is smaller than required length "+toString(getMaxLength()));
	}

	// end of ED window (= first column in file)
	size_t j = 0, lastJ = 0;
	while ( ! inStream.fail() && j < edValues.size2() ) {
		// read first column = end of window = j
		if ( inStream >> j ) {
			// check if lines are consecutive
			if ( j != lastJ+1 ) {
				throw std::runtime_error("AccessibilityFromStream::parsePu_RNAplfold_Text() : non-consecutive line i="+toString(j)+" was preceeded by "+toString(lastJ));
			}
			// check if we line exceeds targeted length
			if ( j > getSequence().size() ) {
				LOG(INFO) <<"AccessibilityFromStream::parsePu_RNAplfold_Text() : more lines found than sequence is long.. sure this is the correct file for this sequence?";
				// stop parsing
				break;
			}
			if ( j == lastJ ) {
				throw std::runtime_error("AccessibilityFromStream::parsePu_RNAplfold_Text() : duplicate for i="
						+ toString(lastJ));
			}
		} else {
			throw std::runtime_error("AccessibilityFromStream::parsePu_RNAplfold_Text() : could not read next line after parsing "
					+ toString(lastJ)+" lines");
		}

		// parse probabilities for this line and store
		double curProb;
		size_t minI = j - std::min( j, getMaxLength() );
		for ( size_t i = j; i>minI; i--) {
			if ( inStream >>curProb ) {
				if (curProb < 0.0 || curProb > 1.0) {
					throw std::runtime_error("AccessibilityFromStream::parsePu_RNAplfold_Text() : in line i="+toString(j)
							+" : the "+toString(j+1-i)+". value = "+toString(curProb)+" is no probability in [0,1]");
				}
				edValues( i-1, j-1 ) = curProb > 0
										? std::min<E_type>(ED_UPPER_BOUND, - RT * std::log( curProb ))
										: ED_UPPER_BOUND;
			} else {
				throw std::runtime_error("AccessibilityFromStream::parsePu_RNAplfold_Text() : in line i="+toString(j)
						+" : could not parse the "+toString(j+1-i)+". probability");
			}
		}
		// check if full line was already parsed
		if (j < maxAvailLength || minI > 0) {
			// skip rest till end of line
			inStream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		}

		// update parsing information for next run
		lastJ = j;
	}

	// check if all needed data was parsed
	if (lastJ < edValues.size2()) {
		throw std::runtime_error("AccessibilityFromStream::parsePu_RNAplfold_Text() : could only parse "
				+toString(lastJ)+" lines, but "+toString(edValues.size2())
				+" expected (length of sequence "+getSequence().getId()+")");
	}

}

/////////////////////////////////////////////////////////////////////////
