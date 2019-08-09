
#include "IntaRNA/AccessibilityFromStream.h"

#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>

namespace IntaRNA {

/////////////////////////////////////////////////////////////////////////

AccessibilityFromStream::
AccessibilityFromStream(
		const RnaSequence& sequence
		, const size_t maxLength
		, const AccessibilityConstraint * const accConstraint
		, std::istream & inStream
		, const InStreamType inStreamType
		, const Z_type RT
		)
 :	Accessibility( sequence, maxLength, accConstraint )
	, edValues()
	, availMaxLength( Accessibility::getMaxLength() )
{
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) <<"reading accessibility from input..."; }

	if (accConstraint != NULL && !accConstraint->isEmpty()) {
		INTARNA_NOT_IMPLEMENTED("AccessibilityFromStream: accessibility constraints not supported for direct accessibility input");
	}
	switch( inStreamType ) {

	case Pu_RNAplfold_Text :
		parsePu_RNAplfold_text( inStream, RT );
		break;

	case ED_RNAplfold_Text :
		parseED_RNAplfold_text( inStream );
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
parseRNAplfold_text( std::istream & inStream, const Z_type RT, const bool parseProbs )
{
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) <<"parsing "<<(parseProbs?"unpaired probabilities":"accessibility values")<<" from RNAplfold"<<(parseProbs?"":"-like")<<" input ..."; }
	// time logging
	TIMED_FUNC_IF(timerObj, VLOG_IS_ON(9));

	// assume VRNA v2* style = matrix with maxLength rows

	// skip leading white spaces
	inStream >>std::skipws;

	// parse first comment line
	std::string line;
	if ( !std::getline( inStream, line ) ) {
		throw std::runtime_error("AccessibilityFromStream::parseRNAplfold_text() : nothing readable");
	}
	if ( ! boost::regex_match(line,boost::regex("^#[\\w\\s]+$"), boost::match_perl) ) {
		throw std::runtime_error("AccessibilityFromStream::parseRNAplfold_text() : first line != expected header line starting with '#'");
	}

	// parse second line = available lengths
	if ( !std::getline( inStream, line ) ) {
		throw std::runtime_error("AccessibilityFromStream::parseRNAplfold_text() : length header (2nd line) not found");
	}
	if ( ! boost::regex_match(line,boost::regex("^\\s*#i.\\s+l=1(\\s+\\d+)*\\s*$"), boost::match_perl) ) {
		throw std::runtime_error("AccessibilityFromStream::parseRNAplfold_text() : second line is no proper lengths header");
	}
	// check if maxLength <= max available length
	size_t cutEnd   = line.find_last_of("1234567890");
	size_t cutStart = line.find_last_not_of("1234567890", cutEnd );
	size_t maxAvailLength = boost::lexical_cast<size_t>( line.substr(cutStart+1,cutEnd-cutStart));
	if (maxAvailLength < getMaxLength()) {
#if INTARNA_MULITHREADING
		#pragma omp critical(intarna_omp_logOutput)
#endif
		{ LOG(INFO) <<"initializing ED data for sequence '"<<getSequence().getId()<<" : available maximal window length "
				<<maxAvailLength<<" is smaller than maximal interaction length "<<getMaxLength()
				<<" : reducing maximal interaction length to "<<maxAvailLength;}
		// reducing maximal interaction length
		availMaxLength = maxAvailLength;
		assert( getMaxLength() == maxAvailLength ); // ensure overwrite is working
	}

	// resize data structure to fill
	edValues.resize( getSequence().size(), getSequence().size(), 0, getMaxLength() );

	// TODO rewrite to support "nan" and "inf" parsing via boost::spirit::qi
	// http://stackoverflow.com/questions/11420263/is-it-possible-to-read-infinity-or-nan-values-using-input-streams

	// end of ED window (= first column in file)
	size_t j = 0, lastJ = 0;
	while ( ! inStream.fail() && j < edValues.size2() ) {
		// read first column = end of window = j
		if ( inStream >> j ) {
			// check if lines are consecutive
			if ( j != lastJ+1 ) {
				throw std::runtime_error("AccessibilityFromStream::parseRNAplfold_text() : non-consecutive line i="+toString(j)+" was preceeded by "+toString(lastJ));
			}
			// check if we line exceeds targeted length
			if ( j > getSequence().size() ) {
#if INTARNA_MULITHREADING
				#pragma omp critical(intarna_omp_logOutput)
#endif
				{ LOG(INFO) <<"AccessibilityFromStream::parseRNAplfold_text() : more lines found than sequence is long.. sure this is the correct file for this sequence?"; }
				// stop parsing
				break;
			}
			if ( j == lastJ ) {
				throw std::runtime_error("AccessibilityFromStream::parseRNAplfold_text() : duplicate for i="
						+ toString(lastJ));
			}
		} else {
			throw std::runtime_error("AccessibilityFromStream::parseRNAplfold_text() : could not read next line start (integer i) after parsing "
					+ toString(lastJ)+" lines of values");
		}

		// parse probabilities or EDs for this line and store
		double curVal;
		size_t minI = j - std::min( j, getMaxLength() );
		for ( size_t i = j; i>minI; i--) {
			if ( inStream >>curVal ) {
				// check if we parse probabilities
				if (parseProbs) {
					if (curVal < 0.0 || curVal > 1.0) {
						throw std::runtime_error("AccessibilityFromStream::parseRNAplfold_Text(Pu) : in line i="+toString(j)
								+" : the "+toString(j+1-i)+". value = "+toString(curVal)+" is no probability in [0,1]");
					}
					edValues( i-1, j-1 ) = curVal > 0
											? std::min<E_type>(ED_UPPER_BOUND, Z_2_E( - RT * Z_log( Z_type(curVal) ) ))
											: ED_UPPER_BOUND;
				}
				// or ED values (in kcal/mol)
				else {
					if (curVal < 0.0) {
						throw std::runtime_error("AccessibilityFromStream::parseRNAplfold_Text(ED) : in line i="+toString(j)
								+" : the "+toString(j+1-i)+". value = "+toString(curVal)+" is no ED value >= 0");
					}
					edValues( i-1, j-1 ) = std::min<E_type>(ED_UPPER_BOUND, Ekcal_2_E(curVal));
				}
			} else {
				throw std::runtime_error("AccessibilityFromStream::parseRNAplfold_text() : in line i="+toString(j)
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
		throw std::runtime_error("AccessibilityFromStream::parseRNAplfold_text() : could only parse "
				+toString(lastJ)+" lines, but "+toString(edValues.size2())
				+" expected (length of sequence "+getSequence().getId()+")");
	}

}

/////////////////////////////////////////////////////////////////////////

} // namespace
