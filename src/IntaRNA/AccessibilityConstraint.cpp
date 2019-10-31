
#include "IntaRNA/AccessibilityConstraint.h"

namespace IntaRNA {

////////////////////////////////////////////////////////////////////////

// the marker for unconstrained positions in dot-bracket notation
const char AccessibilityConstraint::dotBracket_unconstrained = '.';
// the marker for blocked positions in dot-bracket notation
const char AccessibilityConstraint::dotBracket_blocked = 'b';
// the marker for accessible positions in dot-bracket notation
const char AccessibilityConstraint::dotBracket_accessible = 'x';
// the marker for paired positions in dot-bracket notation
const char AccessibilityConstraint::dotBracket_paired = 'p';

const std::string AccessibilityConstraint::dotBracket_constraints =
					 toString(AccessibilityConstraint::dotBracket_accessible)
					+toString(AccessibilityConstraint::dotBracket_blocked)
					+toString(AccessibilityConstraint::dotBracket_paired)
					;

const std::string AccessibilityConstraint::dotBracketAlphabet =
					 toString(AccessibilityConstraint::dotBracket_unconstrained)
					+AccessibilityConstraint::dotBracket_constraints
					;
const std::string AccessibilityConstraint::regionIndexList =
					"["+AccessibilityConstraint::dotBracket_constraints+"]:"
					+ R"((-?[123456789]\d*--?[123456789]\d*,)*-?[123456789]\d*--?[123456789]\d*)";

const boost::regex AccessibilityConstraint::regex(
					"(["
					+ AccessibilityConstraint::dotBracketAlphabet
					+"]+)|"+
					+"("+AccessibilityConstraint::regionIndexList
						+"(,"+AccessibilityConstraint::regionIndexList+")*"
					+")"
					);

// D | Dm1.9 | Db-0.7 | Dm1.9b-0.7
// Z | Zb0.8
// W
const boost::regex AccessibilityConstraint::regexShapeMethod(R"(^((D(m-?(\d+|\d+.\d+))?(b-?(\d+|\d+.\d+))?)|(Z(b(\d+|\d+.\d+))?)|(W))$)");

// M
// C0.2
// S
// Ls0.5i0.5
// Os0.5i-0.5
const boost::regex AccessibilityConstraint::regexShapeConversion(R"(^((M)|(C(\d+|\d+.\d+)?)|(S)|((L|O)(s-?(\d+|\d+.\d+)i-?(\d+|\d+.\d+))?))$)");


////////////////////////////////////////////////////////////////////////

AccessibilityConstraint::
AccessibilityConstraint(
			  const RnaSequence & seq
			, const std::string& stringEncoding
			, const size_t maxBpSpan_
			, const std::string & shapeFile_
			, const std::string & shapeMethod_
			, const std::string & shapeConversion_
			)
	:
	length(seq.size()),
	maxBpSpan( maxBpSpan_==0 ? length : std::min(maxBpSpan_,length) ),
	shapeFile(shapeFile_),
	shapeMethod(shapeFile_.empty() ? "" : shapeMethod_),
	shapeConversion(shapeFile_.empty() ? "" : shapeConversion_),
	blocked(),
	accessible(),
	paired()
{
#if INTARNA_IN_DEBUG_MODE
	if (!boost::regex_match( stringEncoding, AccessibilityConstraint::regex, boost::match_perl )) {
		throw std::runtime_error("AccessibilityConstraint("+stringEncoding+") does not match its encoding regular expression");
	}
	if (!shapeFile.empty()) {
		if (!boost::regex_match( shapeMethod, AccessibilityConstraint::regexShapeMethod, boost::match_perl )) {
			throw std::runtime_error("AccessibilityConstraint(shapeMethod="+shapeMethod+") does not match its encoding regular expression");
		}
		if (!boost::regex_match( shapeConversion, AccessibilityConstraint::regexShapeConversion, boost::match_perl )) {
			throw std::runtime_error("AccessibilityConstraint(shapeConversion="+shapeConversion+") does not match its encoding regular expression");
		}
	}
#endif

	// check if dot-bracket-encoding
	if (stringEncoding.size() == length && stringEncoding.find_first_not_of(dotBracketAlphabet) == std::string::npos) {
		// screen for blocked regions
		screenDotBracket( stringEncoding, dotBracket_blocked, blocked );
		// screen for accessible regions
		screenDotBracket( stringEncoding, dotBracket_accessible, accessible );
		// screen for accessible regions
		screenDotBracket( stringEncoding, dotBracket_paired, paired );
	} else
		// if non-empty : assume region-based encoding
		if (stringEncoding.size() > 0 && stringEncoding.find_first_not_of(dotBracketAlphabet) != std::string::npos)
	{
		// parse region encoding
		size_t start=0, end=stringEncoding.find_first_of(dotBracket_constraints, start+1);
		while( start < stringEncoding.size() ) {
			// deal with search result
			end = (end == std::string::npos ? stringEncoding.size()+1 : end);
			// parse index range encoding
			IndexRangeList curRangeList(stringEncoding.substr(start+2,end-start-3), false, &seq);
			// store range
			switch( stringEncoding.at(start) ) {
			case dotBracket_accessible :
				for ( auto r=curRangeList.begin(); r!=curRangeList.end(); r++) {
					accessible.insert( *r );
				}
				break;
			case dotBracket_blocked :
				for ( auto r=curRangeList.begin(); r!=curRangeList.end(); r++) {
					blocked.insert( *r );
				}
				break;
			case dotBracket_paired :
				for ( auto r=curRangeList.begin(); r!=curRangeList.end(); r++) {
					paired.insert( *r );
				}
				break;
			default:
				throw std::runtime_error("AccessibilityConstraint() : unexpected constraint encoding "+toString(stringEncoding.at(start)));
			}
			// next region encoding
			start = end;
			end = stringEncoding.find_first_of(dotBracket_constraints, start+1);
		}
	}

}

////////////////////////////////////////////////////////////////////////

void
AccessibilityConstraint::
screenDotBracket( const std::string& dotBracket
				, const char marker
				, IndexRangeList & storage )
{
	// temporary variable holding the start of the current region
	size_t lastRegionStart = std::string::npos;
	// screen for consecutive marker regions
	for (size_t i=0; i<dotBracket.size(); i++) {
		// check if marked as blocked
		if (dotBracket.at(i)==marker) {
			// check if start of a region
			if (lastRegionStart == std::string::npos) {
				// store start
				lastRegionStart = i;
			}
		} else {
			// check if we just left a region
			if (lastRegionStart != std::string::npos) {
				// store region
				storage.push_back( IndexRange( lastRegionStart, i-1 ) );
				// reset marker
				lastRegionStart = std::string::npos;
			}
		}
	}
	// check if last region spans till the end of the strign
	if (lastRegionStart != std::string::npos) {
		// store region
		storage.push_back( IndexRange( lastRegionStart, dotBracket.size()-1 ) );
		// reset marker
		lastRegionStart = std::string::npos;
	}
}


////////////////////////////////////////////////////////////////////////

} // namespace
