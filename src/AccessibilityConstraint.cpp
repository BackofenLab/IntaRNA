
#include "AccessibilityConstraint.h"

////////////////////////////////////////////////////////////////////////

// the marker for blocked positions in dot-bracket notation
const char AccessibilityConstraint::dotBracket_blocked = 'b';
// the marker for accessible positions in dot-bracket notation
const char AccessibilityConstraint::dotBracket_accessible = 'x';

const std::string AccessibilityConstraint::dotBracketAlphabet = ".()"
					+toString(AccessibilityConstraint::dotBracket_accessible)
					+toString(AccessibilityConstraint::dotBracket_blocked);


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
