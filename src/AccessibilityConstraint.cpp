
#include "AccessibilityConstraint.h"

#include <algorithm>

////////////////////////////////////////////////////////////////////////

// the marker for blocked positions in dot-bracket notation
const char AccessibilityConstraint::dotBracket_blocked = 'b';
// the marker for accessible positions in dot-bracket notation
const char AccessibilityConstraint::dotBracket_accessible = 'x';

const std::string AccessibilityConstraint::dotBracketAlphabet = ".()"
					+toString(AccessibilityConstraint::dotBracket_accessible)
					+toString(AccessibilityConstraint::dotBracket_blocked);


////////////////////////////////////////////////////////////////////////

AccessibilityConstraint::
AccessibilityConstraint( const size_t length )
	:
	length(length),
	blocked(),
	accessible()
{
}

////////////////////////////////////////////////////////////////////////

AccessibilityConstraint::
AccessibilityConstraint( const std::string& dotBracket )
	:
	length(dotBracket.size()),
	blocked(),
	accessible()
{
#if IN_DEBUG_MODE
	if (dotBracket.find_first_not_of(dotBracketAlphabet)!=std::string::npos) {
		throw std::runtime_error("AccessibilityConstraint("+dotBracket+") contains unsupported characters");
	}
#endif

	// check for base pairs (not implemented yet)
	if (dotBracket.find_first_of("()") != std::string::npos) {
		NOTIMPLEMENTED("AccessibilityConstraint(dotBracket) contains base pairs... currently only '."+toString(dotBracket_accessible)+toString(dotBracket_blocked)+"' implemented");
	}

	// screen for blocked regions
	screenDotBracket( dotBracket, dotBracket_blocked, blocked );
	// screen for accessible regions
	screenDotBracket( dotBracket, dotBracket_accessible, accessible );
}

////////////////////////////////////////////////////////////////////////

AccessibilityConstraint::
AccessibilityConstraint( const AccessibilityConstraint& toCopy
			, const bool reverseIndices)
	:
	length(toCopy.length)
	, blocked(toCopy.blocked)
	, accessible(toCopy.accessible)
{
	// TODO copy structure constraints etc.

	if (reverseIndices) {

		// reverse blocked
		for (IndexRangeList::iterator it=blocked.begin(); it!= blocked.end(); it++) {
			it->from = length - (it->from +1);
			it->to = length - (it->to +1);
		}
		// sort again
		std::sort(accessible.begin(), accessible.end());

		// reverse accessible
		for (IndexRangeList::iterator it=accessible.begin(); it!= accessible.end(); it++) {
			it->from = length - (it->from +1);
			it->to = length - (it->to +1);
		}
		// sort again
		std::sort(accessible.begin(), accessible.end());

		// TODO reverse structure constraints
	}
}

////////////////////////////////////////////////////////////////////////

AccessibilityConstraint::~AccessibilityConstraint()
{
}

////////////////////////////////////////////////////////////////////////

bool
AccessibilityConstraint::
isMarkedBlocked(const size_t i) const
{
	return blocked.covers(i);
}

////////////////////////////////////////////////////////////////////////

bool
AccessibilityConstraint::
isMarkedAccessible(const size_t i) const
{
	return accessible.covers(i);
}

////////////////////////////////////////////////////////////////////////

bool
AccessibilityConstraint::
isUnconstrained( const size_t i ) const
{
	return isEmpty()
	// TODO handle base pairing constraints etc..
			|| !(isMarkedAccessible(i) || isMarkedBlocked(i));
}

////////////////////////////////////////////////////////////////////////

bool
AccessibilityConstraint::
isAccessible( const size_t i ) const
{
	// TODO handle base pairing constraints etc.
	return isEmpty()
			|| !isMarkedBlocked(i);
}

////////////////////////////////////////////////////////////////////////

bool
AccessibilityConstraint::
isEmpty() const
{
	// TODO CHECK FOR BASE PAIRS ETC
	// check if any constrained regions given
	return (accessible.size() + blocked.size()) == 0;
}

////////////////////////////////////////////////////////////////////////

char
AccessibilityConstraint::
getVrnaDotBracket(const size_t i) const
{
	// check if to be accessible or blocked (==unstructured)
	if (isMarkedAccessible(i) || isMarkedBlocked(i)) {
		return 'x';
	}

	// TODO add base pair handling etc.

	return '.';
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

AccessibilityConstraint &
AccessibilityConstraint::
operator= ( const AccessibilityConstraint & c )
{
	// copy data
	length = c.length;
	blocked = c.blocked;
	accessible = c.accessible;
	// TODO copy structure constraints etc.

	return *this;
}

////////////////////////////////////////////////////////////////////////
