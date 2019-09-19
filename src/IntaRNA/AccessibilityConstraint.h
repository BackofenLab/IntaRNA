
#ifndef INTARNA_ACCESSIBILITYCONSTRAINT_H_
#define INTARNA_ACCESSIBILITYCONSTRAINT_H_

#include "IntaRNA/IndexRangeList.h"
#include "IntaRNA/RnaSequence.h"

#include <boost/regex.hpp>

#include <utility>
#include <vector>

namespace IntaRNA {

/**
 * Represents the constraints for accessibility computation, ie. sequence
 * regions that are known to be blocked (not available to interaction) or
 * for sure accessible.
 *
 * These constraints are than incorporated into the computation of the
 * accessibility (ED) values to be used for interaction prediction.
 *
 * @author Martin Mann
 *
 */
class AccessibilityConstraint {

public:

	//! the marker for unconstrained positions in dot-bracket notation
	static const char dotBracket_unconstrained;
	//! the marker for blocked positions in dot-bracket notation
	static const char dotBracket_blocked;
	//! the marker for accessible positions in dot-bracket notation
	static const char dotBracket_accessible;
	//! the marker for paired positions in dot-bracket notation
	static const char dotBracket_paired;

	//! collection of all allowed encoding of constraints
	static const std::string dotBracket_constraints;

	//! the alphabet to encode accessibility constraints in dot-bracket notation
	static const std::string dotBracketAlphabet;

	//! regular expression that encodes a single region encoding as index list
	static const std::string regionIndexList;

	//! the regular expression that marks a non-empty valid constraint encoding
	static const boost::regex regex;

	//! the regular expression to be matched by shapeMethod encodings
	static const boost::regex regexShapeMethod;

	//! the regular expression to be matched by shapeConversion encodings
	static const boost::regex regexShapeConversion;

public:

	/**
	 * Empty constraint construction
	 * @param length length of the constrained sequence
	 * @param maxBpSpan the maximal base pair span to be used for accessibility
	 *        computation; set to 0 for full sequence length
	 * @param shapeFile if non-empty: name of the SHAPE reactivity data file to
	 *        be used for accessibility computation
	 * @param shapeMethod if non-empty: encoding of the method to be used to
	 *        convert SHAPE reactivity data to pseudo energies
	 * @param shapeConversion if non-empty: encoding of the method to be used to
	 *        convert SHAPE reactivity data to pairing probabilities
	 */
	AccessibilityConstraint(  const size_t length
							, const size_t maxBpSpan
							, const std::string & shapeFile
							, const std::string & shapeMethod
							, const std::string & shapeConversion
							);

	/**
	 * Copy construction
	 * @param toCopy the constraint to copy
	 * @param reverseIndices whether or not to reverse indexing (eg to be used
	 *        for ReverseAccessibility)
	 */
	AccessibilityConstraint( const AccessibilityConstraint& toCopy
						, const bool reverseIndices = false );

	/**
	 * Constraint construction based on VRNA-like dot-bracket encoding
	 * @param seq the constrained sequence
	 * @param dotBracket the constraint encoding in VRNA-like dot-bracket encoding
	 * @param maxBpSpan the maximal base pair span to be used for accessibility
	 *        computation; set to 0 for full sequence length
	 * @param shapeFile if non-empty: name of the SHAPE reactivity data file to
	 *        be used for accessibility computation
	 * @param shapeMethod if non-empty: encoding of the method to be used to
	 *        convert SHAPE reactivity data to pseudo energies
	 * @param shapeConversion if non-empty: encoding of the method to be used to
	 *        convert SHAPE reactivity data to pairing probabilities
	 */
	AccessibilityConstraint( const RnaSequence & seq
							, const std::string& dotBracket
							, const size_t maxBpSpan
							, const std::string & shapeFile
							, const std::string & shapeMethod
							, const std::string & shapeConversion
							);

	virtual ~AccessibilityConstraint();

	/**
	 * Checks whether or not a sequence position is marked as blocked or not
	 * @param i the position of interest
	 * @return true if position i is marked blocked
	 */
	bool
	isMarkedBlocked( const size_t i ) const;

	/**
	 * Checks whether or not a sequence position is marked as accessible or not
	 * @param i the position of interest
	 * @return true if position i is marked accessible
	 */
	bool
	isMarkedAccessible( const size_t i ) const;

	/**
	 * Checks whether or not a sequence position is marked as intramolecularly
	 * paired or not
	 * @param i the position of interest
	 * @return true if position i is marked accessible
	 */
	bool
	isMarkedPaired( const size_t i ) const;

	/**
	 * Checks whether or not a sequence position is not constrained
	 * @param i the position of interest
	 * @return true if position i is not constrained
	 */
	bool
	isUnconstrained( const size_t i ) const;

	/**
	 * Checks whether or not a position is available for interaction
	 * @param i the position of interest
	 * @return true if the position @p i is available interaction; false otherwise
	 */
	bool
	isAccessible( const size_t i ) const;

	/**
	 * Checks whether or not any accessibility constraints (base pairs, blocked,
	 * accessible, etc.) are given
	 * @return true if no structural constraints are present; false otherwise
	 */
	bool
	isEmpty() const;


	/**
	 * Provides the VRNA dot-bracket notation of the constraint for position i
	 * @param i the position of interest
	 * @return the VRNA conform constraint encoding for position i
	 */
	char
	getVrnaDotBracket( const size_t i ) const;

	/**
	 * Provides the maximal base pair span to be considered for accessibility
	 * computation.
	 * @return the maximal base pair span for accessibility computation
	 */
	size_t
	getMaxBpSpan() const;

	/**
	 * Get filename from where to read SHAPE reactivity data.
	 * @return the filename or an empty string if no SHAPE data available
	 */
	const std::string &
	getShapeFile() const;

	/**
	 * Get method how to convert SHAPE reactivity data to pseudo energies.
	 * @return the encoding or an empty string if no SHAPE data available
	 */
	const std::string &
	getShapeMethod() const;

	/**
	 * Get method how to convert SHAPE reactivity data to pairing probabilities.
	 * @return the encoding or an empty string if no SHAPE data available
	 */
	const std::string &
	getShapeConversion() const;

	/**
	 * Assignment of constraints for the same rna sequence.
	 *
	 * @param c the interaction to make this a copy of
	 */
	AccessibilityConstraint & operator= ( const AccessibilityConstraint & c );


	friend std::ostream& operator<<(std::ostream& out, const AccessibilityConstraint& c);

protected:

	//! the overall sequence length this constraint is about
	size_t length;

	//! the maximal base pair span to be used for accessibility computation
	size_t maxBpSpan;

	//! filename of SHAPE reactivity data or empty if no SHAPE data to be used
	std::string shapeFile;

	//! method for converting SHAPE reactivity data to pseudo energies
	std::string shapeMethod;

	//! method for conversion of SHAPE reactivity data to paired probabilities
	std::string shapeConversion;

	//! sorted list of ranges that are marked as blocked
	IndexRangeList blocked;

	//! sorted list of ranges that are marked as accessible
	IndexRangeList accessible;

	//! sorted list of ranges that are marked as paired (intramolecular)
	IndexRangeList paired;

protected:

	/**
	 * screens the given dot-bracket string for consecutive regions of
	 * marker characters and pushes the according regions to the storage
	 * @param dotBracket the dot-bracket encoding to screen
	 * @param marker the marker character to screen for
	 * @param storage the container to push the identified regions to
	 */
	static
	void
	screenDotBracket( const std::string& dotBracket
					, const char marker
					, IndexRangeList & storage );

};


/**
 * Prints the accessibility constraint in range-based encoding to stream
 * @param out the ostream to write to
 * @param c the AccessibilityConstraint object to add
 * @return the altered stream out
 */
std::ostream& operator<<(std::ostream& out, const AccessibilityConstraint& c);


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

inline
AccessibilityConstraint::
AccessibilityConstraint( const size_t length_
					, const size_t maxBpSpan_
					, const std::string & shapeFile_
					, const std::string & shapeMethod_
					, const std::string & shapeConversion_
					)
	:
	length(length_),
	maxBpSpan( maxBpSpan_==0 ? length : std::min(maxBpSpan_,length) ),
	shapeFile(shapeFile_),
	shapeMethod(shapeFile_.empty() ? "" : shapeMethod_),
	shapeConversion(shapeFile_.empty() ? "" : shapeConversion_),
	blocked(),
	accessible(),
	paired()
{
#if INTARNA_IN_DEBUG_MODE
	if (!shapeFile.empty()) {
		if (!boost::regex_match( shapeMethod, AccessibilityConstraint::regexShapeMethod, boost::match_perl )) {
			throw std::runtime_error("AccessibilityConstraint(shapeMethod="+shapeMethod+") does not match its encoding regular expression");
		}
		if (!boost::regex_match( shapeConversion, AccessibilityConstraint::regexShapeConversion, boost::match_perl )) {
			throw std::runtime_error("AccessibilityConstraint(shapeConversion="+shapeConversion+") does not match its encoding regular expression");
		}
	}
#endif
}

////////////////////////////////////////////////////////////////////////

inline
AccessibilityConstraint::
AccessibilityConstraint( const AccessibilityConstraint& toCopy
			, const bool reverseIndices)
	:
	length(toCopy.length)
	, maxBpSpan(toCopy.maxBpSpan)
	, shapeFile(toCopy.shapeFile)
	, shapeMethod(toCopy.shapeMethod)
	, shapeConversion(toCopy.shapeConversion)
	, blocked(toCopy.blocked)
	, accessible(toCopy.accessible)
	, paired(toCopy.paired)
{

	if (reverseIndices) {

		// reverse blocked
		blocked.reverse(length);

		// reverse accessible
		accessible.reverse(length);

		// reverse accessible
		paired.reverse(length);

	}
}

////////////////////////////////////////////////////////////////////////

inline
AccessibilityConstraint::~AccessibilityConstraint()
{
}

////////////////////////////////////////////////////////////////////////

inline
bool
AccessibilityConstraint::
isMarkedBlocked(const size_t i) const
{
	return blocked.covers(i);
}

////////////////////////////////////////////////////////////////////////

inline
bool
AccessibilityConstraint::
isMarkedAccessible(const size_t i) const
{
	return accessible.covers(i);
}

////////////////////////////////////////////////////////////////////////

inline
bool
AccessibilityConstraint::
isMarkedPaired(const size_t i) const
{
	return paired.covers(i);
}

////////////////////////////////////////////////////////////////////////

inline
bool
AccessibilityConstraint::
isUnconstrained( const size_t i ) const
{
	return isEmpty()
			|| !(isMarkedAccessible(i) || isMarkedBlocked(i) || isMarkedPaired(i));
}

////////////////////////////////////////////////////////////////////////

inline
bool
AccessibilityConstraint::
isAccessible( const size_t i ) const
{
	return isEmpty()
			|| !(isMarkedBlocked(i) || isMarkedPaired(i));
}

////////////////////////////////////////////////////////////////////////

inline
bool
AccessibilityConstraint::
isEmpty() const
{
	// check if any constrained regions given
	return (accessible.size() + blocked.size() + paired.size() + shapeFile.size()) == 0;
}

////////////////////////////////////////////////////////////////////////

inline
char
AccessibilityConstraint::
getVrnaDotBracket(const size_t i) const
{
	// check if to be accessible or blocked (==unstructured)
	if (isMarkedAccessible(i) || isMarkedBlocked(i)) {
		return 'x';
	}

	// check if intramolecularly paired
	if (isMarkedPaired(i)) {
		return '|';
	}

	return '.';
}

////////////////////////////////////////////////////////////////////////

inline
size_t
AccessibilityConstraint::
getMaxBpSpan() const
{
	return maxBpSpan;
}

////////////////////////////////////////////////////////////////////////

inline
const std::string &
AccessibilityConstraint::
getShapeFile() const
{
	return shapeFile;
}

////////////////////////////////////////////////////////////////////////

inline
const std::string &
AccessibilityConstraint::
getShapeMethod() const
{
	return shapeMethod;
}

////////////////////////////////////////////////////////////////////////

inline
const std::string &
AccessibilityConstraint::
getShapeConversion() const
{
	return shapeConversion;
}

////////////////////////////////////////////////////////////////////////

inline
AccessibilityConstraint &
AccessibilityConstraint::
operator= ( const AccessibilityConstraint & c )
{
	// copy data
	length = c.length;
	maxBpSpan = c.maxBpSpan;
	shapeFile = c.shapeFile;
	shapeMethod = c.shapeMethod;
	shapeConversion = c.shapeConversion;
	blocked = c.blocked;
	accessible = c.accessible;
	paired = c.paired;

	return *this;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

inline
std::ostream& operator<<(std::ostream& out, const AccessibilityConstraint& c)
{
	bool alreadyPrinted = false;
	// print ranges if non-empty
	if (!c.accessible.empty()) {
		out <<(alreadyPrinted?",":"") <<AccessibilityConstraint::dotBracket_accessible <<':' <<c.accessible;
		alreadyPrinted = true;
	}
	// print ranges if non-empty
	if (!c.blocked.empty()) {
		out <<(alreadyPrinted?",":"") <<AccessibilityConstraint::dotBracket_blocked <<':' <<c.blocked;
		alreadyPrinted = true;
	}
	// print ranges if non-empty
	if (!c.paired.empty()) {
		out <<(alreadyPrinted?",":"") <<AccessibilityConstraint::dotBracket_paired <<':' <<c.paired;
		alreadyPrinted = true;
	}
	// print shape file if non-empty
	if (!c.shapeFile.empty()) {
		out <<(alreadyPrinted?",":"") <<"shapeFile" <<':' <<c.shapeFile;
		alreadyPrinted = true;
	}
	return out;
}

////////////////////////////////////////////////////////////////////////

} // namespace

#endif /* ACCESSIBILITYCONSTRAINT_H_ */
