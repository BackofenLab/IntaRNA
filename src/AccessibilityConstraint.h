
#ifndef ACCESSIBILITYCONSTRAINT_H_
#define ACCESSIBILITYCONSTRAINT_H_

#include "IndexRangeList.h"
#include "RnaSequence.h"

#include <utility>
#include <vector>

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

	//! the marker for blocked positions in dot-bracket notation
	static const char dotBracket_blocked;
	//! the marker for accessible positions in dot-bracket notation
	static const char dotBracket_accessible;

	//! the alphabet to encode accessibility constraints in dot-bracket notation
	static const std::string dotBracketAlphabet;


public:

	/**
	 * Empty constraint construction
	 */
	AccessibilityConstraint( const size_t length );

	/**
	 * Copy construction
	 * @param toCopy the constraint to copy
	 * @param revereseIndices whether or not to reverse indexing (eg to be used
	 *        for ReverseAccessibility)
	 */
	AccessibilityConstraint( const AccessibilityConstraint& toCopy
						, const bool reverseIndices = false );

	/**
	 * Constraint construction based on VRNA-like dot-bracket encoding
	 * @param dotBracket the constraint encoding in VRNA-like dot-bracket encoding
	 */
	AccessibilityConstraint( const std::string& dotBracket );

	virtual ~AccessibilityConstraint();

	/**
	 * Checks whether or not a sequence position is marked as blocked or not
	 * @return true if position i is marked blocked
	 */
	bool
	isBlocked( const size_t i ) const;

	/**
	 * Checks whether or not a sequence position is marked as accessible or not
	 * @return true if position i is marked accessible
	 */
	bool
	isAccessible( const size_t i ) const;

	/**
	 * Checks whether or not a sequence position is not constrained
	 * @return true if position i is not constrained
	 */
	bool
	isUnconstrained( const size_t i ) const;

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
	 * Assignment of constraints for the same rna sequence.
	 *
	 * @param c the interaction to make this a copy of
	 */
	AccessibilityConstraint & operator= ( const AccessibilityConstraint & c );


protected:

	//! the overall sequence length this constraint is about
	size_t length;

	//! sorted list of ranges that are marked as blocked
	IndexRangeList blocked;

	//! sorted list of ranges that are marked as accessible
	IndexRangeList accessible;

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

#endif /* ACCESSIBILITYCONSTRAINT_H_ */
