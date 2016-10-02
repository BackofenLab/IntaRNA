/*
 * AccessibilityDisabled.h
 *
 *  Created on: 25.06.2014
 *      Author: Mmann
 */

#ifndef ACCESSIBILITYDISABLED_H_
#define ACCESSIBILITYDISABLED_H_

#include "Accessibility.h"

/**
 * Implements the Accessibility interface but disables ED value computation,
 * i.e. all ED values are set to zero.
 *
 * @author Martin Mann 2014
 */
class AccessibilityDisabled: public Accessibility {

public:

	/**
	 * Construction
	 * @param sequence the sequence the accessibility data belongs to
	 * @param maxLength the maximal length of accessible regions to be
	 *          considered. 0 defaults to the sequence's length.
	 * @param accConstr accessibility constraint: all positions marked not as
	 * 			unconstrained (.) are ommited from interaction prediction, thus
	 * 			resulting in ED_UPPER_BOUND. Empty string input is equivalent
	 * 			to fully unconstrained accessibility computation.
	 */
	AccessibilityDisabled( const RnaSequence& sequence
							, const size_t maxLength = 0
							, const std::string  * const accConstr = NULL
						);

	/**
	 * destruction
	 */
	virtual ~AccessibilityDisabled();

	/**
	 * Always returns a zero accessibility energy value.
	 *
	 * @param from the start index of the regions (from <= to)
	 * @param to the end index of the regions (to <= seq.length())
	 *
	 * @return 0  if (j-1+1) <= maxLength or ED_UPPER_BOUND otherwise
	 */
	virtual
	E_type
	getED( const size_t from, const size_t to ) const;




};



///////////////////////////////////////////////////////////////////////////////

inline
AccessibilityDisabled::AccessibilityDisabled(const RnaSequence& seq
				, const size_t maxLength
				, const std::string * const accConstr)
 :
	Accessibility(seq, maxLength, accConstr)
{
}

///////////////////////////////////////////////////////////////////////////////

inline
AccessibilityDisabled::~AccessibilityDisabled()
{
}

///////////////////////////////////////////////////////////////////////////////

inline
E_type
AccessibilityDisabled::
getED( const size_t from, const size_t to ) const
{
	// input check
	checkIndices(from,to);

	if ((to-from+1) <= getMaxLength()) {
		// check for constrained positions within region
		for (size_t i=from; i<=to; ++i) {
			if (getAccConstraint().size() > 0 && accConstraint[i] != '.') {
				// regions covers a constrained position --> omit accessibility
				return ED_UPPER_BOUND;
			}
		}
		// else: no accessibility computation done --> always zero
		return (E_type)0;
	} else {
		// region length exceeds maximally allowed length -> no value
		return ED_UPPER_BOUND;
	}
}

///////////////////////////////////////////////////////////////////////////////


#endif /* ACCESSIBILITYDISABLED_H_ */
