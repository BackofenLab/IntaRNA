/*
 * Accessibility.h
 *
 *  Created on: 25.06.2014
 *      Author: Mmann
 */

#ifndef ACCESSIBILITY_H_
#define ACCESSIBILITY_H_


#include "general.h"
#include "RnaSequence.h"
#include <stdexcept>

/**
 * Abstract interface that represents accessibility data for a given RNA
 * sequence.
 *
 * @author Martin Mann 2014
 */
class Accessibility {

public:

	//! upper bound for all ED return values
	const static E_type ED_UPPER_BOUND;

	//! allowed alphabet for accessibility constraint encoding
	const static std::string AccessibilityConstraintAlphabet;

public:

	/**
	 * Construction
	 * @param sequence the sequence the accessibility data belongs to
	 * @param maxLength the maximal length of accessible regions (>0) to be
	 *          considered. 0 defaults to the sequence's length, otherwise
	 *          is is internally set to min(maxLength,seq.length).
	 * @param accConstr accessibility constraint: all positions marked not as
	 * 			unconstrained (.) are ommited from interaction prediction, thus
	 * 			resulting in ED_UPPER_BOUND. Empty string input is equivalent
	 * 			to fully unconstrained accessibility computation.
	 */
	Accessibility( const RnaSequence& sequence
					, const size_t maxLength
					, const std::string& accConstr
				);

	/**
	 * destruction
	 */
	virtual ~Accessibility();

	/**
	 * Returns the accessibility energy value for the given range in the
	 * sequence, i.e. the energy difference (ED) to make the region accessible.
	 *
	 * @param from the start index of the regions (from <= to)
	 * @param to the end index of the regions (to < seq.length)
	 *
	 * @return the ED value if (j-1+1) <= maxLength or ED_UPPER_BOUND otherwise
	 *
	 * @throw std::runtime_error in case it does not hold 0 <= from <= to < seq.length
	 */
	virtual
	E_type
	getED( const size_t from, const size_t to ) const = 0;


	/**
	 * Access to the RnaSequence this accessibility values are accounting for.
	 * @return the underlying sequence for this accessibility object.
	 */
	virtual
	const RnaSequence &
	getSequence() const;

	/**
	 * Access to the maximal length of accessible regions (>0) to be considered.
	 * @return the maximal length of accessible regions considered
	 */
	virtual
	const size_t
	getMaxLength() const;

	/**
	 * Access to the globally enforced accessibility constraint. Here '.'
	 * denotes unconstrained positions and 'x' positions that have to be
	 * unstructured. Regions covering constrained positions will result in
	 * ED_UPPER_BOUND accessibility values.
	 * @return the global accessibility constraint applied
	 */
	virtual
	const std::string&
	getAccConstraint() const;

protected:

	//! the RNA sequence the accessibilities correspond to
	const RnaSequence & seq;

	//! the maximal length of a unpaired regions to be considered
	const size_t maxLength;

	//! accessibility constraint: all positions marked not as unconstrained (.)
	//! are ommited from interaction prediction, thus resulting in
	//! ED_UPPER_BOUND
	std::string accConstraint;

	/**
	 * Checks the given indices to be in the range 0 <= from <= to < seq.length
	 * and throws a std::runtime_error if the constraint is not met.
	 * @param from the start index of the regions
	 * @param to the end index of the regions
	 *
	 * @throw std::runtime_error in case it does not hold 0 <= from <= to < seq.length
	 */
	virtual
	void
	checkIndices( const size_t from, const size_t to ) const;

};



/////////////////////////////////////////////////////////////////////////////

inline
Accessibility::Accessibility( const RnaSequence& seq
							, const size_t maxLength
							, const std::string& accConstraint_ )
 :
	seq(seq)
	, maxLength( maxLength==0 ? seq.size() : std::min(maxLength,seq.size()) )
	, accConstraint(accConstraint_)
{
	// fill accessibility constraint if not initialized
	if (accConstraint.empty()) {
		accConstraint = std::string(seq.size(), '.');
	}

#ifdef NDEBUG           /* required by ANSI standard */
	// no check
#else
	// debug checks
	if(accConstraint.size() != getSequence().size()) {
		throw std::runtime_error("Accessibility : sequence and accConstr differ in length");
	}
	if (accConstraint.find_first_not_of(AccessibilityConstraintAlphabet)!=std::string::npos) {
		throw std::runtime_error("Accessibility : accConstr '"+accConstraint+"' contains unsupported characters");
	}
#endif
}

/////////////////////////////////////////////////////////////////////////////

inline
Accessibility::~Accessibility()
{
}

/////////////////////////////////////////////////////////////////////////////

inline
void
Accessibility::
checkIndices( const size_t from, const size_t to ) const
{
	if (from > to || to >= getSequence().size()) {
		throw std::runtime_error("Accessibility::checkIndices : region ["+toString(from)+","+toString(to)+"] do not fulfill 0 <= from <= to < seq.length");
	}
}

/////////////////////////////////////////////////////////////////////////////

inline
const RnaSequence &
Accessibility::
getSequence() const
{
	return seq;
}

/////////////////////////////////////////////////////////////////////////////

inline
const size_t
Accessibility::
getMaxLength() const
{
	return maxLength;
}

/////////////////////////////////////////////////////////////////////////////

inline
const std::string&
Accessibility::
getAccConstraint() const
{
	return accConstraint;
}

/////////////////////////////////////////////////////////////////////////////


#endif /* ACCESSIBILITY_H_ */
