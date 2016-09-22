/*
 * ReverseAccessibility.h
 *
 *  Created on: 30.06.2014
 *      Author: Mmann
 */

#ifndef REVERSEACCESSIBILITY_H_
#define REVERSEACCESSIBILITY_H_

#include "Accessibility.h"

/**
 * Defines an accessibility meta object with reversed index access to sequence
 * and accessibility.
 */
class ReverseAccessibility: public Accessibility {
public:

	/**
	 * Construction based on a given accessibility object to represent reversed
	 *
	 * @param origAcc Access to the accessibility object to reverse
	 */
	ReverseAccessibility( Accessibility & origAcc );

	/**
	 * destruction
	 */
	virtual ~ReverseAccessibility();

	/**
	 * Access to the accessibility energy value (according to reversed indices).
	 *
	 * @param from the start index of the regions (from <= to)
	 * @param to the end index of the regions (to <= seq.length())
	 *
	 * @return 0  if (j-1+1) <= maxLength or ED_UPPER_BOUND otherwise
	 */
	virtual
	E_type
	getED( const size_t from, const size_t to ) const;


	/**
	 * Access to the RnaSequence this accessibility values are accounting for.
	 * @return the underlying sequence for this accessibility object.
	 */
	virtual
	const RnaSequence &
	getSequence() const;


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

	/**
	 * Access to the original not-reversed accessibility object.
	 * @return the original (index reversed) accessibility object.
	 */
	virtual
	const Accessibility &
	getAccessibilityOrigin() const;


protected:


	//! Access to the accessibility object to reverse
	Accessibility & origAcc;

	//! reversed sequence object
	RnaSequence seqReversed;

	//! reversed accessibility constraint object
	std::string accConstrReversed;


	/**
	 * Computes a reversed string representation of a sequence
	 * @param seq the sequence to reverse
	 * @return the reversed string representation
	 */
	static
	std::string
	getReversedString( const RnaSequence& seq );

	/**
	 * Computes a reversed string representation of a given string
	 * @param str the string to reverse
	 * @return the reversed string
	 */
	static
	std::string
	getReversedString( const std::string& str );


};

#endif /* REVERSEACCESSIBILITY_H_ */
