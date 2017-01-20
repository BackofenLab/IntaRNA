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
	 * @param accConstr optional accessibility constraint
	 */
	AccessibilityDisabled( const RnaSequence& sequence
							, const size_t maxLength = 0
							, const AccessibilityConstraint * const accConstr = NULL
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


	/**
	 * Provides the ensemble energy (ES) of all intramolecular substructures
	 * that can be formed within a given region of the sequence under the
	 * assumption that the region is part of an (intermolecular) multiloop,
	 * i.e. at least one base pair is formed by each substructure, given by
	 *
	 *   ES(i,j) = -RT * log( Qm[i,j] )
	 *
	 * where Qm denotes the according partition function.
	 *
	 * If no structure can be formed within the region (Qm==0), E_INF is returned.
	 *
	 * @param i the start of the structured region
	 * @param j the end of the structured region
	 * @return E_INF
	 */
	virtual
	E_type
	getES( const size_t i, const size_t j ) const;



};



///////////////////////////////////////////////////////////////////////////////

inline
AccessibilityDisabled::AccessibilityDisabled(const RnaSequence& seq
				, const size_t maxLength
				, const AccessibilityConstraint * const accConstr)
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
			if (!getAccConstraint().isAccessible(i)) {
				// position covers a blocked position --> omit accessibility
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

inline
E_type
AccessibilityDisabled::
getES( const size_t from, const size_t to ) const
{
	throw std::runtime_error("accessibility computation was disabled, which also disables ES computation..");
	return ED_UPPER_BOUND;
}

///////////////////////////////////////////////////////////////////////////////


#endif /* ACCESSIBILITYDISABLED_H_ */
