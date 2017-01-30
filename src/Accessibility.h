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
#include "AccessibilityConstraint.h"

#include <stdexcept>

/**
 * Abstract interface that represents accessibility data for a given RNA
 * sequence.
 *
 * TODO : init function to trigger accessibility computation for a certain region
 *
 * @author Martin Mann 2014
 */
class Accessibility {

public:

	//! upper bound for all ED return values
	const static E_type ED_UPPER_BOUND;

public:

	/**
	 * Construction
	 * @param sequence the sequence the accessibility data belongs to
	 * @param maxLength the maximal length of accessible regions (>0) to be
	 *          considered. 0 defaults to the full sequence's length, otherwise
	 *          is is internally set to min(maxLength,seq.length).
	 * @param accConstr optional accessibility constraint
	 */
	Accessibility( const RnaSequence& sequence
					, const size_t maxLength
					, const AccessibilityConstraint * const accConstr
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
	 * @return the ES value for [i,j] or E_INF if no intramolecular
	 *         structure can be formed
	 */
	virtual
	E_type
	getES( const size_t i, const size_t j ) const = 0;


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
	size_t
	getMaxLength() const;

	/**
	 * Access to the globally enforced accessibility constraint. Here '.'
	 * denotes unconstrained positions and 'x' positions that have to be
	 * unstructured. Regions covering constrained positions will result in
	 * ED_UPPER_BOUND accessibility values.
	 * @return the global accessibility constraint applied
	 */
	virtual
	const AccessibilityConstraint&
	getAccConstraint() const;

	/**
	 * Writes the ED values as unpaired probabilities in RNAplfold style to
	 * stream.
	 *
	 * @param out the output stream to write to
	 * @param RT the scaled temperature value to be used for conversion of
	 *        ED to Pu :  Pu = exp( -ED/RT )
	 */
	void
	writeRNAplfold_Pu_text( std::ostream& out, const E_type RT ) const;

	/**
	 * Writes the ED values in RNAplfold style to stream.
	 *
	 * @param out the output stream to write to
	 */
	void
	writeRNAplfold_ED_text( std::ostream& out ) const;

	/**
	 * Prints the accessibility values to stream as upper triangular matrix
	 * @param out the ostream to write to
	 * @param acc the Accessibility object to add
	 * @return the altered stream out
	 */
	friend std::ostream& operator<<(std::ostream& out, const Accessibility& acc);


protected:

	//! the RNA sequence the accessibilities correspond to
	const RnaSequence & seq;

	//! the maximal length of an unpaired regions to be considered
	const size_t maxLength;

	//! accessibility constraint
	AccessibilityConstraint accConstraint;

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

	/**
	 * Writes the ED values as unpaired probabilities in RNAplfold style to
	 * stream.
	 *
	 * @param out the output stream to write to
	 * @param RT the scaled temperature value to be used for conversion of
	 *        ED to Pu :  Pu = exp( -ED/RT )
	 * @param writeProbs (true) write unpaired probabilities; (false) write ED
	 */
	void
	writeRNAplfold_text( std::ostream& out, const E_type RT, const bool writeProbs ) const;

};



/////////////////////////////////////////////////////////////////////////////

inline
Accessibility::Accessibility( const RnaSequence& seq
							, const size_t maxLength
							, const AccessibilityConstraint * const accConstraint_ )
 :
	seq(seq)
	// set maxLength to appropriate value
	, maxLength( maxLength==0 ? seq.size() : std::min(maxLength,seq.size()) )
	, accConstraint( seq.size() )
{
	// set constraint if needed
	if (accConstraint_ != NULL) {
		accConstraint = *accConstraint_;
	}
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
#if IN_DEBUG_MODE
	if (from > to || to >= getSequence().size()) {
		throw std::runtime_error("Accessibility::checkIndices : region ["+toString(from)+","+toString(to)+"] does not fulfill 0 <= from <= to < seq.length");
	}
#endif
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
size_t
Accessibility::
getMaxLength() const
{
	return maxLength;
}

/////////////////////////////////////////////////////////////////////////////

inline
const AccessibilityConstraint&
Accessibility::
getAccConstraint() const
{
	return accConstraint;
}

/////////////////////////////////////////////////////////////////////////////

inline
void
Accessibility::
writeRNAplfold_ED_text( std::ostream& out ) const
{
	writeRNAplfold_text( out, 1.0, false );
}

/////////////////////////////////////////////////////////////////////////////

inline
void
Accessibility::
writeRNAplfold_Pu_text( std::ostream& out, const E_type RT ) const
{
	writeRNAplfold_text( out, RT, true );
}

/////////////////////////////////////////////////////////////////////////////


#endif /* ACCESSIBILITY_H_ */
