
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
	 * Provides the ensemble energy (ES) (according to reversed indices).
	 *
	 * @param i the start of the structured region
	 * @param j the end of the structured region
	 * @return the ES value for [i,j] or E_INF if no intramolecular
	 *         structure can be formed
	 */
	virtual
	E_type
	getES( const size_t i, const size_t j ) const;


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
	const AccessibilityConstraint&
	getAccConstraint() const;

	/**
	 * Access to the original not-reversed accessibility object.
	 * @return the original (index reversed) accessibility object.
	 */
	virtual
	const Accessibility &
	getAccessibilityOrigin() const;


	/**
	 * Provides the reverse index of a given sequence position.
	 * @param i the index of interest
	 * @return the reverse index, i.e. (seq.size()-i-1)
	 */
	size_t
	getReversedIndex( const size_t i ) const;

	/**
	 * Provides the reverse index range of a given sequence position.
	 * @param i the index range of interest
	 * @return the reversed index range, i.e. for index i : (seq.size()-i-1)
	 */
	IndexRange
	getReversedIndexRange( const IndexRange & r ) const;


protected:


	//! Access to the accessibility object to reverse
	Accessibility & origAcc;

	//! reversed sequence object
	RnaSequence seqReversed;

	//! reversed accessibility constraint object
	AccessibilityConstraint accConstrReversed;


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


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

inline
ReverseAccessibility::ReverseAccessibility( Accessibility & origAcc )
 :
	Accessibility( origAcc.getSequence(), origAcc.getMaxLength(), &origAcc.getAccConstraint() )
	, origAcc(origAcc)
	, seqReversed( seq.getId(), getReversedString(seq) )
	, accConstrReversed( origAcc.getAccConstraint(), true )
{
}

////////////////////////////////////////////////////////////////////////////

inline
ReverseAccessibility::~ReverseAccessibility() {
}

////////////////////////////////////////////////////////////////////////////

inline
const RnaSequence &
ReverseAccessibility::
getSequence() const
{
	// access reversed sequence
	return seqReversed;
}


////////////////////////////////////////////////////////////////////////////

inline
const AccessibilityConstraint&
ReverseAccessibility::
getAccConstraint() const
{
	// access reversed accessibility constraint
	return accConstrReversed;
}


////////////////////////////////////////////////////////////////////////////


inline
E_type
ReverseAccessibility::
getED( const size_t from, const size_t to ) const
{
	// check indices
	checkIndices(from,to);
	// reversed ED access
	return origAcc.getED( seq.size()-to-1, seq.size()-from-1 );
}


////////////////////////////////////////////////////////////////////////////


inline
E_type
ReverseAccessibility::
getES( const size_t from, const size_t to ) const
{
	// check indices
	checkIndices(from,to);
	// reversed ES access
	return origAcc.getES( seq.size()-to-1, seq.size()-from-1 );
}


////////////////////////////////////////////////////////////////////////////

inline
const Accessibility &
ReverseAccessibility::
getAccessibilityOrigin() const
{
	return origAcc;
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
ReverseAccessibility::
getReversedIndex( const size_t i ) const
{
	// check indices
	checkIndices(i,i);
	// compute reverse index
	return this->seq.size() -i -1;
}

////////////////////////////////////////////////////////////////////////////

inline
IndexRange
ReverseAccessibility::
getReversedIndexRange( const IndexRange & r ) const
{
	// reverse indices and their order to keep them ascending/descending
	return IndexRange( getReversedIndex(r.to), getReversedIndex(r.from) );
}

////////////////////////////////////////////////////////////////////////////



#endif /* REVERSEACCESSIBILITY_H_ */
