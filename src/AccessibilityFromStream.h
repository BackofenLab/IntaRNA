
#ifndef ACCESSIBILITYFROMSTREAM_H_
#define ACCESSIBILITYFROMSTREAM_H_

#include "Accessibility.h"

#include <iostream>

#include <boost/numeric/ublas/banded.hpp>

/**
 * Reads accessibility information from a data stream, e.g. from file or STDIN
 *
 */
class AccessibilityFromStream: public Accessibility
{
public:

	enum InStreamType {
		Pu_RNAplfold_Text
	};

public:

	/**
	 * construction
	 * @param sequence the sequence the accessibility data is about
	 * @param maxLength the maximal length of accessible regions (>0) to be
	 *          considered. 0 defaults to the full sequence's length, otherwise
	 *          is is internally set to min(maxLength,seq.length).
	 * @param accConstraint if not NULL, accessibility constraint that enforces some regions
	 *        to be unstructured both in sequence and interaction
	 * @param inStream the input stream to read the accessibility data from
	 * @param inStreamType inStream data type to be expected
	 * @param RT the RT constant to be used to transform the probabilities to
	 *        ED values
	 */
	AccessibilityFromStream(
			const RnaSequence& sequence
			, const size_t maxLength
			, const AccessibilityConstraint * const accConstraint
			, std::istream & inStream
			, const InStreamType inStreamType
			, const E_type RT
			);


	/**
	 * destruction
	 */
	virtual ~AccessibilityFromStream();


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
	getED( const size_t from, const size_t to ) const;

	/**
	 * Not available for this subclass implementation.
	 *
	 * @param i the start of the structured region
	 * @param j the end of the structured region
	 * @return E_INF since not implemented
	 *
	 * @throws std::runtime_error not implemented
	 */
	virtual
	E_type
	getES( const size_t i, const size_t j ) const;


protected:

	//! type for the ED value matrix (upper triangular matrix banded by maxLength)
	typedef boost::numeric::ublas::banded_matrix<E_type> EdMatrix;

	//! the ED values for the given sequence
	EdMatrix edValues;


	/**
	 * Parses a VRNA unpaired probability file and fills the ED data
	 *
	 * @param inStream the stream to read the probabilities from
	 * @param RT the RT constant to be used to transform the probabilities to
	 *        ED values
	 */
	void
	parsePu_RNAplfold_Text( std::istream & inStream, const E_type RT );


};


/////////////////////////////////////////////////////////////////////////

inline
E_type
AccessibilityFromStream::
getES( const size_t i, const size_t j ) const
{
	NOTIMPLEMENTED("AccessibilityFromStream::getES() not supported");
	return E_INF;
}

/////////////////////////////////////////////////////////////////////////

inline
E_type
AccessibilityFromStream::
getED( const size_t from, const size_t to ) const
{
	// input range check
	checkIndices(from,to);

	if ((to-from+1) <= getMaxLength()) {
		// check for constrained positions within region
		if (!getAccConstraint().isAccessible(from, to)) {
			// position covers a blocked position --> omit accessibility
			return ED_UPPER_BOUND;
		}
		// return according ED value from the precomputed matrix
		return edValues (from,to);
	} else {
		// region length exceeds maximally allowed length -> no value
		return ED_UPPER_BOUND;
	}
}

/////////////////////////////////////////////////////////////////////////



#endif /* ACCESSIBILITYFROMSTREAM_H_ */
