
#ifndef INTARNA_ACCESSIBILITYFROMSTREAM_H_
#define INTARNA_ACCESSIBILITYFROMSTREAM_H_

#include "IntaRNA/Accessibility.h"

#include <iostream>

#include <boost/numeric/ublas/banded.hpp>

namespace IntaRNA {

/**
 * Reads accessibility information from a data stream, e.g. from file or STDIN
 *
 */
class AccessibilityFromStream: public Accessibility
{
public:

	enum InStreamType {
		Pu_RNAplfold_Text //! Pu values in RNAplfold text format
		, ED_RNAplfold_Text //!< ED values in RNAplfold text Pu format
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
			, const Z_type RT
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
	 * Access to the maximal length of accessible regions (>0) to be considered.
	 *
	 * Here, it returns the minimum of the originally targeted interaction range
	 * and the from the data parsed maximal window size.
	 *
	 * @return the maximal length of accessible regions considered
	 */
	virtual
	size_t
	getMaxLength() const;


protected:

	//! type for the ED value matrix (upper triangular matrix banded by maxLength)
	typedef boost::numeric::ublas::banded_matrix<E_type> EdMatrix;

	//! the ED values for the given sequence
	EdMatrix edValues;

	//! maximal available window size
	size_t availMaxLength;

	/**
	 * Parses a VRNA unpaired probability file and fills the ED data
	 *
	 * @param inStream the stream to read the probabilities from
	 * @param RT the RT constant to be used to transform the probabilities to
	 *        ED values
	 */
	void
	parsePu_RNAplfold_text( std::istream & inStream, const Z_type RT );


	/**
	 * Parses ED values from a VRNA unpaired probability file styled stream
	 *
	 * @param inStream the stream to read the ED values from
	 */
	void
	parseED_RNAplfold_text( std::istream & inStream );

	/**
	 * Parses ED values from a VRNA unpaired probability file styled stream
	 *
	 * @param inStream the stream to read the ED values from
	 * @param RT the RT constant to be used to transform the probabilities to
	 *        ED values
	 * @param parseProbs whether or not to expect unpaired probabilities (true)
	 *        or ED values within the file
	 */
	void
	parseRNAplfold_text( std::istream & inStream
						, const Z_type RT
						, const bool parseProbs );


};

/////////////////////////////////////////////////////////////////////////

inline
E_type
AccessibilityFromStream::
getED( const size_t from, const size_t to ) const
{
	// input range check
	checkIndices(from,to);

	if ((to-from+1) <= getMaxLength()) {
		// check for constrained end positions
		if (!getAccConstraint().isAccessible(from) || !getAccConstraint().isAccessible(to)) {
			// end position blocked --> omit accessibility
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

inline
size_t
AccessibilityFromStream::
getMaxLength() const
{
	return availMaxLength;
}

/////////////////////////////////////////////////////////////////////////

inline
void
AccessibilityFromStream::
parsePu_RNAplfold_text( std::istream & inStream, const Z_type RT )
{
	parseRNAplfold_text( inStream, RT, true );
}

/////////////////////////////////////////////////////////////////////////

inline
void
AccessibilityFromStream::
parseED_RNAplfold_text( std::istream & inStream )
{
	parseRNAplfold_text( inStream, 1.0, false );
}

/////////////////////////////////////////////////////////////////////////

} // namespace

#endif /* ACCESSIBILITYFROMSTREAM_H_ */
