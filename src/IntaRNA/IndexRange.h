
#ifndef INTARNA_INDEXRANGE_H_
#define INTARNA_INDEXRANGE_H_


#include "IntaRNA/general.h"
#include "IntaRNA/RnaSequence.h"

#include <stdexcept>

#include <boost/regex.hpp>

namespace IntaRNA {

/**
 * Defines an index region with functionality similar to pair
 *
 * @author Martin Mann
 */
class IndexRange {

public:

	//! placeholder to define the whole range (without explicit naming the last
	//! index) is defined
	static const size_t LAST_INDEX;

	//! placeholder for not-defined values
	static const size_t NA_INDEX;

	//! the start of the index range
	size_t from;

	//! the end of the index range
	size_t to;

	//! regular expression that matches valid IndexRange string encodings
	static constexpr const char* regexString = R"((0|-?[123456789]\d*)-(0|-?[123456789]\d*))";

	//! regular expression that matches valid IndexRange string encodings
	static const boost::regex regex;

public:

	/**
	 * Creates a range
	 * @param from the start index (default 0)
	 * @param to the end index (default NA_INDEX)
	 */
	IndexRange(const size_t from = 0,
			const size_t to = NA_INDEX)
		: from(from), to(to)
	{
	}

	/**
	 * Creates a range from a string encoding
	 * @param stringEncoding string encoding of the range as produced by the
	 *  ostream operator
	 * @param seq if not NULL, the RnaSequence to be used to shift the indices
	 *            of the input string from in/output positions to internal
	 *            indices
	 *
	 */
	IndexRange(const std::string & stringEncoding, const RnaSequence * seq )
		: from(0), to(NA_INDEX)
	{
		fromString(stringEncoding,seq);
	}

	/**
	 * destruction
	 */
	virtual ~IndexRange() {}

	/**
	 * Checks whether or not the range encoding is ascending.
	 * @return from <= to
	 */
	bool isAscending() const
	{
		return from <= to;
	}

	/**
	 * Checks whether or not the range encoding is descending.
	 * @return from >= to
	 */
	bool isDescending() const
	{
		return from >= to;
	}

	/**
	 * Adds the given shift to the index range but ensure that the minimal value
	 * is 0.
	 * @param r the range to be added to
	 * @param shift the shift to be added to r
	 * @return an altered range or (NA_INDEX,NA_INDEX) if the range falls
	 *    completely below zero or changing from to the lower bound of 0.
	 */
	IndexRange operator + ( const int shift ) const {
		if (shift == 0) {
			return *this;
		}
		if (shift > 0) {
			return IndexRange(from+shift, to+shift);
		}
		if (to < std::abs(shift)) {
			return IndexRange( NA_INDEX, NA_INDEX);
		}
		return IndexRange( from - std::min(from, (size_t)std::abs(shift)), to + shift );
	}

	/**
	 * Substracts the given shift to the index range but ensure that the minimal
	 * value is 0.
	 * @param r the range to be altered
	 * @param shift the shift to be substracted from r
	 * @return an altered range or (NA_INDEX,NA_INDEX) if the range falls
	 *    completely below zero or changing from to the lower bound of 0.
	 */
	IndexRange operator - ( const int shift ) const {
		// forward implementation using inverted shift
		return this->operator+(-shift);
	}

	/**
	 * Checks whether or not the start of this range preceeds the given range.
	 * @param r the range to compare to
	 * @return (from<r.from) || (from==r.from && to<r.to)
	 */
	const bool operator < ( const IndexRange &r ) const {
		return ( from < r.from || (from==r.from && to<r.to) );
	}

	/**
	 * Checks whether or not two ranges are equivalent
	 * @param r the range to compare to
	 * @return ( from == r.from && to == r.to )
	 */
	const bool operator == ( const IndexRange &r ) const {
		return ( from == r.from && to == r.to );
	}

	/**
	 * Checks whether or not two ranges are different
	 * @param r the range to compare to
	 * @return !( this == r)
	 */
	const bool operator != ( const IndexRange &r ) const {
		return !( this->operator ==(r) );
	}


	/**
	 * Prints the range's boundaries to stream
	 * @param out the ostream to write to
	 * @param range the IndexRange object to add
	 * @return the altered stream out
	 */
	friend std::ostream& operator<<(std::ostream& out, const IndexRange& range)
	{
		// has to be in accordance to this.regex
		return (out <<range.from<<"-"<<range.to);
	}

	/**
	 * updates the range data from a valid string encoding (matching regex)
	 * @param stringEncoding the interval string encoding
	 * @throws std::runtime_error if stringEncoding does not match regex
	 */
	void
	fromString( const std::string & stringEncoding, const RnaSequence * seq )
	{
		if( ! boost::regex_match(stringEncoding, regex, boost::match_perl) ) {
			throw std::runtime_error("IndexRange::fromString("+stringEncoding+") uses no valid index range string encoding");
		}
		// find split position assuming second index to be negative
		size_t splitPos = stringEncoding.find("--");
		if (splitPos == std::string::npos) {
			// last index should not be negative, thus, split at last '-'
			splitPos = stringEncoding.find_last_of('-');
		}
		// parse interval boundaries
		if (seq == NULL) {
			from = boost::lexical_cast<size_t>(stringEncoding.substr(0,splitPos));
			to = boost::lexical_cast<size_t>(stringEncoding.substr(splitPos+1));
		} else {
			// correct for shifted indexing within the input sequence
			from = seq->getIndex(boost::lexical_cast<long>(stringEncoding.substr(0,splitPos)));
			to = seq->getIndex(boost::lexical_cast<long>(stringEncoding.substr(splitPos+1)));
		}
	}
	
	/**
	 * Computes an overlapping window decomposition of this IndexRange
	 * @param windowWidth the width of a window
	 * @param windowOverlap the amount of overlap between two windows
	 * @return a vector of overlapping IndexRanges covering this IndexRange
	 */
	std::vector<IndexRange>
	overlappingWindows(const size_t windowWidth, const size_t windowOverlap) const;
	
	
};

} // namespace

#endif /* INDEXRANGE_H_ */
