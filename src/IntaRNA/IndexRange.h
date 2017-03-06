
#ifndef INDEXRANGE_H_
#define INDEXRANGE_H_


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

	//! the start of the index range
	size_t from;

	//! the end of the index range
	size_t to;

	//! regular expression that matches valid IndexRange string encodings
	static const boost::regex regex;

public:

	/**
	 * Creates a range
	 * @param from the start index (default 0)
	 * @param to the end index (default max())
	 */
	IndexRange(const size_t from = 0,
			const size_t to = RnaSequence::lastPos)
		: from(from), to(to)
	{
	}

	/**
	 * Creates a range from a string encoding
	 * @param stringEncoding string encoding of the range as produced by the
	 *  ostream operator
	 */
	IndexRange(const std::string & stringEncoding)
		: from(0), to(RnaSequence::lastPos)
	{
		fromString(stringEncoding);
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
	fromString( const std::string & stringEncoding )
	{
		if( ! boost::regex_match(stringEncoding, regex, boost::match_perl) ) {
			throw std::runtime_error("IndexRange::fromString("+stringEncoding+") uses no valid index range string encoding");
		}
		// find split position
		const size_t splitPos = stringEncoding.find('-');
		// parse interval boundaries
		from = boost::lexical_cast<size_t>(stringEncoding.substr(0,splitPos));
		to = boost::lexical_cast<size_t>(stringEncoding.substr(splitPos+1));
	}

};

} // namespace

#endif /* INDEXRANGE_H_ */
