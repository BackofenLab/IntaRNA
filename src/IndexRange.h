
#ifndef INDEXRANGE_H_
#define INDEXRANGE_H_


#include "general.h"

#include <stdexcept>

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

public:

	/**
	 * Creates a range
	 * @param from the start index (default 0)
	 * @param to the end index (default max())
	 */
	IndexRange(const size_t from = 0,
			const size_t to = std::numeric_limits<size_t>::max())
		: from(from), to(to)
	{
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
	 * Prints the range's boundaries to stream
	 * @param out the ostream to write to
	 * @param range the IndexRange object to add
	 * @return the altered stream out
	 */
	friend std::ostream& operator<<(std::ostream& out, const IndexRange& range)
	{
		return (out <<"["<<range.from<<","<<range.to<<"]");
	}

	/**
	 * Checks whether or not the start of this range preceeds the given range.
	 * @param r the range to compare to
	 * @return (from<r.from) || (from==r.from && to<r.to)
	 */
	const bool operator < ( const IndexRange &r ) const{
		return ( from < r.from || (from==r.from && to<r.to) );
	}


};

#endif /* INDEXRANGE_H_ */
