
#ifndef INTARNA_INDEXRANGELIST_H_
#define INTARNA_INDEXRANGELIST_H_

#include "IntaRNA/IndexRange.h"

#include <list>

#include <boost/regex.hpp>

namespace IntaRNA {

/**
 * Sorted list of ascending ranges that can be constrained to be non-overlapping
 *
 * @author Martin Mann
 *
 */
class IndexRangeList
{
protected:

	//! List of ranges
	typedef std::list< IndexRange > List;

public:

	//! regular expression string (without start/end matching)
	static const std::string regexString;

	//! regular expression that matches valid IndexRangeList string encodings
	static const boost::regex regex;

public:

	typedef List::iterator iterator;
	typedef List::const_iterator const_iterator;
	typedef List::reverse_iterator reverse_iterator;
	typedef List::const_reverse_iterator const_reverse_iterator;

public:

	/**
	 * empty construction
	 * @param allowOverlap whether or not overlapping ranges are allowed
	 */
	IndexRangeList( const bool allowOverlap = false );

	/**
	 * String encoding based construction
	 *
	 * @param stringEncoding the string encoding to be parsed
	 * @param allowOverlap whether or not overlapping ranges are allowed
	 * @param seq if not NULL, the RnaSequence to be used to shift in/output
	 *         indices to 0-based internal index positions
	 *
	 * @throws std::runtime_error if stringEncoding does not match regex
	 */
	IndexRangeList( const std::string & stringEncoding
					, const bool allowOverlap = false
					, const RnaSequence * seq = NULL);

	/**
	 * copy construction
	 * @param toCopy the list to make this a copy of
	 */
	IndexRangeList( const IndexRangeList & toCopy );


	/**
	 * destruction
	 */
	virtual ~IndexRangeList();


	/**
	 * checks whether or not a given index is covered by one of the stored
	 * index ranges
	 * @param index the index to check
	 * @return true if @p index is within one of the index ranges (including
	 *         boundaries); false otherwise
	 */
	bool covers( const size_t index ) const;

	/**
	 * checks whether or not a given index range is completely covered by one
	 * of the stored index ranges
	 * @param from the start index of the range of interest
	 * @param to the end index of the range of interest
	 * @return true if [from,to] is within one of the index ranges (including
	 *         boundaries); false otherwise
	 */
	bool covers( const size_t from, const size_t to ) const;

	/**
	 * checks whether or not a given index range is completely covered by one
	 * of the stored index ranges
	 * @param range the index range of interest
	 * @return true if range is within one of the index ranges (including
	 *         boundaries); false otherwise
	 */
	bool covers( const IndexRange & range ) const;

	/**
	 * checks whether or not a given ascending index range is overlapping with at least
	 * one of the stored index ranges
	 * @param range the range to check
	 * @return true if @p range is within one of the index ranges (including
	 *         boundaries); false otherwise
	 */
	bool overlaps( const IndexRange& range ) const;

	/**
	 * adds an ascending range to the end of the list.
	 * NOTE: if the insertion would violate range sorting, an exception is raised
	 * @param range the ascending index range to add (from should be > rbegin()->to)
	 */
	void push_back( const IndexRange& range );

	/**
	 * adds an ascending range to the according position of the sorted list.
	 * Insertion of duplicates is avoided.
	 * @param range the ascending index range to add
	 * @return iterator to the inserted (or already present) element
	 */
	iterator insert( const IndexRange& range );

	/**
	 * access to the stored range with the given index
	 * @param idx the index of the range within the list to access
	 * @return the range for index idx
	 * @throws std::runtime_error if idx >= this.size()
	 */
	IndexRange & get( const size_t idx );

	/**
	 * Constant access to the stored range with the given index
	 * @param idx the index of the range within the list to access
	 * @return the range for index idx (const access)
	 * @throws std::runtime_error if idx >= this.size()
	 */
	const IndexRange & get( const size_t idx ) const;

	/**
	 * removes an element from the list
	 * @param i the iterator pointing to the element to delete
	 */
	iterator erase( iterator i );

	/**
	 * Access to the first index range within the list
	 * @return the begin of the list
	 */
	iterator begin();

	/**
	 * Access to the end of the index range iteration
	 * @return the end of the list
	 */
	iterator end();

	/**
	 * Constant access to the first index range within the list
	 * @return the begin of the list
	 */
	const_iterator begin() const;

	/**
	 * Constant access to the end of the index range iteration
	 * @return the end of the list
	 */
	const_iterator end() const;

	/**
	 * Access to the last index range within the list
	 * @return the begin of the list
	 */
	reverse_iterator rbegin();

	/**
	 * Access to the end of the reversed index range iteration
	 * @return the end of the list
	 */
	reverse_iterator rend();

	/**
	 * Constant access to the last index range within the list
	 * @return the begin of the list
	 */
	const_reverse_iterator rbegin() const;

	/**
	 * Constant access to the end of the reversed index range iteration
	 * @return the end of the list
	 */
	const_reverse_iterator rend() const;

	/**
	 * Whether or not an index range is present
	 * @return true if no index range is stored; false otherwise
	 */
	bool empty() const;

	/**
	 * The number of index ranges stored
	 * @return the number of stored index ranges
	 */
	size_t size() const;

	/**
	 * Removes all stored elements
	 */
	void clear();

	/**
	 * Whether or not index ranges are allowed to overlap
	 * @return true if ranges are allowed to overlap; false otherwise
	 */
	bool isAllowingOverlap() const;

	/**
	 * updates the range list data from a valid string encoding (matching regex)
	 * @param stringEncoding the interval list string encoding
	 * @param seq if not NULL, the RnaSequence to be used to shift in/output
	 *         indices to 0-based internal index positions
	 * @throws std::runtime_error if stringEncoding does not match regex
	 */
	void
	fromString( const std::string & stringEncoding, const RnaSequence * seq = NULL );

	/**
	 * shifts all indices by the given value and returns all intervals within
	 * the boundaries [0,indexMax] including indexMax
	 * @param indexShift the shift to be applied to all index range boundaries
	 * @param indexMax the maximal value any index in the returned list can have
	 * @return a new range list with shifted ranges where all ranges shifted to
	 *    indices below 0 are (I) removed if range'.to < 0 or (II) cut to
	 *    range'.from = 0; the same holds respectively for upper bound
	 *    violations exceeding indexMax.
	 */
	IndexRangeList
	shift( const int indexShift, const size_t indexMax ) const;

	/**
	 * Reverses all indices for the given sequence length, i.e.
	 * (newIdx = seqLength-1-oldIdx)
	 */
	IndexRangeList &
	reverseInplace( const size_t seqLength );

	/**
	 * Reverses all indices for the given sequence length, i.e.
	 * (newIdx = seqLength-1-oldIdx)
	 */
	IndexRangeList
	reverse( const size_t seqLength ) const;

	/**
	 * Checks whether or not two range lists are equivalent
	 * @param r the range list to compare to
	 * @return list == r.list
	 */
	const bool operator == ( const IndexRangeList &r ) const {
		return this->list == r.list;
	}

	/**
	 * Checks whether or not two range lists are inequivalent
	 * @param r the range list to compare to
	 * @return list != r.list
	 */
	const bool operator != ( const IndexRangeList &r ) const {
		return this->list != r.list;
	}


	/**
	 * Prints the boundaries of the list's ranges to stream
	 * @param out the ostream to write to
	 * @param l the IndexRangeList object to add
	 * @return the altered stream out
	 */
	friend std::ostream& operator<<(std::ostream& out, const IndexRangeList& l);


protected:

	//! whether or not overlapping ranges are allowed
	bool allowOverlap;

	//! the list of indices
	List list;

};



//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

inline
IndexRangeList::IndexRangeList( const bool allowOverlap_ )
: allowOverlap(allowOverlap_)
, list()
{
}

//////////////////////////////////////////////////////////////////////

inline
IndexRangeList::IndexRangeList( const std::string & stringEncoding
							, const bool allowOverlap_
							, const RnaSequence * seq )
: allowOverlap(allowOverlap_)
, list()
{
	fromString(stringEncoding, seq);
}

//////////////////////////////////////////////////////////////////////

inline
IndexRangeList::IndexRangeList( const IndexRangeList & toCopy )
: allowOverlap(toCopy.allowOverlap)
, list(toCopy.list)
{
}

//////////////////////////////////////////////////////////////////////

inline
IndexRangeList::~IndexRangeList()
{
}

//////////////////////////////////////////////////////////////////////

inline
bool
IndexRangeList::
covers( const size_t from, const size_t to ) const
{
	return covers( IndexRange(from, to) );
}

//////////////////////////////////////////////////////////////////////

inline
IndexRangeList
IndexRangeList::
reverse( const size_t seqLength ) const
{
	// create copy
	IndexRangeList tmp(*this);
	// reverse and return copy
	return tmp.reverseInplace( seqLength );
}

//////////////////////////////////////////////////////////////////////

inline
IndexRangeList::iterator IndexRangeList::erase( IndexRangeList::iterator i ) { return list.erase( i ); }

//////////////////////////////////////////////////////////////////////

inline
IndexRangeList::iterator IndexRangeList::begin() { return list.begin(); }

//////////////////////////////////////////////////////////////////////

inline
IndexRangeList::iterator IndexRangeList::end() { return list.end(); }

//////////////////////////////////////////////////////////////////////

inline
IndexRangeList::const_iterator IndexRangeList::begin() const { return list.begin(); }

//////////////////////////////////////////////////////////////////////

inline
IndexRangeList::const_iterator IndexRangeList::end() const { return list.end(); }

//////////////////////////////////////////////////////////////////////

inline
IndexRangeList::reverse_iterator IndexRangeList::rbegin() { return list.rbegin(); }

//////////////////////////////////////////////////////////////////////

inline
IndexRangeList::reverse_iterator IndexRangeList::rend() { return list.rend(); }

//////////////////////////////////////////////////////////////////////

inline
IndexRangeList::const_reverse_iterator IndexRangeList::rbegin() const { return list.rbegin(); }

//////////////////////////////////////////////////////////////////////

inline
IndexRangeList::const_reverse_iterator IndexRangeList::rend() const { return list.rend(); }

//////////////////////////////////////////////////////////////////////

inline
bool IndexRangeList::empty() const { return list.empty(); }

//////////////////////////////////////////////////////////////////////

inline
bool IndexRangeList::isAllowingOverlap() const { return allowOverlap; }

//////////////////////////////////////////////////////////////////////

inline
size_t IndexRangeList::size() const { return list.size(); }

//////////////////////////////////////////////////////////////////////

inline
void IndexRangeList::clear() { return list.clear(); }

//////////////////////////////////////////////////////////////////////

inline
std::ostream& operator<<(std::ostream& out, const IndexRangeList& l)
{
	// output according to regex
	for (IndexRangeList::const_iterator i=l.begin(); i!=l.end(); i++)
		out <<(i==l.begin()?"":",") <<*i;
	return out;
}

//////////////////////////////////////////////////////////////////////


} // namespace

#endif /* INDEXRANGELIST_H_ */
