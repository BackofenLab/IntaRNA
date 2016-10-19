
#ifndef INDEXRANGELIST_H_
#define INDEXRANGELIST_H_

#include "IndexRange.h"

#include <vector>

class IndexRangeList
{
protected:

	//! List of ranges
	typedef std::vector< IndexRange > List;

public:

	typedef List::iterator iterator;
	typedef List::const_iterator const_iterator;
	typedef List::reverse_iterator reverse_iterator;
	typedef List::const_reverse_iterator const_reverse_iterator;

public:

	/**
	 * empty construction
	 */
	IndexRangeList();

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
	 * adds an ascending range to the end of the list.
	 * NOTE: if the insertion would violate range sorting, an exception is raised
	 * @param range the ascending index range to add (from should be > rbegin()->to)
	 */
	void push_back( const IndexRange& range );

	/**
	 * adds an ascending range to the according position of the sorted list.
	 * @param range the ascending index range to add
	 */
	iterator insert( const IndexRange& range );

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

protected:

	//! the list of indices
	List list;

};

#endif /* INDEXRANGELIST_H_ */
