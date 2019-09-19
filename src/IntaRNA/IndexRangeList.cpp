
#include "IntaRNA/IndexRangeList.h"

#include <algorithm>
#include <exception>

namespace IntaRNA {

//////////////////////////////////////////////////////////////////////

const std::string IndexRangeList::regexString("("+std::string(IndexRange::regexString)+",)*"+std::string(IndexRange::regexString));

const boost::regex IndexRangeList::regex(IndexRangeList::regexString);

//////////////////////////////////////////////////////////////////////

bool
IndexRangeList::
covers( const size_t index ) const
{
	// quick check
	if (list.empty()) {
		return false;
	}

	// find first range that with begin > index
	const_iterator r = std::upper_bound( list.begin(), list.end(), IndexRange(index,std::numeric_limits<size_t>::max()) );
	if ( r == list.begin() ) {
		return false;
	} else {
		// go to preceding range and check if <= the end of the blocked range
		return index <= (--r)->to;
	}
}

//////////////////////////////////////////////////////////////////////

bool
IndexRangeList::
covers( const IndexRange & range ) const
{
	// quick check
	if (list.empty()) {
		return false;
	}
	// find first range with (begin > range.from) or (begin == range.from && end > range.from)
	const_iterator r = std::upper_bound( list.begin(), list.end(), range );

	// check if covered by trailing range
	if ( r != list.end() && r->from <= range.from && range.to <= r->to) {
		return true;
	}

	if (allowOverlap) {

		// check ALL preceding ranges
		while( r-- != list.begin() ) {
			if ( r->from <= range.from && r->to >= range.to) {
				return true;
			}
		}

	} else {

		// simply check only the preceding range if any
		if ( r != list.begin() ) {
			// go to preceding range and check
			--r;
			if (r->from <= range.from && range.to <= r->to) {
				return true;
			}
		}
	}
	// default assumption
	return false;
}

//////////////////////////////////////////////////////////////////////

bool
IndexRangeList::
overlaps( const IndexRange& range ) const
{
#if INTARNA_IN_DEBUG_MODE
	if (!range.isAscending())  {
		throw std::runtime_error("IndexRangeList::overlaps("+toString(range)+") range is not ascending");
	}
#endif

	// quick check
	if (list.empty()) {
		return false;
	}

	// find first range that with begin > range.from
	const_iterator r = std::upper_bound( list.begin(), list.end(), range );

	// check if overlapping with trailing range
	if ( r != list.end() && r->from <= range.to) {
		return true;
	}

	if (allowOverlap) {

		// check ALL ranges
		for ( auto & r : list ) {
			if ( (r.from >= range.from && r.from <= range.to)
				|| (r.to >= range.from && r.to <= range.to)
				|| (r.from <= range.from && r.to >= range.to)	)
			{
				return true;
			}
			// check if subsequent ranges can not overlap anymore
			if (r.from > range.to) {
				break;
			}
		}

	} else { // non-overlapping ranges

		// check if not overlapping with predecessor
		if ( r != list.begin()) {
			// go to predecessor
			--r;
			// check for overlap
			if ( r->to >= range.from ) {
				return true;
			}
		}
	}

	// default assumption
	return false;

}

//////////////////////////////////////////////////////////////////////

void
IndexRangeList::
push_back( const IndexRange& range )
{
	if (!range.isAscending())  {
		throw std::runtime_error("IndexRangeList::push_back("+toString(range)+") range is not ascending");
	}
	if (!list.empty() && list.rbegin()->from >= range.from) {
		throw std::runtime_error("IndexRangeList::push_back("+toString(range)+") violates order given last range = "+toString(*(list.rbegin())));
	}
	if (!allowOverlap && !list.empty() && list.rbegin()->to >= range.from) {
		throw std::runtime_error("IndexRangeList::push_back() : not allowed overlap with existing range");
	}
	// sorting should be OK (in debug mode.. ;) )
	list.push_back( range );
}

//////////////////////////////////////////////////////////////////////

IndexRangeList::iterator
IndexRangeList::
insert( const IndexRange& range )
{
	if (!range.isAscending())  {
		throw std::runtime_error("IndexRangeList::insert("+toString(range)+") range is not ascending");
	}
	// init list
	if (list.empty()) {
		list.push_back( range );
		// iterator to first element
		return begin();
	} else
	// init list or add at the end to list
	if ( *(list.rbegin()) < range ) {
		list.push_back( range );
		// iterator to last element
		iterator idxIter = list.begin();
		std::advance( idxIter, list.size()-1 );
		return idxIter;
	} else
	// insert accordingly
	{
		// find first range that with begin < from or begin==from && end>to
		List::iterator r = std::upper_bound( list.begin(), list.end(), range );
		if (r != list.end()) {
			// check for overlap
			if (range.to >= r->from && !allowOverlap) {
				throw std::runtime_error("IndexRangeList::insert("+toString(range)+") range is overlapping with "+toString(*r)+", which is not allowed");
			}
			return list.insert( r, range );
		}
		// check if already existing (predecessor)
		if (r != list.begin()){
			--r;
			// check if already present
			if (*r == range) {
				// return iterator to already present element
				return r;
			}
			// check for overlap
			if (r->to >= range.from && !allowOverlap) {
				throw std::runtime_error("IndexRangeList::insert("+toString(range)+") range is overlapping with "+toString(*r)+", which is not allowed");
			}
			++r;
		}
		// insert accordingly preserving sorting
		return list.insert( r, range );
	}
}

//////////////////////////////////////////////////////////////////////

IndexRange &
IndexRangeList::
get( const size_t idx )
{
	if (idx >= size()) {
		throw std::runtime_error("IndexRangeList::get() : index "+toString(idx)+" out of range (>= "+toString(size())+")");
	}
	// access according element via iterator and return
	iterator idxIter = list.begin();
	std::advance( idxIter, idx );
	return *(idxIter);
}

//////////////////////////////////////////////////////////////////////

const IndexRange &
IndexRangeList::
get( const size_t idx ) const
{
	if (idx >= size()) {
		throw std::runtime_error("IndexRangeList::get() : index "+toString(idx)+" out of range (>= "+toString(size())+")");
	}
	// access according element via iterator and return
	const_iterator idxIter = list.begin();
	std::advance( idxIter, idx );
	return *(idxIter);
}

//////////////////////////////////////////////////////////////////////

IndexRangeList
IndexRangeList::
shift( const int indexShift, const size_t indexMax ) const
{
	IndexRangeList l;
	IndexRange r2;
	for (IndexRangeList::const_iterator r=begin(); r!=end(); r++) {
		// skip ranges leaving the valid interval
		if ( (indexShift+((int)r->to)) < 0 || (int)indexMax < (indexShift+((int)r->from))) {
			continue;
		}
		// get shifted range boundaries
		r2.from = (size_t)std::max(0,indexShift+((int)r->from));
		r2.to = (size_t)std::min((int)indexMax,indexShift+((int)r->to));
		// store shifted and cut range
		l.insert(r2);
	}

	// final updated range list
	return l;
}

//////////////////////////////////////////////////////////////////////

IndexRangeList &
IndexRangeList::
reverseInplace( const size_t seqLength )
{
	// reverse each entry
	size_t tmpFrom;
	for (IndexRangeList::iterator r=begin(); r!=end(); r++) {
		// check if reversal is possible
		if (r->from >= seqLength || r->to >= seqLength) throw std::runtime_error("IndexRangeList::reverse("+toString(seqLength)+") = range "+toString(*r)+" exceeds seqLength");
		// reverse boundaries
		tmpFrom = r->from;
		r->from = seqLength -1 - r->to;
		r->to = seqLength -1 - tmpFrom;
	}
	// reverse order of list entries
	list.reverse();
	// return access to altered element
	return *this;
}

//////////////////////////////////////////////////////////////////////

void
IndexRangeList::
fromString( const std::string & stringEncoding, const RnaSequence * seq )
{
	// clear current data
	this->clear();
	// check if something to be parsed
	if (!stringEncoding.empty()) {
		// check if parsable
		if( ! boost::regex_match(stringEncoding, IndexRangeList::regex, boost::match_perl) ) {
			throw std::runtime_error("IndexRangeList::fromString("+stringEncoding+") uses no valid index range string encoding matching '"+regex.str()+"'");
		}
		// find split position
		size_t startPos = 0, splitPos = std::string::npos;
		while (startPos != splitPos) {
			splitPos = stringEncoding.find(',',startPos);
			// insert interval
			this->insert( IndexRange(stringEncoding.substr(startPos,splitPos-(splitPos==std::string::npos?0:startPos)),seq));
			// update start of next interval encoding to parse
			startPos = splitPos + (splitPos != std::string::npos ? 1 : 0);
		}
	}
}

//////////////////////////////////////////////////////////////////////



} // namespace
