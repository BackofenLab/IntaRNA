
#ifndef OUTPUTHANDLERINTERACTIONLIST_H_
#define OUTPUTHANDLERINTERACTIONLIST_H_

#include "OutputHandler.h"

#include <list>

namespace IntaRNA
{

/**
 * OutputHandler that stores all reported interactions in a list
 *
 * NOTE: this handler does NOT support InteractionRange instances
 *
 */
class OutputHandlerInteractionList: public OutputHandler
{
protected:

	//! the container used internally for storing interactions
	typedef std::list<Interaction*> StorageContainer;

public:

	//! const iterator on stored interactions
	typedef StorageContainer::const_iterator const_iterator;
	//! iterator on stored interactions
	typedef StorageContainer::iterator iterator;

public:

	/**
	 * construction of empty container
	 * @param maxToStore the maximal number of (best) elements to be stored
	 * @param outConstraint the output constraint applied to find the reported
	 *        interaction
	 */
	OutputHandlerInteractionList( const OutputConstraint & outConstraint
								, const size_t maxToStore);

	//! destruction
	virtual ~OutputHandlerInteractionList();

	/**
	 * Adds a given RNA-RNA interaction to the storage/output.
	 *
	 * @param interaction the interaction to add
	 */
	virtual
	void
	add( const Interaction & interaction );

	//! whether or not the container is empty
	//! @return true if no elements are stored
	bool empty() const;

	//! constant iterator to first stored interaction or end() if empty
	//! @return first stored interaction (pointer) or end() if empty
	const_iterator begin() const;

	//! constant iterator pointing AFTER last interaction (non-valid)
	//! @return iterator pointing after last element
	const_iterator end() const;

	//! iterator to first stored interaction or end() if empty
	//! @return first stored interaction (pointer) or end() if empty
	iterator begin();

	//! iterator pointing AFTER last interaction (non-valid)
	//! @return iterator pointing after last element
	iterator end();

protected:

	//! counter of reported interactions
	using OutputHandler::reportedInteractions;

	//! container where interactions are stored
	StorageContainer storage;

	//! number of elements to store within storage (only maxToStore best)
	const size_t maxToStore;

	//! 'less-than' comparison struct dereferencing storage pointer
	static
	bool lessThan_StorageContainer( const Interaction * const a
			, const Interaction * const b )
	    { return *a < *b; }
};

////////////////////////////////////////////////////////////////////////////

inline
bool
OutputHandlerInteractionList::
empty() const {
	return storage.empty();
}

////////////////////////////////////////////////////////////////////////////

inline
OutputHandlerInteractionList::
const_iterator
OutputHandlerInteractionList::
begin() const {
	return storage.begin();
}

////////////////////////////////////////////////////////////////////////////

inline
OutputHandlerInteractionList::
const_iterator
OutputHandlerInteractionList::
end() const {
	return storage.end();
}

////////////////////////////////////////////////////////////////////////////

inline
OutputHandlerInteractionList::
iterator
OutputHandlerInteractionList::
begin() {
	return storage.begin();
}

////////////////////////////////////////////////////////////////////////////

inline
OutputHandlerInteractionList::
iterator
OutputHandlerInteractionList::
end() {
	return storage.end();
}

////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////

} /* namespace IntaRNA */

#endif /* OUTPUTHANDLERINTERACTIONLIST_H_ */
