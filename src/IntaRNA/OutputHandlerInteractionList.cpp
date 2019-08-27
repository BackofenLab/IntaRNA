
#include "OutputHandlerInteractionList.h"

#include <algorithm>

namespace IntaRNA
{

/////////////////////////////////////////////////////////////////////////////

OutputHandlerInteractionList::
OutputHandlerInteractionList(const OutputConstraint & outConstraint
							, const size_t maxToStore)
 :	OutputHandler(outConstraint)
	, storage()
	, maxToStore(maxToStore)
{
}

/////////////////////////////////////////////////////////////////////////////

OutputHandlerInteractionList::
~OutputHandlerInteractionList()
{
	// cleanup stored interactions
	for( auto it = storage.begin(); it != storage.end(); it++ ) {
		delete (*it);
	}
	storage.clear();
}

/////////////////////////////////////////////////////////////////////////////

void
OutputHandlerInteractionList::
add( const Interaction & interaction )
{
	if (interaction.isEmpty()) {
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_OutputHandlerInteractionListUpdate)
#endif
		{
		// count interaction
		reportedInteractions++;
		}
		return;
	}

#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_OutputHandlerInteractionListUpdate)
#endif
	{
		// count interaction
		reportedInteractions++;
		if (storage.size() < maxToStore || lessThan_StorageContainer( &interaction, *(storage.rbegin()) )) {
			// find where to insert this interaction
			StorageContainer::iterator insertPos = std::lower_bound( storage.begin(), storage.end(), &interaction, lessThan_StorageContainer );
			// check if interaction is NOT already part of the list
			if ( insertPos == storage.end() || lessThan_StorageContainer( &interaction, *(insertPos) ) ) {
				// insert current interaction
				storage.insert( insertPos, new Interaction(interaction) );
				// remove last element if needed
				if (storage.size() > maxToStore) {
					// delete object
					delete (*(storage.rbegin()));
					// remove pointer
					storage.resize(storage.size()-1);
				}
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////////

} /* namespace IntaRNA */
