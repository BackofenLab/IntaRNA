
#include "OutputHandlerInteractionList.h"

#include <algorithm>

namespace IntaRNA
{

/////////////////////////////////////////////////////////////////////////////

OutputHandlerInteractionList::
OutputHandlerInteractionList(const size_t maxToStore)
 :	storage()
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
			// remove last element if needed
			if (storage.size() >= maxToStore) {
				// delete object
				delete (*(storage.rbegin()));
				// remove pointer
				storage.resize(storage.size()-1);
			}
			// find where to insert this interaction
			StorageContainer::iterator insertPos = std::lower_bound( storage.begin(), storage.end(), &interaction, lessThan_StorageContainer );
			// insert current interaction
			storage.insert( insertPos, new Interaction(interaction) );
		}
	}
}

/////////////////////////////////////////////////////////////////////////////

void
OutputHandlerInteractionList::
add( const InteractionRange & range )
{
	INTARNA_NOT_IMPLEMENTED("OutputHandlerInteractionList::add( const InteractionRange & range )");
}

/////////////////////////////////////////////////////////////////////////////

} /* namespace IntaRNA */
