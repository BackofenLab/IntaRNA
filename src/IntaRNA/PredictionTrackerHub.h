
#ifndef INTARNA_PREDICTIONTRACKERHUB_H_
#define INTARNA_PREDICTIONTRACKERHUB_H_

#include "IntaRNA/PredictionTracker.h"

#include <list>

namespace IntaRNA {

/**
 * Hub implementation of an PredictionTracker that forwards each reported interaction
 * to each of the registered PredictionTracker instances.
 *
 */
class PredictionTrackerHub: public PredictionTracker
{
protected:

	//! list of PredictionTracker objects to forward to
	std::list< PredictionTracker * > trackList;

	//! whether or not the elements in trackList are to be deleted when this
	//! object is destroyed
	bool deleteTrackListOnDestruction;

public:

	/**
	 * construction
	 * @param deleteTrackListOnDestruction whether or not the objects the hub is
	 *        forwarding to are to be deleted on destruction
	 */
	PredictionTrackerHub( const bool deleteTrackListOnDestruction = true );

	/**
	 * copy construction (will not delete the copied PredictionTracker list!)
	 *
	 * @param toCopy the hub to copy
	 */
	PredictionTrackerHub( const PredictionTrackerHub & toCopy );


	/**
	 * destruction
	 */
	virtual ~PredictionTrackerHub();

	/**
	 * Informs all registered trackers that Predictor.updateOptima() was called
	 * with the given data.
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1
	 * @param energy the overall energy of the interaction site
	 */
	virtual
	void
	updateOptimumCalled( const size_t i1, const size_t j1
						, const size_t i2, const size_t j2
						, const E_type energy
						);

	/**
	 * Adds a new PredictionTracker to the forwarding list.
	 * @param tracker pointer to the tracker to forward to
	 */
	virtual
	void
	addPredictionTracker( PredictionTracker * tracker );

	/**
	 * Defines whether or not the PredictionTracker instances this hub forwards to
	 * are deleted when this hub is deleted.
	 * @param deleteOnDestruction true: forward objects are deleted; false: not
	 */
	virtual
	void
	setDeleteOnDestruction( const bool deleteOnDestruction );

	/**
	 * Returns whether or not the PredictionTracker instances this hub forwards to
	 * are deleted when this hub is deleted.
	 * @return true: forward objects are deleted; false: not
	 */
	virtual
	bool
	isDeleteOnDestruction() const;


	/**
	 * Copies the trackList of toCopy but wont delete it later.
	 *
	 * NOTE: this (copy) hub will NOT delete the trackList to avoid duplicated
	 * deletions. This should be done by the original hub (toCopy)
	 *
	 * @param toCopy the tracker to make this a copy of
	 * @return the altered range object (*this)
	 */
	PredictionTrackerHub &
	operator= ( const PredictionTrackerHub & toCopy);


	/**
	 * Whether or not a tracker is registered for forwarding.
	 * @return true if any tracker is registered; false otherwise
	 */
	bool
	empty() const;

};


/////////////////////////////////////////////////////////////////////////

inline
PredictionTrackerHub
::PredictionTrackerHub( const bool deleteOnDestruction )
 : PredictionTracker()
	, trackList()
	, deleteTrackListOnDestruction(deleteOnDestruction)
{
}

/////////////////////////////////////////////////////////////////////////

inline
PredictionTrackerHub
::PredictionTrackerHub( const PredictionTrackerHub & toCopy )
 : PredictionTracker()
	, trackList( toCopy.trackList )
	// dont delete on destruction (should be done by original) to avoid double cleanup
	, deleteTrackListOnDestruction(false)
{
}

/////////////////////////////////////////////////////////////////////////

inline
PredictionTrackerHub::
~PredictionTrackerHub()
{
	// check if something is to be done
	if (isDeleteOnDestruction()) {
		// delete all tracker
		for( auto it = trackList.begin(); it != trackList.end(); it++ ) {
			 INTARNA_CLEANUP( (*it) );
		}
	}
}



/////////////////////////////////////////////////////////////////////////

inline
void
PredictionTrackerHub::
addPredictionTracker( PredictionTracker * tracker )
{
#if INTARNA_IN_DEBUG_MODE
	if (tracker != NULL) {
		// check if already present
		for (auto it=trackList.begin(); it!=trackList.end(); it++) {
			if (tracker == (*it)) {
				throw std::runtime_error("PredictionTrackerHub::addPredictionTracker() : repeated adding of the same PredictionTracker pointer");
			}
		}
	}
#endif
	// add if not null
	if (tracker != NULL) {
		trackList.push_back(tracker);
	}
}

/////////////////////////////////////////////////////////////////////////

inline
void
PredictionTrackerHub::
updateOptimumCalled( const size_t i1, const size_t j1
					, const size_t i2, const size_t j2
					, const E_type energy
					)
{
	// forward to all in list
	for (auto trackIt=trackList.begin(); trackIt!=trackList.end(); trackIt++) {
		(*trackIt)->updateOptimumCalled(i1,j1,i2,j2,energy);
	}
}

/////////////////////////////////////////////////////////////////////////

inline
void
PredictionTrackerHub::
setDeleteOnDestruction( const bool deleteOnDestruction )
{
	this->deleteTrackListOnDestruction = deleteOnDestruction;
}

/////////////////////////////////////////////////////////////////////////

inline
bool
PredictionTrackerHub::
isDeleteOnDestruction() const
{
	return this->deleteTrackListOnDestruction;
}

/////////////////////////////////////////////////////////////////////////

inline
PredictionTrackerHub &
PredictionTrackerHub::
operator= ( const PredictionTrackerHub & toCopy)
{
	// clean up current content
	if (this->deleteTrackListOnDestruction) {
		for (auto it=trackList.begin(); it!=trackList.end(); it++) {
			 INTARNA_CLEANUP( (*it) );
		}
	}
	trackList.clear();

	// copy data
	trackList = toCopy.trackList;
	// ensure this copy doesnt delete the list content
	deleteTrackListOnDestruction = false;

	// modified object access
	return *this;
}

/////////////////////////////////////////////////////////////////////////

inline
bool
PredictionTrackerHub::
empty() const
{
	// return whether or not the tracker list is empty
	return trackList.empty();
}

/////////////////////////////////////////////////////////////////////////

} // namespace

#endif /* PREDICTIONTRACKERHUB_H_ */
