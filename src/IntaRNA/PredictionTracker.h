
#ifndef INTARNA_PREDICTIONTRACKER_H_
#define INTARNA_PREDICTIONTRACKER_H_


#include "IntaRNA/general.h"

namespace IntaRNA {

/**
 * Generic interface to track prediction progress of Predictor instances.
 *
 */
class PredictionTracker
{

public:

	/**
	 * construction
	 */
	PredictionTracker();

	/**
	 * destruction
	 */
	virtual ~PredictionTracker();

	/**
	 * Informs the tracker that Predictor.updateOptima() was called with the
	 * given data.
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
						) = 0;

};



///////////////////////////////////////////////////////////////////////////

inline
PredictionTracker::PredictionTracker()
{
}

///////////////////////////////////////////////////////////////////////////

inline
PredictionTracker::~PredictionTracker()
{
}

///////////////////////////////////////////////////////////////////////////


} // namespace

#endif /* PREDICTIONTRACKER_H_ */
