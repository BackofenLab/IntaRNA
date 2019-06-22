
#ifndef INTARNA_PREDICTORMFEENS_H_
#define INTARNA_PREDICTORMFEENS_H_


#include "IntaRNA/PredictorMfe.h"

#include "IntaRNA/IndexRangeList.h"

#include <list>
#include <utility>

namespace IntaRNA {

/**
 * Generic Predictor interface for ensemble-based MFE interaction computation
 *
 * @author Martin Raden
 *
 */
class PredictorMfeEns : public PredictorMfe {


public:

	/**
	 * Construction call the super constructor
	 *
	 * @param energy the interaction energy handler
	 * @param output the output handler to report optimal interactions to
	 * @param predTracker the prediction tracker to be used or NULL if no
	 *         tracking is to be done; if non-NULL, the tracker gets deleted
	 *         on this->destruction.
	 */
	PredictorMfeEns( const InteractionEnergy & energy
				, OutputHandler & output
				, PredictionTracker * predTracker );

	virtual ~PredictorMfeEns();

	/**
	 * Access to the current overall partition function covering
	 * all interactions of the last predict() call.
	 *
	 * @return the overall hybridization partition function
	 */
	Z_type
	getOverallZ() const;

	/**
	 * Report Z information to the prediction trackers
	 */
	void
	reportZ() const;

protected:

// TODO move to subclass ?!
	//! data container to encode a site with respective partition function
	struct ZPartition {
		size_t i1;
		size_t j1;
		size_t i2;
		size_t j2;
		Z_type partZ;
	};

	//! access to the interaction energy handler of the super class
	using PredictorMfe::energy;

	//! access to the output handler of the super class
	using PredictorMfe::output;

	//! access to the prediction tracker of the super class
	using PredictorMfe::predTracker;

	//! the overall partition function since initZ() was last called.
	//! its value is updated by updateZ()
	Z_type overallZ;

// TODO move to subclass...
	//! map storing Z partitions for a given interaction
	std::unordered_map<size_t, ZPartition> Z_partitions;


	/**
	 * Initializes the overall hybridization partition function as well as
	 * the global energy minimum storage. Will be called by predict().
	 *
	 * @param outConstraint constrains the interactions reported to the output handler
	 */
	virtual
	void
	initZ( const OutputConstraint & outConstraint );

	/**
	 * Check if energy maxLength exceeds allowed limit for key generation
	 * Throw runtime error if exceeding limit
	 * @param maxLength the maximal length of considered subsequences
	 */
	static
	void
	checkKeyBoundaries( const size_t maxLength );

	/**
	 * Updates the overall hybridization partition function.
	 *
	 * Note: if called multiple time for the same boundaries then the
	 * reported partition functions have to represent disjoint interaction sets!
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1
	 * @param partFunct the partition function of the interaction
	 * @param isHybridZ whether or not the given Z is only covering
	 *        hybridization energy terms (init+loops) or the total
	 *        interaction energy
	 */
	virtual
	void
	updateZ( const size_t i1, const size_t j1
				, const size_t i2, const size_t j2
				, const Z_type partFunct
				, const bool isHybridZ );


	/**
	 * Reports the overall partition function and calls
	 * PredictorMfe::reportOptima().
	 *
	 * @param outConstraint constrains the interactions reported to the output handler
	 */
	virtual
	void
	reportOptima( const OutputConstraint & outConstraint );

};

} // namespace

#endif /* INTARNA_PREDICTORMFEENS_H_ */
