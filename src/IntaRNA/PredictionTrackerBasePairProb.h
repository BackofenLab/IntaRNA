
#ifndef INTARNA_PREDICTIONTRACKERBASEPAIRPROB_H_
#define INTARNA_PREDICTIONTRACKERBASEPAIRPROB_H_

#include "IntaRNA/PredictionTracker.h"
#include "IntaRNA/InteractionEnergy.h"
#include "IntaRNA/Interaction.h"

#include <iostream>

#include <boost/algorithm/string.hpp>

namespace IntaRNA {

/**
 * Collects partition function parts Z(i,j), computes base-pair probabilities
 * and prints them as a dot plot
 *
 * The information is written to stream on destruction.
 */
class PredictionTrackerBasePairProb: public PredictionTracker
{

public:

	/**
	 * Constructs a PredictionTracker that collects probability information
	 * for a set of interaction spots by computing the Boltzmann probability
	 * of any interaction enclosing the positions.
	 *
	 * @param energy the energy function used for energy calculation
	 * @param spots the interaction spots to track for their probabilities
	 * @param outStreamName the stream name where the probability data
	 *        is to be written to. use STDOUT/STDERR for the respective stream.
	 *        Otherwise, an according file is created
	 */
	PredictionTrackerBasePairProb(
				const InteractionEnergy & energy
				, const std::string & outStreamName
			);

	/**
	 * Constructs a PredictionTracker that collects probability information
	 * for a set of interaction spots by computing the Boltzmann probability
	 * of any interaction enclosing the positions.
	 *
	 * @param energy the energy function used for energy calculation
	 * @param spots the interaction spots to track for their probabilities
	 * @param outStream the stream where the probability data
	 *        is to be written to.
	 */
	PredictionTrackerBasePairProb(
				const InteractionEnergy & energy
				, std::ostream & outStream
			);


	/**
	 * destruction: write the probabilities to stream.
	 */
	virtual ~PredictionTrackerBasePairProb();


	/**
	 * Updates the probability information for each Predictor.updateOptima() call.
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
	 * Updates the probability information.
	 *
	 * @param predictor the predictor providing the probability information
	 */
	virtual
	void
	updateZ( const PredictorMfeEns *predictor ) override;

private:

	/**
	 * Generates key for storing values in map
	 */
	virtual
	size_t
	generateMapKey( const size_t i1, const size_t j1
						, const size_t i2, const size_t j2 );

	/**
	 * Recursively computes structure probabilities
	 */
	virtual
	void
	computeStructureProb( const size_t i1, const size_t j1
						, const size_t i2, const size_t j2
						, const size_t i1last, const size_t j1last
						, const size_t i2last, const size_t j2last, const size_t key, const size_t lastKey );

protected:

	struct ZPartition {
		size_t i1;
		size_t j1;
		size_t i2;
		size_t j2;
		Z_type partZ;
	};

	struct StructureProb {
		size_t i1;
		size_t j1;
		size_t i2;
		size_t j2;
		float prob;
	};

	//! energy handler used for predictions
	const InteractionEnergy & energy;

	//! the stream to write the probabilities to
	std::ostream * outStream;

	//! whether or not outStream is to be deleted on destruction
	const bool deleteOutStream;

	//! overall partition function of all reported interactions
	Z_type overallZ;

	//! map storing hybridZ partitions for a given interaction
	std::unordered_map<size_t, ZPartition> Z_partitions;

	//! map storing structure probabilities
	std::unordered_map<size_t, StructureProb> structureProbs;

};

//////////////////////////////////////////////////////////////////////

} // namespace

#endif /* INTARNA_PREDICTIONTRACKERBASEPAIRPROB_H_ */
