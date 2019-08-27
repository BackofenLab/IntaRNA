
#ifndef INTARNA_PREDICTIONTRACKERBASEPAIRPROB_H_
#define INTARNA_PREDICTIONTRACKERBASEPAIRPROB_H_

#include "IntaRNA/PredictionTracker.h"
#include "IntaRNA/InteractionEnergy.h"
#include "IntaRNA/Interaction.h"
#include "IntaRNA/PredictorMfeEns.h"

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
	 * for an interaction by computing the Boltzmann probabilities and
	 * generates a basepair-probabilities dotplot
	 *
	 * @param energy the energy function used for energy calculation
	 * @param fileName the name of the generated postscript file containing
	 * the dotplot
	 */
	PredictionTrackerBasePairProb(
				const InteractionEnergy & energy
				, const std::string & fileName
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
	 * @param seedHandler the seedHandler of the predictor (NULL if noseed predictior)
	 */
	virtual
	void
	updateZ( PredictorMfeEns *predictor, SeedHandler* seedHandler ) override;

	/**
	 * Compute partition function of the
	 * the subregion (i1, j1, i2, j2)
	 * @param i1 region index
	 * @param j1 region index
	 * @param i2 region index
	 * @param j2 region index
	 * @param predictor the predictor providing the probability information
	 * @param seedHandler the seedHandler of the predictor
	 */
	void
	computeMissingZ( const size_t i1, const size_t j1
								, const size_t i2, const size_t j2
								, PredictorMfeEns *predictor
								, SeedHandler* seedHandler );

	/**
	 * Check if full seed exists in the region (i1, j1, i2, j2)
	 * @param i1 region index
	 * @param j1 region index
	 * @param i2 region index
	 * @param j2 region index
	 * @param seedHandler the seedHandler of the predictor
	 *
	 * @return true if region contains a full seed
	 */
	bool
	isFullSeedinRegion( const size_t i1, const size_t j1
					        	, const size_t i2, const size_t j2
					        	, SeedHandler* seedHandler );

	/**
	 * Get leftmost seeds containing basepair (k1, k2)
	 * @param k1 base pair index
	 * @param k2 base pair index
	 * @param seedHandler the seedHandler of the predictor
	 *
	 * @return vector of pairs containing the left-most seed of each
	 *         loop-overlapping cluster of seeds containing base pair k
	 */
	std::vector< std::pair <size_t, size_t> >
	getLeftMostSeedsAtK( const size_t k1, const size_t k2
					           , SeedHandler* seedHandler );

	/**
	 * Get partial energy of seed (si1, si2) at region (i1, j1, i2, j2)
	 * @param i1 region index
	 * @param j1 region index
	 * @param i2 region index
	 * @param j2 region index
	 * @param seedHandler the seedHandler of the predictor
	 *
	 * @return energy of seed at given region
	 */
	E_type
	getPartialSeedEnergy( const size_t si1, const size_t si2
		                  , const size_t i1, const size_t j1
											, const size_t i2, const size_t j2
											, SeedHandler* seedHandler );

	/**
	 * Access to the current partition function covering
	 * the interaction at region (i1, j1, i2, j2).
	 * @param i1 region index
	 * @param j1 region index
	 * @param i2 region index
	 * @param j2 region index
	 * @param predictor the predictor providing the probability information
	 *
	 * @return the hybridization partition function at given region
	 */
	Z_type
	getHybridZ( const size_t i1, const size_t j1
	          , const size_t i2, const size_t j2
	          , PredictorMfeEns *predictor);

	/**
	 * Set the current partition function covering
	 * the interaction at region (i1, j1, i2, j2).
	 * @param i1 region index
	 * @param j1 region index
	 * @param i2 region index
	 * @param j2 region index
	 * @param partZ the new partition function for givent region
	 */
	void
	updateHybridZ( const size_t i1, const size_t j1
						 	 , const size_t i2, const size_t j2
							 , const Z_type partZ );

protected:

	struct StructureProb {
		size_t i1;
		size_t j1;
		size_t i2;
		size_t j2;
		float prob;
	};

	//! data container to encode a site with respective partition function
	struct ZPartition {
		size_t i1;
		size_t j1;
		size_t i2;
		size_t j2;
		Z_type partZ;
	};

	//! energy handler used for predictions
	const InteractionEnergy & energy;

	//! filename of the generated dotplot
	const std::string fileName;

	//! threshold used to draw probabilities in dotplot
	const Z_type probabilityThreshold;

	//! map storing structure probabilities
	std::unordered_map<size_t, StructureProb> structureProbs;

	//! map storing missing Z partitions for a given interaction
	std::unordered_map<size_t, ZPartition> Z_partitions;

	/**
	 * Generates key for storing values in map
	 */
	virtual
	size_t
	generateMapKey( const size_t i1, const size_t j1
						, const size_t i2, const size_t j2 ) const;

};

//////////////////////////////////////////////////////////////////////

} // namespace

#endif /* INTARNA_PREDICTIONTRACKERBASEPAIRPROB_H_ */
