
#ifndef INTARNA_PREDICTORMFEENS2DSEEDEXTENSION_H_
#define INTARNA_PREDICTORMFEENS2DSEEDEXTENSION_H_

#include "IntaRNA/PredictorMfeEns.h"
#include "IntaRNA/BasePairProbabilities.h"
#include "IntaRNA/Matrix.h"
#include "IntaRNA/SeedHandlerIdxOffset.h"
#include <memory>

namespace IntaRNA {

/**
 * Implements seed-based space-efficient interaction prediction
 * based on minimizing ensemble free energy of interaction sites.
 *
 * Note, for each seed start (i1,i2) only the mfe seed is considered for the
 * overall interaction computation instead of considering all possible seeds
 * starting at (i1,i2).
 *
 * @author Frank Gelhausen
 * @author Martin Raden
 *
 */
class PredictorMfeEns2dSeedExtension: public PredictorMfeEns {

protected:

	//! matrix type to hold the partition functions for interaction site starts
	typedef Matrix<Z_type> Z2dMatrix;

public:


	/**
	 * Constructs a predictor and stores the energy and output handler
	 *
	 * @param energy the interaction energy handler
	 * @param output the output handler to report mfe interactions to
	 * @param predTracker the prediction tracker to be used or NULL if no
	 *         tracking is to be done; if non-NULL, the tracker gets deleted
	 *         on this->destruction.
	 * @param pairProbabilities optional non-owning sink; must outlive predictions
	 *        and remain pending until all disjoint regions have succeeded
	 * @param seedHandler the seed handler to be used
	 */
	PredictorMfeEns2dSeedExtension(
			const InteractionEnergy & energy
			, OutputHandler & output
			, PredictionTracker * predTracker
			, SeedHandler * seedHandler
			, BasePairProbabilities * pairProbabilities = nullptr );


	/**
	 * data cleanup
	 */
	virtual ~PredictorMfeEns2dSeedExtension();


	/**
	 * Computes the mfe for the given sequence ranges (i1-j1) in the first
	 * sequence and (i2-j2) in the second sequence and reports it to the output
	 * handler.
	 *
	 * Each considered interaction contains a seed according to the seed handler
	 * constraints.
	 *
	 * @param r1 the index range of the first sequence interacting with r2
	 * @param r2 the index range of the second sequence interacting with r1
	 *
	 */
	virtual
	void
	predict( const IndexRange & r1 = IndexRange(0,RnaSequence::lastPos)
			, const IndexRange & r2 = IndexRange(0,RnaSequence::lastPos) );


protected:

	/** Shared forward/outside objective coefficient. Override here for supported
	 * reweighting; return zero for forbidden boundaries. Coordinates are local.
	 * @return complete boundary Boltzmann factor (ED, dangles, ends, energyAdd)
	 */
	virtual Z_type exactBoundaryWeight(size_t i1,size_t j1,size_t i2,size_t j2) const;
	/** Run the selected disjoint stack-seed partition backend. */
	void predictStackSeeds(const IndexRange & r1,const IndexRange & r2);

	//! Optional signed reverse trace for the legacy anchored extension recurrences.
	struct ExtensionProbabilityTrace;
	std::unique_ptr<ExtensionProbabilityTrace> extensionProbabilityTrace;
	/** Begin a regional trace; no allocation if no probability sink was requested. */
	void initExtensionProbabilities(size_t n,size_t m);
	/** Store the current anchor seed's actual pairs, including its endpoints. */
	void beginExtensionProbabilitySeed(size_t i1,size_t i2);
	/** Add one accepted left/seed/right objective to the optional reverse trace.
	 * All coordinates are local, with seed boundaries in increasing DP order.
	 * @param left whether the left matrix is a factor (otherwise use initiation)
	 * @param right whether the right matrix is a factor
	 */
	void addExtensionProbabilityRoot(size_t i1,size_t j1,size_t i2,size_t j2,
			Z_type partZ,bool left,bool right);
	/** Reverse both extension traces and accumulate actual pair numerators. */
	void finishExtensionProbabilitySeed();
	/** Validate and atomically commit this region in original sequence order. */
	void commitExtensionProbabilities();
	/** Actual seed pairs strictly before the given local target coordinate. */
	std::vector<Interaction::BasePair> extensionSeedPrefix(size_t i1,size_t i2,size_t before) const;

	//! optional sequence-pair owner, never deleted or finalized by the predictor
	BasePairProbabilities * pairProbabilities;

	//! access to the interaction energy handler of the super class
	using PredictorMfeEns::energy;

	//! access to the output handler of the super class
	using PredictorMfeEns::output;

	//! partition function of all interaction hybrids that start on the left side of the seed including E_init
	Z2dMatrix hybridZ_left;

	//! the seed handler (with idx offset)
	SeedHandlerIdxOffset seedHandler;

	//! partition function of all interaction hybrids that start on the right side of the seed excluding E_init
	Z2dMatrix hybridZ_right;

protected:

	/**
	 * Computes all entries of the hybridE matrix for interactions ending in
	 * p=j1 and q=j2 and report all valid interactions to updateOptima()
	 *
	 * @param j1 end of the interaction within seq 1
	 * @param j2 end of the interaction within seq 2
	 *
	 */
	virtual
	void
	fillHybridZ_left( const size_t j1, const size_t j2 );

	/**
	 * Computes all entries of the hybridE matrix for interactions starting in
	 * i1 and i2 and report all valid interactions to updateOptima()
	 *
	 * Note: (i1,i2) have to be complementary (right-most base pair of seed)
	 *
	 * @param i1 end of the interaction within seq 1
	 * @param i2 end of the interaction within seq 2
	 *
	 */
	virtual
	void
	fillHybridZ_right( const size_t i1, const size_t i2 );

	/**
	 * adds seed information and calls traceBack() of super class
	 * @param interaction IN/OUT the interaction to fill
	 */
	virtual
	void
	traceBack( Interaction & interaction );

	/**
	 * Returns the hybridization energy of the non overlapping part of seeds
	 * starting at si and sj
	 *
	 * @param si1 the index of seed1 in the first sequence
	 * @param si2 the index of seed1 in the second sequence
	 * @param sj1 the index of seed2 in the first sequence
	 * @param sj2 the index of seed2 in the second sequence
	 */
	virtual
	E_type
	getNonOverlappingEnergy( const size_t si1, const size_t si2, const size_t sj1, const size_t sj2 );

	// debug function
	void
	printMatrix( const Z2dMatrix & matrix );

};

} // namespace

#endif /* INTARNA_PREDICTORMFEENS2DSEEDEXTENSION_H_ */
