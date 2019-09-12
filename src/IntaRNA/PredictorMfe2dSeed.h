/*
 * PredictorMfe2dSeed.h
 *
 *  Created on: 13.10.2016
 *      Author: Mmann
 */

#ifndef INTARNA_PREDICTORMFE2DSEED_H_
#define INTARNA_PREDICTORMFE2DSEED_H_

#include "IntaRNA/PredictorMfe2d.h"
#include "IntaRNA/SeedHandlerIdxOffset.h"

namespace IntaRNA {

/**
 * Implements seed-based space-efficient interaction prediction.
 *
 * Note, for each seed start (i1,i2) only the mfe seed is considered for the
 * overall interaction computation instead of considering all possible seeds
 * starting at (i1,i2).
 *
 * @author Martin Mann
 *
 */
class PredictorMfe2dSeed: public PredictorMfe2d {

protected:

	//! matrix type to hold the mfe energies for interaction site starts
	typedef PredictorMfe2d::E2dMatrix E2dMatrix;

public:


	/**
	 * Constructs a predictor and stores the energy and output handler
	 *
	 * @param energy the interaction energy handler
	 * @param output the output handler to report mfe interactions to
	 * @param predTracker the prediction tracker to be used or NULL if no
	 *         tracking is to be done; if non-NULL, the tracker gets deleted
	 *         on this->destruction.
	 * @param seedHandler the seed handler to be used
	 */
	PredictorMfe2dSeed(
			const InteractionEnergy & energy
			, OutputHandler & output
			, PredictionTracker * predTracker
			, SeedHandler * seedHandler );


	/**
	 * data cleanup
	 */
	virtual ~PredictorMfe2dSeed();


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


	//! access to the interaction energy handler of the super class
	using PredictorMfe2d::energy;

	//! access to the output handler of the super class
	using PredictorMfe2d::output;

	//! energy of all interaction hybrids that end in position p (seq1) and
	//! q (seq2) and do not necessarily contain a seed interaction
	using PredictorMfe2d::hybridE_pq;

	//! the seed handler (with idx offset)
	SeedHandlerIdxOffset seedHandler;

	//! for fixed interaction end p=j1,q=j2: each cell (i1,i2) provides the mfe
	//! for the interaction i1..j1 with i2..j2 given that the range contains
	//! a valid seed interaction
	E2dMatrix hybridE_pq_seed;

protected:

	/**
	 * Computes all entries of both hybridE matrices for interactions ending in
	 * p=j1 and q=j2 and reports all valid interactions via updateOptima()
	 *
	 * @param j1 end of the interaction within seq 1
	 * @param j2 end of the interaction within seq 2
	 * @param i1init smallest value for i1
	 * @param i2init smallest value for i2
	 * @param callUpdateOptima whether or not to call updateOptima()
	 */
	virtual
	void
	fillHybridE( const size_t j1, const size_t j2
				, const size_t i1init, const size_t i2init
				, const bool callUpdateOptima
				);

	/**
	 * Fills a given interaction (boundaries given) with the according
	 * hybridizing base pairs using hybridE_seed.
	 * @param interaction IN/OUT the interaction to fill
	 */
	virtual
	void
	traceBack( Interaction & interaction );

};

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

} // namespace

#endif /* INTARNA_PREDICTORMFE2DSEED_H_ */
