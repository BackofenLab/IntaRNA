
#ifndef INTARNA_PREDICTORMFEENS2DSEEDEXTENSION_H_
#define INTARNA_PREDICTORMFEENS2DSEEDEXTENSION_H_

#include "IntaRNA/PredictorMfeEns.h"
#include "IntaRNA/SeedHandlerIdxOffset.h"

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
	typedef boost::numeric::ublas::matrix<Z_type> Z2dMatrix;

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
	PredictorMfeEns2dSeedExtension(
			const InteractionEnergy & energy
			, OutputHandler & output
			, PredictionTracker * predTracker
			, SeedHandler * seedHandler );


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
