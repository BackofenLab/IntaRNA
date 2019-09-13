
#ifndef INTARNA_PREDICTORMFEENS2D_H_
#define INTARNA_PREDICTORMFEENS2D_H_

#include "IntaRNA/PredictorMfeEns.h"
#include "IntaRNA/Interaction.h"

#include <boost/numeric/ublas/matrix.hpp>

namespace IntaRNA {

/**
 * Memory efficient ensemble predictor for RNAup-like computation, i.e. full
 * DP-implementation without seed-heuristic, using 2D matrices
 *
 * @author Martin Raden
 * @author Frank Gelhausen
 *
 */
class PredictorMfeEns2d: public PredictorMfeEns {

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
	 */
	PredictorMfeEns2d( const InteractionEnergy & energy
					, OutputHandler & output
					, PredictionTracker * predTracker );

	virtual ~PredictorMfeEns2d();

	/**
	 * Computes the mfe for the given sequence ranges (i1-j1) in the first
	 * sequence and (i2-j2) in the second sequence and reports it to the output
	 * handler.
	 *
	 * @param r1 the index range of the first sequence interacting with r2
	 * @param r2 the index range of the second sequence interacting with r1
	 *
	 */
	virtual
	void
	predict( const IndexRange & r1 = IndexRange(0,RnaSequence::lastPos)
			, const IndexRange & r2 = IndexRange(0,RnaSequence::lastPos)
			);

protected:

	//! access to the interaction energy handler of the super class
	using PredictorMfeEns::energy;

	//! access to the output handler of the super class
	using PredictorMfeEns::output;

	//! energy of all interaction hybrids that end in position p (seq1) and
	//! q (seq2)
	Z2dMatrix hybridZ;

protected:

	/**
	 * Computes all entries of the hybridE matrix for interactions ending in
	 * p=j1 and q=j2 and report all valid interactions to updateZ()
	 *
	 * @param j1 end of the interaction within seq 1
	 * @param j2 end of the interaction within seq 2
	 * @param i1init smallest value for i1
	 * @param i2init smallest value for i2
	 * @param callUpdateZ whether or not updateZ() is to be called
	 *
	 */
	virtual
	void
	fillHybridZ( const size_t j1, const size_t j2
				, const size_t i1init, const size_t i2init
				, const bool callUpdateZ
				);

};

} // namespace

#endif /* INTARNA_PREDICTORMFEENS2D_H_ */
