
#ifndef INTARNA_PREDICTORMFEENS2DHEURISTIC_H_
#define INTARNA_PREDICTORMFEENS2DHEURISTIC_H_

#include "IntaRNA/PredictorMfeEns2d.h"
#include "IntaRNA/Interaction.h"

#include "IntaRNA/Matrix.h"

namespace IntaRNA {

/**
 * Memory efficient ensemble interaction predictor that uses a heuristic to
 * find the mfe or a close-to-mfe interaction.
 *
 * To this end, for each interaction start i1,i2 only the optimal right side
 * interaction with boundaries j1,j2 is considered in the recursion instead of
 * all possible interaction ranges.
 *
 * This yields a quadratic time and space complexity.
 *
 * @author Martin Raden
 * @author Frank Gelhausen
 *
 */
class PredictorMfeEns2dHeuristic: public PredictorMfeEns2d {

protected:

	//! matrix type to hold the mfe energies and boundaries for interaction site starts
	typedef Matrix<BestInteractionZ> Z2dMatrix;

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
	PredictorMfeEns2dHeuristic( const InteractionEnergy & energy
							, OutputHandler & output
							, PredictionTracker * predTracker );

	virtual ~PredictorMfeEns2dHeuristic();

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
			, const IndexRange & r2 = IndexRange(0,RnaSequence::lastPos) );

protected:

	//! access to the interaction energy handler of the super class
	using PredictorMfeEns2d::energy;

	//! access to the output handler of the super class
	using PredictorMfeEns2d::output;

	//! access to the list of reported interaction ranges of the super class
	using PredictorMfeEns2d::reportedInteractions;

	//! energy of all interaction hybrids starting in i1,i2
	Z2dMatrix hybridZ;

protected:

	/**
	 * Computes all entries of the hybridE matrix
	 * and reports all valid interactions via updateOptima()
	 */
	virtual
	void
	fillHybridZ();

	// Restricted output uses PredictorMfe's best finalized site per left boundary.
	// Intermediate hybridZ cells do not contain complete site ensemble energies.

};

} // namespace

#endif /* INTARNA_PREDICTORMFEENS2DHEURISTIC_H_ */
