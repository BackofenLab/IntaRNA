
#ifndef INTARNA_PREDICTORMFE2DHELIXHEURISTIC_H_
#define INTARNA_PREDICTORMFE2DHELIXHEURISTIC_H_

#include "IntaRNA/PredictorMfe2dHeuristic.h"
#include "IntaRNA/Interaction.h"
#include "IntaRNA/HelixHandlerIdxOffset.h"

#include <boost/numeric/ublas/matrix.hpp>

namespace IntaRNA {

/**
 * Memory efficient interaction predictor that uses a heuristic to
 * find the mfe or a close-to-mfe interaction.
 *
 * To this end, for each interaction start i1,i2 only the optimal right side
 * interaction with boundaries j1,j2 is considered in the recursion instead of
 * all possible interaction ranges.
 *
 * Helices (continues stackings of base pairs) are restricted to a given length.
 *
 * This yields a quadratic time and space complexity.
 *
 * @author Rick Gelhausen
 *
 */
class PredictorMfe2dHelixHeuristic: public PredictorMfe2dHeuristic {

protected:

	//! matrix type to hold the mfe energies and boundaries for interaction site starts
	typedef PredictorMfe2dHeuristic::E2dMatrix E2dMatrix;

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
	PredictorMfe2dHelixHeuristic( const InteractionEnergy & energy
			, OutputHandler & output
			, PredictionTracker * predTracker
			, const HelixConstraint & helixConstraint);

	virtual ~PredictorMfe2dHelixHeuristic();

	/**
	 * Computes the mfe for the given sequence ranges (i1-j1) in the first
	 * sequence and (i2-j2) in the second sequence and reports it to the output
	 * handler.
	 *
	 * @param r1 the index range of the first sequence interacting with r2
	 * @param r2 the index range of the second sequence interacting with r1
	 * @param outConstraint constrains the interactions reported to the output handler
	 *
	 */
	virtual
	void
	predict( const IndexRange & r1 = IndexRange(0,RnaSequence::lastPos)
			, const IndexRange & r2 = IndexRange(0,RnaSequence::lastPos)
			, const OutputConstraint & outConstraint = OutputConstraint() );

protected:

	//! access to the interaction energy handler of the super class
	using PredictorMfe2dHeuristic::energy;

	//! access to the output handler of the super class
	using PredictorMfe2dHeuristic::output;

	//! access to the list of reported interaction ranges of the super class
	using PredictorMfe2dHeuristic::reportedInteractions;

	//! energy of all interaction hybrids starting in i1,i2
	using PredictorMfe2dHeuristic::hybridE;

	//! helixHandler used only for creating HelixHandlerOffset
	//HelixHandler helixHandler;
	HelixHandlerIdxOffset helixHandler;

protected:

	/**
	 * Computes all entries of the hybridE matrix
	 */
	virtual
	void
	fillHybridE();

	/**
	 * Fills a given interaction (boundaries given) with the according
	 * hybridizing base pairs.
	 * @param interaction IN/OUT the interaction to fill
	 */
	virtual
	void
	traceBack( Interaction & interaction, const OutputConstraint & outConstraint );

	/**
	 * Identifies the next best interaction with an energy equal to or higher
	 * than the given interaction. The new interaction will not overlap any
	 * index range stored in reportedInteractions.
	 *
	 * @param curBest IN/OUT the current best interaction to be replaced with one
	 *        of equal or higher energy not overlapping with any reported
	 *        interaction so far; an interaction with energy E_INF is set, if
	 *        there is no better interaction left
	 */
	virtual
	void
	getNextBest( Interaction & curBest );
};

} // namespace

#endif /* INTARNA_PREDICTORMFE2DHELIXHEURISTIC_H_ */
