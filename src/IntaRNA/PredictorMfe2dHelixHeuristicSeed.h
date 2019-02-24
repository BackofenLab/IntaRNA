
#ifndef INTARNA_PREDICTORMFE2DHELIXHEURISTICSEED_H
#define INTARNA_PREDICTORMFE2DHELIXHEURISTICSEED_H

#include "IntaRNA/PredictorMfe2dHelixHeuristic.h"
#include "IntaRNA/SeedHandlerIdxOffset.h"

namespace IntaRNA {


/**
 * Memory efficient interaction predictor that uses both qualitative heuristics
 * (interactions have to have a seed interaction) and performance heuristics
 * (not all possible interactions considered)
 *
 * To this end, for each interaction start i1,i2 only the optimal right side
 * interaction with boundaries j1,j2 is considered in the recursion instead of
 * all possible interaction ranges.
 *
 * This yields a quadratic time and space complexity.
 *
 * @author Rick Gelhausen
 *
 */
class PredictorMfe2dHelixHeuristicSeed: public PredictorMfe2dHelixHeuristic {


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
	 * @param seedConstraint the seed constraint to be applied
	 */
	PredictorMfe2dHelixHeuristicSeed( const InteractionEnergy & energy
			, OutputHandler & output
			, PredictionTracker * predTracker
		 	, const HelixConstraint & helixConstraint
			, SeedHandler * seedHandlerInstance
	);

	virtual ~PredictorMfe2dHelixHeuristicSeed();

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
			, const OutputConstraint & outConstraint = OutputConstraint());

protected:

	//! access to the interaction energy handler of the super class
	using PredictorMfe2dHelixHeuristic::energy;

	//! access to the output handler of the super class
	using PredictorMfe2dHelixHeuristic::output;

	//! access to the list of reported interaction ranges of the super class
	using PredictorMfe2dHelixHeuristic::reportedInteractions;

	//! energy of all interaction hybrids that end in position p (seq1) and
	//! q (seq2)
	using PredictorMfe2dHelixHeuristic::hybridE;

	//! handler to generate and access helix information with idx offset
	using PredictorMfe2dHelixHeuristic::helixHandler;

	//! the best hybridization energy including a seed for start i1,i2
	E2dMatrix hybridE_seed;

	//
	SeedHandlerIdxOffset seedHandler;

protected:

	/**
	 * Fills a given interaction (boundaries given) with the according
	 * hybridizing base pairs.
	 * @param interaction IN/OUT the interaction to fill
	 */
	virtual
	void
	traceBack( Interaction & interaction, const OutputConstraint & outConstraint );


	/**
	 * Identifies the next best interaction (containing a seed)
	 * with an energy equal to or higher
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

	/**
	 * does nothing but to ignore the calls from fillHybridE()
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1
	 * @param energy ignored
	 * @param isHybridE ignored
	 */
	virtual
	void
	updateOptima( const size_t i1, const size_t j1
			, const size_t i2, const size_t j2
			, const E_type energy
			, const bool isHybridE );


};

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

inline
void
PredictorMfe2dHelixHeuristicSeed::
updateOptima( const size_t i1, const size_t j1
		, const size_t i2, const size_t j2
		, const E_type energy
		, const bool isHybridE )
{
	// do nothing and ignore calls from fillHybridE()
}

//////////////////////////////////////////////////////////////////////////

} // namespace
#endif //INTARNA_PREDICTORMFE2DHELIXHEURISTICSEED_H
