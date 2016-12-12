
#ifndef PREDICTORMFE2DHEURISTICSEED_H_
#define PREDICTORMFE2DHEURISTICSEED_H_

#include "PredictorMfe2dHeuristic.h"
#include "SeedHandler.h"


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
 * @author Martin Mann
 *
 */
class PredictorMfe2dHeuristicSeed: public PredictorMfe2dHeuristic {


	//! matrix type to hold the mfe energies and boundaries for interaction site starts
	typedef PredictorMfe2dHeuristic::E2dMatrix E2dMatrix;

public:

	/**
	 * Constructs a predictor and stores the energy and output handler
	 *
	 * @param energy the interaction energy handler
	 * @param output the output handler to report mfe interactions to
	 * @param seedConstraint the seed constraint to be applied
	 */
	PredictorMfe2dHeuristicSeed( const InteractionEnergy & energy
			, OutputHandler & output
			, const SeedConstraint & seedConstraint );

	virtual ~PredictorMfe2dHeuristicSeed();

	/**
	 * Computes the mfe for the given sequence ranges (i1-j1) in the first
	 * sequence and (i2-j2) in the second sequence and reports it to the output
	 * handler.
	 *
	 * @param r1 the index range of the first sequence interacting with r2
	 * @param r2 the index range of the second sequence interacting with r1
	 * @param reportMax the maximal number of (sub)optimal interactions to be
	 *            reported to the output handler
	 * @param reportNonOverlapping whether or not the reported interactions
	 *            should be non-overlapping or not
	 *
	 */
	virtual
	void
	predict( const IndexRange & r1 = IndexRange(0,RnaSequence::lastPos)
			, const IndexRange & r2 = IndexRange(0,RnaSequence::lastPos)
			, const size_t reportMax = 1
			, const bool reportNonOverlapping = true );

protected:

	//! access to the interaction energy handler of the super class
	using PredictorMfe2dHeuristic::energy;

	//! access to the output handler of the super class
	using PredictorMfe2dHeuristic::output;

	// TODO provide all data structures as arguments to make predict() call threadsafe

	//! energy of all interaction hybrids that end in position p (seq1) and
	//! q (seq2)
	using PredictorMfe2dHeuristic::hybridE;

	//! the best hybridization energy including a seed for start i1,i2
	E2dMatrix hybridE_seed;

	//! handler to generate and access seed information
	SeedHandler seedHandler;

protected:

	/**
	 * Fills a given interaction (boundaries given) with the according
	 * hybridizing base pairs.
	 * @param interaction IN/OUT the interaction to fill
	 */
	virtual
	void
	traceBack( Interaction & interaction );


};

#endif /* PREDICTORMFE2DHEURISTICSEED_H_ */
