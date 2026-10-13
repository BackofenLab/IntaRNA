
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
 * Optional pair probabilities refer to all admitted candidate chains considered
 * by this pruning rule, with the same Zall as ordinary heuristic prediction.
 * They approximate the unrestricted interaction ensemble. Collection adds
 * quadratic storage and a reverse traversal of the retained continuations.
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
	 * @param pairProbabilities optional non-owning sink for probabilities of the
	 *         pruned candidate ensemble; requires needZall
	 */
	PredictorMfeEns2dHeuristic( const InteractionEnergy & energy
							, OutputHandler & output
							, PredictionTracker * predTracker
							, BasePairProbabilities * pairProbabilities = nullptr );

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

	//! The selected chain owns its first pair and, optionally, its next stack pair.
	struct ProbabilityContinuation {
		size_t next1 = RnaSequence::lastPos, next2 = RnaSequence::lastPos;
		bool extraPair = false;
	};
	//! Allocated only for probability output; selected continuations form a DAG.
	Matrix<ProbabilityContinuation> probabilityContinuation;
	//! Complete candidate weight flowing into each retained continuation.
	Matrix<Z_type> probabilityFlow;
	//! Numerators in local target/reversed-query coordinates until region commit.
	Matrix<Z_type> probabilityMass;
	//! Independently accumulated objective, required to agree exactly with Zall.
	Z_type probabilityDenominator = 0;

protected:

	/** Compute one region under predict()'s failure guard. */
	void predictRegionHeuristic(const IndexRange & r1,const IndexRange & r2);

	/** Record precisely the accepted updateZ candidate, crediting newly prepended
	 * pairs and routing its full weight to the retained child chain. Coordinates
	 * are local and the query is reversed; lastPos denotes an initial candidate.
	 */
	void recordProbabilityCandidate(size_t i1,size_t j1,size_t i2,size_t j2,
			Z_type hybrid,const ProbabilityContinuation & continuation);

	/** Propagate candidate weights through selected chains and commit this region. */
	void commitProbabilityRegion();

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
