
#ifndef INTARNA_PREDICTORMFE2DHEURISTIC_H_
#define INTARNA_PREDICTORMFE2DHEURISTIC_H_

#include "IntaRNA/PredictorMfe.h"
#include "IntaRNA/Interaction.h"

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
 * This yields a quadratic time and space complexity.
 *
 * @author Martin Mann
 *
 */
class PredictorMfe2dHeuristic: public PredictorMfe {

protected:

	//! matrix type to hold the mfe energies and boundaries for interaction site starts
	typedef boost::numeric::ublas::matrix<BestInteractionE> E2dMatrix;

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
	PredictorMfe2dHeuristic( const InteractionEnergy & energy
							, OutputHandler & output
							, PredictionTracker * predTracker );

	virtual ~PredictorMfe2dHeuristic();

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
	using PredictorMfe::energy;

	//! access to the output handler of the super class
	using PredictorMfe::output;

	//! access to the list of reported interaction ranges of the super class
	using PredictorMfe::reportedInteractions;

	//! energy of all interaction hybrids starting in i1,i2
	E2dMatrix hybridE;

protected:

	/**
	 * Computes all entries of the hybridE matrix
	 * and reports all valid interactions via updateOptima()
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
	traceBack( Interaction & interaction );

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

	/**
	 * Overwrites function of super class to surpress the update.
	 *
	 * @param i1 interaction start in seq1
	 * @param j1 interaction end in seq1
	 * @param i2 interaction start in seq2
	 * @param i2 interaction end in seq2
	 * @param curInteraction the interaction information to be used for update
	 */
	virtual
	void
	updateMfe4leftEnd(const size_t i1, const size_t j1
					, const size_t i2, const size_t j2
					, const Interaction & curInteraction );

};

} // namespace

#endif /* INTARNA_PREDICTORMFE2DHEURISTIC_H_ */
