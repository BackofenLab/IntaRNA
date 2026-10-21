
#ifndef INTARNA_PREDICTORMFEENS2DHEURISTICSEEDEXTENSION_H_
#define INTARNA_PREDICTORMFEENS2DHEURISTICSEEDEXTENSION_H_

#include "IntaRNA/PredictorMfeEns2dSeedExtension.h"
#include "IntaRNA/SeedHandlerIdxOffset.h"
#include <unordered_map>

namespace IntaRNA {

/**
 * Implements heuristic, seed-based, space-efficient interaction prediction
 * based on minimizing ensemble free energy of interaction sites.
 *
 * Note, for each seed start (i1,i2) only the mfe seed is considered for the
 * overall interaction computation instead of considering all possible seeds
 * starting at (i1,i2).
 * Stack-only seeds share the exact predictor's left-most anchor tables. The
 * first completed seed uniquely owns each chain, equivalently the left-most
 * seed after removing containing occurrences.
 * Retain its one-sided extensions and the left extensions ending at its
 * preferred right endpoint. Pair probabilities describe this retained ensemble.
 *
 * @author Frank Gelhausen
 * @author Martin Raden
 *
 */
class PredictorMfeEns2dHeuristicSeedExtension: public PredictorMfeEns2dSeedExtension {

public:


	/**
	 * Constructs a predictor and stores the energy and output handler
	 *
	 * @param energy the interaction energy handler
	 * @param output the output handler to report mfe interactions to
	 * @param predTracker the prediction tracker to be used or NULL if no
	 *         tracking is to be done; if non-NULL, the tracker gets deleted
	 *         on this->destruction.
	 * @param pairProbabilities optional non-owning sink for this approximate ensemble
	 * @param seedHandler the seed handler to be used
	 */
	PredictorMfeEns2dHeuristicSeedExtension(
			const InteractionEnergy & energy
			, OutputHandler & output
			, PredictionTracker * predTracker
			, SeedHandler * seedHandler
			, BasePairProbabilities * pairProbabilities = nullptr );


	/**
	 * data cleanup
	 */
	virtual ~PredictorMfeEns2dHeuristicSeedExtension() override;


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
			, const IndexRange & r2 = IndexRange(0,RnaSequence::lastPos) ) override;


protected:

	/** Freeze one preferred right endpoint for each admitted stack anchor.
	 * The existing virtual right-extension hooks determine the choice; partition
	 * contributions are deferred to the shared left-most-anchor recurrence.
	 * @param n first active range length
	 * @param m second active range length
	 * @param seeds contained stack occurrences in the active range
	 */
	void prepareStackSeedExtensions(size_t n,size_t m,const StackSeedDomain & seeds) override;
	/** @return the frozen right-endpoint selector for the approximate ensemble */
	std::function<Interaction::BasePair(Interaction::BasePair,Interaction::BasePair)>
	stackSeedRightSelector() const override;

	//! Only rank right extensions while preparing canonical stack-seed choices.
	bool selectingStackSeedRight = false;
	//! Preferred local endpoint keyed by the local first pair of each seed.
	std::unordered_map<Interaction::BasePair,Interaction::BasePair,Interaction::BasePair::Hash> stackSeedRight;


	//! access to the interaction energy handler of the super class
	using PredictorMfeEns2dSeedExtension::energy;

	//! access to the output handler of the super class
	using PredictorMfeEns2dSeedExtension::output;

	//! the seed handler (with idx offset)
	using PredictorMfeEns2dSeedExtension::seedHandler;

	//! energy of all interaction hybrids that end in position p (seq1) and
	//! q (seq2) and do not necessarily contain a seed interaction
	using PredictorMfeEns2dSeedExtension::hybridZ_left;

	//! energy of all interaction hybrids that start in position p (seq1) and
	//! q (seq2)
	using PredictorMfeEns2dSeedExtension::hybridZ_right;

	//! optimal energy of the right extension
	Z_type E_right_opt = E_INF;

	//! boundaries of the right extension with the optimal energy
	size_t j1opt, j2opt;

protected:

	/**
	 * updates E_right_opt, j1opt, and j2opt if needed
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1
	 * @param hybridE ensemble energy of the hybridization
	 */
	virtual
	void
	updateOptRightZ( const size_t i1, const size_t j1
			, const size_t i2, const size_t j2
			, const E_type hybridE );

	/**
	 * Computes all entries of the hybridE matrix for interactions starting in
	 * i1 and i2 and report all valid interactions to updateOptima()
	 *
	 * @param i1 end of the interaction within seq 1
	 * @param i2 end of the interaction within seq 2
	 *
	 */
	virtual
	void
	fillHybridZ_right( const size_t i1, const size_t i2
				, const size_t si1, const size_t si2
				);

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
	fillHybridZ_left( const size_t j1, const size_t j2 ) override;

};

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

inline
void
PredictorMfeEns2dHeuristicSeedExtension::
updateOptRightZ( const size_t i1, const size_t j1
		, const size_t i2, const size_t j2
		, const E_type hybridE )
{
	// store boundaries and energy of the optimal right extension
	E_type fullE = energy.getE( i1, j1, i2, j2, hybridE );
	if (fullE < E_right_opt
			&& (!output.getOutputConstraint().noGUend || !energy.isGU(j1,j2)))
	{
		E_right_opt = fullE;
		j1opt = j1;
		j2opt = j2;
	}
}

//////////////////////////////////////////////////////////////////////////

} // namespace

#endif /* INTARNA_PREDICTORMFEENS2DHEURISTICSEEDEXTENSION_H_ */
