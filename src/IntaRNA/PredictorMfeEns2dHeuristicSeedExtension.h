
#ifndef INTARNA_PREDICTORMFEENS2DHEURISTICSEEDEXTENSION_H_
#define INTARNA_PREDICTORMFEENS2DHEURISTICSEEDEXTENSION_H_

#include "IntaRNA/PredictorMfeEns2dSeedExtension.h"
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
	 * @param seedHandler the seed handler to be used
	 */
	PredictorMfeEns2dHeuristicSeedExtension(
			const InteractionEnergy & energy
			, OutputHandler & output
			, PredictionTracker * predTracker
			, SeedHandler * seedHandler );


	/**
	 * data cleanup
	 */
	virtual ~PredictorMfeEns2dHeuristicSeedExtension();


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
	 * @param outConstraint constrains the interactions reported to the output handler
	 *
	 */
	virtual
	void
	fillHybridZ_right( const size_t i1, const size_t i2
				, const OutputConstraint & outConstraint
				, const size_t si1, const size_t si2
				);

	/**
	 * Computes all entries of the hybridE matrix for interactions ending in
	 * p=j1 and q=j2 and report all valid interactions to updateOptima()
	 *
   * @param j1 end of the interaction within seq 1
	 * @param j2 end of the interaction within seq 2
	 * @param outConstraint constrains the interactions reported to the output handler
	 *
	 */
	virtual
	void
	fillHybridZ_left( const size_t j1, const size_t j2
				, const OutputConstraint & outConstraint );

//	/**
//	 * Fills a given interaction (boundaries given) with the according
//	 * hybridizing base pairs using hybridE_seed.
//	 * @param interaction IN/OUT the interaction to fill
//	 * @param outConstraint constrains the interactions reported to the output handler
//	 */
//	virtual
//	void
//	traceBack( Interaction & interaction, const OutputConstraint & outConstraint  );

	/**
	 * Identifies the next best interaction with an energy equal to or higher
	 * than the given interaction. The new interaction will not overlap any
	 * index range stored in reportedInteractions.
	 *
	 * NOTE: this is not possible for this predictor (unless a full recomputation
	 * of the matrices is done). Thus, calling this method raises an exception.
	 *
	 * @param curBest ignored (see method comment)
	 */
	virtual
	void
	getNextBest( Interaction & curBest );

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
	if (fullE < E_right_opt) {
		E_right_opt = fullE;
		j1opt = j1;
		j2opt = j2;
	}
}

//////////////////////////////////////////////////////////////////////////

} // namespace

#endif /* INTARNA_PREDICTORMFEENS2DHEURISTICSEEDEXTENSION_H_ */
