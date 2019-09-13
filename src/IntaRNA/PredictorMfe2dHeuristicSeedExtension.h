
#ifndef INTARNA_PREDICTORMFE2DHEURISTICSEEDEXTENSION_H_
#define INTARNA_PREDICTORMFE2DHEURISTICSEEDEXTENSION_H_

#include "IntaRNA/PredictorMfe2dSeedExtension.h"
#include "IntaRNA/SeedHandlerIdxOffset.h"

namespace IntaRNA {

/**
 * Implements heuristic, seed-based, space-efficient interaction prediction.
 *
 * Note, for each seed start (i1,i2) only the mfe seed is considered for the
 * overall interaction computation instead of considering all possible seeds
 * starting at (i1,i2).
 *
 * @author Frank Gelhausen
 * @author Martin Raden
 *
 */
class PredictorMfe2dHeuristicSeedExtension: public PredictorMfe2dSeedExtension {

protected:

	//! matrix type to hold the mfe energies for interaction site starts
	typedef PredictorMfe2dSeedExtension::E2dMatrix E2dMatrix;

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
	PredictorMfe2dHeuristicSeedExtension(
			const InteractionEnergy & energy
			, OutputHandler & output
			, PredictionTracker * predTracker
			, SeedHandler * seedHandler );


	/**
	 * data cleanup
	 */
	virtual ~PredictorMfe2dHeuristicSeedExtension();


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
	using PredictorMfe2dSeedExtension::energy;

	//! access to the output handler of the super class
	using PredictorMfe2dSeedExtension::output;

	//! energy of all interaction hybrids that end in position p (seq1) and
	//! q (seq2) and do not necessarily contain a seed interaction
	using PredictorMfe2dSeedExtension::hybridE_left;

	//! the seed handler (with idx offset)
	using PredictorMfe2dSeedExtension::seedHandler;

	//! energy of all interaction hybrids that start in position p (seq1) and
	//! q (seq2)
	using PredictorMfe2dSeedExtension::hybridE_right;

  //! optimal energy of the right extension
	E_type E_right_opt;

	//! boundaries of the right extension with the optimal energy
	size_t j1opt, j2opt;

protected:

	/**
	 * updates j1opt/j2opt
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
	updateOptimalRightExt( const size_t i1, const size_t j1
			, const size_t i2, const size_t j2
			, const E_type energy
			, const bool isHybridE );

	/**
	 * Computes all entries of the hybridE matrix for interactions starting in
	 * i1 and i2 and report all valid interactions to updateOptima()
	 *
	 * @param i1 start of the interaction within seq 1
	 * @param i2 start of the interaction within seq 2
	 *
	 */
	void
	fillHybridE_right( const size_t i1, const size_t i2
				, const size_t si1, const size_t si2 );

	/**
	 * Computes all entries of the hybridE matrix for interactions ending in
	 * p=j1 and q=j2 and report all valid interactions to updateOptima()
	 *
     * @param j1 start of the interaction within seq 1
	 * @param j2 start of the interaction within seq 2
	 *
	 */
	void
	fillHybridE_left( const size_t j1, const size_t j2 );

	void
	printMatrix( const E2dMatrix & matrix );

};

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

inline
void
PredictorMfe2dHeuristicSeedExtension::
updateOptimalRightExt( const size_t i1, const size_t j1
		, const size_t i2, const size_t j2
		, const E_type energy
		, const bool isHybridE )
{
	E_type fullE = isHybridE ? this->energy.getE(i1,j1,i2,j2,energy) : energy;
	// if energy lower
	// if valid right boundary
	if (fullE < E_right_opt
			&& (!output.getOutputConstraint().noGUend || !this->energy.isGU(j1,j2)))
	{
		// store boundaries and energy of the optimal right extension
		E_right_opt = fullE;
		j1opt = j1;
		j2opt = j2;
	}
}

//////////////////////////////////////////////////////////////////////////

} // namespace

#endif /* INTARNA_PREDICTORMFE2DHEURISTICSEEDEXTENSION_H_ */
