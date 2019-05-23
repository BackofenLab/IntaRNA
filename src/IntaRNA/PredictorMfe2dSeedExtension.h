
#ifndef INTARNA_PREDICTORMFE2DSEEDEXTENSION_H_
#define INTARNA_PREDICTORMFE2DSEEDEXTENSION_H_

#include "IntaRNA/PredictorMfe2d.h"
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
class PredictorMfe2dSeedExtension: public PredictorMfe2d {

protected:

	//! matrix type to hold the mfe energies for interaction site starts
	typedef PredictorMfe2d::E2dMatrix E2dMatrix;

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
	PredictorMfe2dSeedExtension(
			const InteractionEnergy & energy
			, OutputHandler & output
			, PredictionTracker * predTracker
			, SeedHandler * seedHandler );


	/**
	 * data cleanup
	 */
	virtual ~PredictorMfe2dSeedExtension();


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
	using PredictorMfe2d::energy;

	//! access to the output handler of the super class
	using PredictorMfe2d::output;

	//! the seed handler (with idx offset)
	SeedHandlerIdxOffset seedHandler;

	//! energy of all interaction hybrids that start on the left side of the seed
	E2dMatrix hybridE_left;

	//! energy of all interaction hybrids that start on the right side of the seed
	E2dMatrix hybridE_right;

protected:

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

	/**
	 * Computes all entries of the hybridE matrix for interactions starting in
	 * i1 and i2 and report all valid interactions to updateOptima()
	 *
	 * @param j1 start of the interaction within seq 1
	 * @param j2 start of the interaction within seq 2
	 * @param outConstraint constrains the interactions reported to the output handler
	 *
	 */
	void
	fillHybridE_left( const size_t j1, const size_t j2
						, const OutputConstraint & outConstraint);

	/**
	 * Computes all entries of the hybridE matrix for interactions starting in
	 * i1 and i2 and report all valid interactions to updateOptima()
	 *
	 * @param i1 start of the interaction within seq 1
	 * @param i2 start of the interaction within seq 2
	 * @param outConstraint constrains the interactions reported to the output handler
	 *
	 */
	void
	fillHybridE_right( const size_t i1, const size_t i2
				, const OutputConstraint & outConstraint);

	/**
	 * Fills a given interaction (boundaries given) with the according
	 * hybridizing base pairs using hybridE_seed.
	 * @param interaction IN/OUT the interaction to fill
	 * @param outConstraint constrains the interactions reported to the output handler
	 */
	virtual
	void
	traceBack( Interaction & interaction, const OutputConstraint & outConstraint  );

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

	void
	printMatrix( const E2dMatrix & matrix );

};

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

inline
void
PredictorMfe2dSeedExtension::
updateOptima( const size_t i1, const size_t j1
		, const size_t i2, const size_t j2
		, const E_type energy
		, const bool isHybridE )
{
	// do nothing and ignore calls from fillHybridE()
}

//////////////////////////////////////////////////////////////////////////

inline
void
PredictorMfe2dSeedExtension::
printMatrix( const E2dMatrix & matrix )
{
	for (int i = 0; i < matrix.size1(); i++) {
		std::cout << "| ";
		for (int j = 0; j < matrix.size2(); j++) {
			std::cout << matrix(i, j) << " | ";
		}
		std::cout << std::endl;
	}
}

} // namespace

#endif /* INTARNA_PREDICTORMFE2DSEEDEXTENSION_H_ */
