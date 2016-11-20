
#ifndef PREDICTORMFE2DHEURISTIC_H_
#define PREDICTORMFE2DHEURISTIC_H_

#include "PredictorMfe.h"
#include "Interaction.h"

#include <boost/numeric/ublas/matrix.hpp>

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

	/**
	 * Describes the currently best interaction found for a left interaction
	 * boundary i1,i2
	 */
	class BestInteraction {
	public:

		//! init data
		BestInteraction( const E_type E=E_INF, const size_t j1=RnaSequence::lastPos, const size_t j2=RnaSequence::lastPos )
			: E(E), j1(j1), j2(j2)
		{}

	public:
		//! energy of the interaction
		E_type E;
		//! right end of the interaction in seq1
		size_t j1;
		//! right end of the interaction in seq2
		size_t j2;
	};

	//! matrix type to hold the mfe energies and boundaries for interaction site starts
	typedef boost::numeric::ublas::matrix<BestInteraction> E2dMatrix;

public:

	/**
	 * Constructs a predictor and stores the energy and output handler
	 *
	 * @param energy the interaction energy handler
	 * @param output the output handler to report mfe interactions to
	 */
	PredictorMfe2dHeuristic( const InteractionEnergy & energy, OutputHandler & output );

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

	//! access to the mfe interaction of the super class
	using PredictorMfe::mfeInteraction;

	// TODO provide all data structures as arguments to make predict() call threadsafe

	//! energy of all interaction hybrids starting in i1,i2
	E2dMatrix hybridE;

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
	traceBack( Interaction & interaction );

};

#endif /* PREDICTORMFE2DHEURISTIC_H_ */
