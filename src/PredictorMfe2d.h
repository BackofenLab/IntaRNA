
#ifndef PREDICTORMFE2D_H_
#define PREDICTORMFE2D_H_

#include "PredictorMfe.h"
#include "Interaction.h"

#include <boost/numeric/ublas/matrix.hpp>

/**
 * Memory efficient predictor for RNAup-like computation, i.e. full
 * DP-implementation without seed-heuristic, using 2D matrices
 *
 * @author Martin Mann
 *
 */
class PredictorMfe2d: public PredictorMfe {

protected:

	//! matrix type to hold the mfe energies for interaction site starts
	typedef boost::numeric::ublas::matrix<E_type> E2dMatrix;

public:

	/**
	 * Constructs a predictor and stores the energy and output handler
	 *
	 * @param energy the interaction energy handler
	 * @param output the output handler to report mfe interactions to
	 */
	PredictorMfe2d( const InteractionEnergy & energy, OutputHandler & output );

	virtual ~PredictorMfe2d();

	/**
	 * Computes the mfe for the given sequence ranges (i1-j1) in the first
	 * sequence and (i2-j2) in the second sequence and reports it to the output
	 * handler.
	 *
	 * @param r1 the index range of the first sequence interacting with r2
	 * @param 22 the index range of the second sequence interacting with r1
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

	//! access to the index offset in seq1 of the super class
	using PredictorMfe::i1offset;

	//! access to the index offset in seq2 of the super class
	using PredictorMfe::i2offset;

	// TODO provide all data structures as arguments to make predict() call threadsafe

	//! energy of all interaction hybrids that end in position p (seq1) and
	//! q (seq2)
	E2dMatrix hybridE_pq;

protected:

	/**
	 * Initializes the hybridE_pq table for the computation for interactions
	 * ending in p=j1 and q=j2
	 *
	 * @param energy the energy function to use
	 * @param j1 end of the interaction within seq 1
	 * @param j2 end of the interaction within seq 2
	 */
	void
	initHybridE( const size_t j1, const size_t j2 );

	/**
	 * Computes all entries of the hybridE matrix for interactions ending in
	 * p=j1 and q=j2
	 *
	 * @param energy the energy function to use
	 * @param j1 end of the interaction within seq 1
	 * @param j2 end of the interaction within seq 2
	 *
	 */
	void
	fillHybridE( const size_t j1, const size_t j2  );

	/**
	 * Fills a given interaction (boundaries given) with the according
	 * hybridizing base pairs.
	 * @param interaction IN/OUT the interaction to fill
	 */
	void
	traceBack( Interaction & interaction );

};

#endif /* PREDICTORMFE2D_H_ */
