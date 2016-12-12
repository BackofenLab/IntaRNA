
#ifndef PREDICTORMAXPROB_H_
#define PREDICTORMAXPROB_H_

#include "Predictor.h"
#include "InteractionRange.h"

#include <boost/numeric/ublas/matrix.hpp>

/**
 * Computes the interaction site with maximal probability among all interaction
 * sites
 *
 * @author Martin Mann
 *
 */
class PredictorMaxProb: public Predictor {

protected:

	//! matrix type to cover the energies for different interaction site widths
	typedef boost::numeric::ublas::matrix<E_type> E2dMatrix;

	//! full 4D DP-matrix for computation to hold all start position combinations
	//! first index = start positions (i1,i2) of (seq1,seq2)
	//! second index = interaction window sizes (w1,w2) or NULL if (i1,i2) not complementary
	typedef boost::numeric::ublas::matrix< E2dMatrix* > E4dMatrix;


public:

	/**
	 * Constructs a predictor and stores the energy and output handler
	 *
	 * @param energy the interaction energy handler
	 * @param output the output handler to report optimal interactions to
	 */
	PredictorMaxProb( const InteractionEnergy & energy, OutputHandler & output );

	virtual ~PredictorMaxProb();

	/**
	 * Computes the interaction site with maximal probability
	 * for the given sequence ranges (i1-j1) in the first
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
	using Predictor::energy;

	//! access to the output handler of the super class
	using Predictor::output;

	//! partition function of all interaction hybrids computed by the recursion with indices
	//! hybridZ(i1,i2)->(w1,w2), with interaction start i1 (seq1) and i2 (seq2) and
	//! ineraction end j1=i1+w1 and j2=j2+w2
	//! NOTE: hybridZ(i1,i2)==NULL if not complementary(seq1[i1],seq2[i2])
	E4dMatrix hybridZ;

	//! the overall partition function = sum or all hybridZ entries
	double Z;

	//! interaction boundaries with maximal probability
	InteractionRange maxProbInteraction;

protected:

	/**
	 * Removes all temporary data structures and resets the predictor
	 */
	void
	clear();

	/**
	 * computes all entries of the hybridZ matrix
	 */
	void
	fillHybridZ( );

	/**
	 * Initializes the interaction site with maximal probability
	 */
	virtual
	void
	initMaxProbInteraction();

	/**
	 * updates the global optimum if needed
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1
	 * @param hybridZ partition function for the interaction only (init+loops)
	 */
	virtual
	void
	updateMaxProbInteraction( const size_t i1, const size_t j1
			, const size_t i2, const size_t j2
			, const E_type hybridZ );

};

#endif /* PREDICTORMAXPROB_H_ */
