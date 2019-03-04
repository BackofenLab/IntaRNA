
#ifndef INTARNA_PREDICTORMAXPROB_H_
#define INTARNA_PREDICTORMAXPROB_H_

#include "IntaRNA/Predictor.h"
#include "IntaRNA/InteractionRange.h"

#include <boost/numeric/ublas/matrix.hpp>

namespace IntaRNA {

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
	typedef boost::numeric::ublas::matrix<Z_type> Z2dMatrix;

	//! full 4D DP-matrix for computation to hold all start position combinations
	//! first index = start positions (i1,i2) of (seq1,seq2)
	//! second index = interaction window sizes (w1,w2) or NULL if (i1,i2) not complementary
	typedef boost::numeric::ublas::matrix< Z2dMatrix* > Z4dMatrix;


public:

	/**
	 * Constructs a predictor and stores the energy and output handler
	 *
	 * @param energy the interaction energy handler
	 * @param output the output handler to report optimal interactions to
	 * @param predTracker the prediction tracker to be used or NULL if no
	 *         tracking is to be done; if non-NULL, the tracker gets deleted
	 *         on this->destruction.
	 */
	PredictorMaxProb( const InteractionEnergy & energy
					, OutputHandler & output
					, PredictionTracker * predTracker );

	virtual ~PredictorMaxProb();

	/**
	 * Computes the interaction site with maximal probability
	 * for the given sequence ranges (i1-j1) in the first
	 * sequence and (i2-j2) in the second sequence and reports it to the output
	 * handler.
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
	using Predictor::energy;

	//! access to the output handler of the super class
	using Predictor::output;

	//! partition function of all interaction hybrids computed by the recursion with indices
	//! hybridZ(i1,i2)->(w1,w2), with interaction start i1 (seq1) and i2 (seq2) and
	//! ineraction end j1=i1+w1 and j2=j2+w2
	//! NOTE: hybridZ(i1,i2)==NULL if not complementary(seq1[i1],seq2[i2])
	Z4dMatrix hybridZ;

	//! the overall partition function = sum or all hybridZ entries
	Z_type Z;

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
	 *
	 * @param outConstraint constrains the interactions reported to the output handler
	 */
	virtual
	void
	initOptima( const OutputConstraint & outConstraint );

	/**
	 * updates the global optimum if needed
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1
	 * @param Z partition function for the interaction
	 * @param isHybridZ whether or not the given hybridZ is only for
	 *        hybridizations (init+loops) or the total interaction energy details
	 */
	virtual
	void
	updateOptima( const size_t i1, const size_t j1
			, const size_t i2, const size_t j2
			, const Z_type Z
			, const bool isHybridZ );


	/**
	 * obsolete .. not used but needed for current interface
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1
	 * @param energy the energy of the interaction site
	 * @param isHybridE whether or not the given energy is only the
	 *        hybridization energy (init+loops) or the total interaction energy
	 */
	virtual
	void
	updateOptima( const size_t i1, const size_t j1
				, const size_t i2, const size_t j2
				, const E_type energy
				, const bool isHybridE );


	/**
	 * Pushes the stored optimal and suboptimal solutions to the output handler.
	 *
	 * @param outConstraint constrains the interactions reported to the output handler
	 */
	virtual
	void
	reportOptima( const OutputConstraint & outConstraint );


};

} // namespace

#endif /* INTARNA_PREDICTORMAXPROB_H_ */
