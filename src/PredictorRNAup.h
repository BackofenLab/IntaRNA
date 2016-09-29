
#ifndef PREDICTORRNAUP_H_
#define PREDICTORRNAUP_H_

#include "Predictor.h"

#include <boost/numeric/ublas/banded.hpp>

/**
 * Predictor for RNAup-like computation, i.e. full DP-implementation without
 * seed-heuristic
 *
 * @author Martin Mann
 *
 */
class PredictorRNAup: public Predictor {

protected:

	//! matrix type to cover the energies for different interaction site widths
	typedef boost::numeric::ublas::matrix<E_type> E2dMatrix;

	//! full 4D DP-matrix for computation to hold all start position combinations
	//! first index = start positions (i1,i2) of (seq1,seq2)
	//! second index = interaction window sizes (w1,w2) or NULL if (i1,i2) not complementary
	typedef boost::numeric::ublas::matrix< E2dMatrix* > E4dMatrix;


public:


	/**
	 * construction
	 *
	 * @param energy the interaction energy handler
	 * @param output the output handler
	 */
	PredictorRNAup( const InteractionEnergy & energy, const OutputHandler & output );

	virtual ~PredictorRNAup();


protected:

	//! access to the interaction energy handler of the super class
	using Predictor::energy;

	//! access to the output handler of the super class
	using Predictor::output;

	//! energy of all interaction hybrids computed by the recursion with indices
	//! hybridE(i1,i2)->(w1,w2), with interaction start i1 (seq1) and i2 (seq2) and
	//! ineraction end j1=i1+w1 and j2=j2+w2
	//! NOTE: hybridE(i1,i2)==NULL if not complementary(seq1[i1],seq2[i2])
	E4dMatrix hybridE;

	//! minimal value within hybridE
	E_type hybridEmin;

protected:


	/**
	 * computes all entries of the hybridE matrix
	 */
	void
	fillHybridE( const InteractionEnergy & energy );


};

#endif /* PREDICTORRNAUP_H_ */
