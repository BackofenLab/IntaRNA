
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

	//! type for a 2D DP-matrix element
	//! (upper triangular matrix banded by acc.maxLength)
	typedef boost::numeric::ublas::banded_matrix<E_type> E2dMatrix;

	//! type of the full 4D DP-matrix for computation (banded by acc.maxLength)
	//! first index = (i,j) of query
	//! second index = (k,l) of target
	typedef boost::numeric::ublas::banded_matrix< E2dMatrix* > E4dMatrix;


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

	//! energy of all interaction hybrids computed by the recursion
	E4dMatrix hybridE;

protected:


	/**
	 * computes all entries of the hybridE matrix
	 */
	void
	fillHybridE( const InteractionEnergy & energy );


};

#endif /* PREDICTORRNAUP_H_ */
