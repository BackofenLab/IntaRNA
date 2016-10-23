
#ifndef PREDICTORMFE_H_
#define PREDICTORMFE_H_


#include "Predictor.h"

/**
 * Generic Predictor interface for MFE interaction computation to avoid
 * code redundancy
 *
 * @author Martin Mann
 *
 */
class PredictorMfe : public Predictor {


public:

	/**
	 * Construction call the super constructor
	 */
	PredictorMfe( const InteractionEnergy & energy, OutputHandler & output );

	virtual ~PredictorMfe();

protected:


	//! access to the interaction energy handler of the super class
	using Predictor::energy;

	//! access to the output handler of the super class
	using Predictor::output;

	// TODO provide all data structures as arguments to make predict() call threadsafe

	//! mfe interaction boundaries
	Interaction mfeInteraction;

	//! minimal stacking energy
	const E_type minStackingEnergy;
	//! minimal interaction initiation energy
	const E_type minInitEnergy;
	//! minimal dangling end energy
	const E_type minDangleEnergy;
	//! minimal interaction end energy
	const E_type minEndEnergy;

	/**
	 * Initializes the global energy minimum
	 */
	virtual
	void
	initMfe();

	/**
	 * updates the global optimum to be the mfe interaction if needed
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1
	 * @param hybridE the energy of the interaction only (init+loops)
	 */
	virtual
	void
	updateMfe( const size_t i1, const size_t j1
			, const size_t i2, const size_t j2
			, const E_type energy );

};

#endif /* PREDICTORMFE_H_ */
