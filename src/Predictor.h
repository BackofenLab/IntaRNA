
#ifndef PREDICTOR_H_
#define PREDICTOR_H_

#include "InteractionEnergy.h"

#include "OutputHandler.h"

/**
 * RNA-RNA interaction prediction handler interface.
 *
 * @author Martin Mann 2014
 *
 */
class Predictor {

public:

	/**
	 * Constructs an RNA-RNA interaction prediction handler and sets the
	 * central data members.
	 *
	 * @param energy the handler for interaction energy computation
	 * @param output the output handler for identified interactions
	 *
	 */
	Predictor( const InteractionEnergy & energy, OutputHandler & output );

	/**
	 * destruction
	 */
	virtual ~Predictor();

	/**
	 * Computes the predictors optimization target for the given sequence
	 * ranges (i1-j1) in the first sequence and (i2-j2) in the second sequence.
	 * The according optimal interaction is given to the output handler.
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1
	 *
	 */
	virtual
	void
	predict( const size_t i1 = 0, const size_t j1 = RnaSequence::lastPos
			, const size_t i2 = 0, const size_t j2 = RnaSequence::lastPos) = 0;

protected:

	//! energy computation handler
	const InteractionEnergy & energy;

	//! interaction output handler
	OutputHandler & output;

};

#endif /* PREDICTOR_H_ */
