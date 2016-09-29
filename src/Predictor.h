
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
	Predictor( const InteractionEnergy & energy, const OutputHandler & output );

	/**
	 * destruction
	 */
	virtual ~Predictor();

protected:

	//! energy computation handler
	const InteractionEnergy & energy;

	//! interaction output handler
	const OutputHandler & output;

};

#endif /* PREDICTOR_H_ */
