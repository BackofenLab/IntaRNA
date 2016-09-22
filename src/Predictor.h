
#ifndef PREDICTOR_H_
#define PREDICTOR_H_

#include "Energy.h"

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
	 */
	Predictor( const Energy & energy, const OutputHandler & output );

	/**
	 * destruction
	 */
	virtual ~Predictor();

protected:

	//! energy computation handler
	const Energy & energy;

	//! interaction output handler
	const OutputHandler & output;

};

#endif /* PREDICTOR_H_ */
