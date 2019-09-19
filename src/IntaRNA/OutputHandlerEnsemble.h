
#ifndef INTARNA_OUTPUTHANDLERENSEMBLE_H_
#define INTARNA_OUTPUTHANDLERENSEMBLE_H_

#include "IntaRNA/OutputConstraint.h"
#include "IntaRNA/OutputHandler.h"
#include "IntaRNA/InteractionEnergy.h"

#include <iostream>

namespace IntaRNA {

/**
 * Interaction output handler that writes ensemble information to stream
 * in a simple linewise, key-value format.
 *
 * @author Martin Mann 2014
 */
class OutputHandlerEnsemble: public OutputHandler {
public:

	/**
	 * Construct a ensemble output handler.
	 *
	 * @param outConstraint the output constraint applied to find the reported
	 *        interaction
	 * @param out the stream to write to
	 * @param energy the interaction energy object used for computation
	 */
	OutputHandlerEnsemble( const OutputConstraint & outConstraint
				, std::ostream & out
				, const InteractionEnergy & energy
				);

	/**
	 * destruction, which
	 * - prints ensemble information to stream and
	 * - enforces a flush on the output stream.
	 */
	virtual ~OutputHandlerEnsemble();

	/**
	 * Does nothing but counts the call
	 * (no individual interactions are reported)
	 *
	 * @param interaction the interaction to output
	 */
	virtual
	void
	add( const Interaction & interaction );

protected:

	//! counter of reported interactions
	using OutputHandler::reportedInteractions;

	//! aggregated overall partition function
	using OutputHandler::Z;

	//! the output stream to write the ensemble information to
	std::ostream & out;

	//! the interaction energy handler used for the energy computations
	const InteractionEnergy & energy;

};


////////////////////////////////////////////////////////////////////////////

} // namespace

#endif /* INTARNA_OUTPUTHANDLERENSEMBLE_H_ */
