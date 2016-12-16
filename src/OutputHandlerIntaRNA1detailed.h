
#ifndef OUTPUTHANDLERINTARNA1DETAILED_H_
#define OUTPUTHANDLERINTARNA1DETAILED_H_

#include "OutputHandler.h"
#include "InteractionEnergy.h"

#include <iostream>

/**
 * Interaction output handler that writes all interactions directly to stream
 * as done by IntaRNA version 1.2.* in detailed output mode.
 *
 * @author Martin Mann 2016
 */
class OutputHandlerIntaRNA1detailed: public OutputHandler {
public:

	/**
	 * Construct a simple text output handler for interaction reporting
	 * similar to the detailed output of IntaRNA version 1.2.*
	 *
	 * @param out the stream to write to
	 * @param energy the interaction energy object used for computation
	 */
	OutputHandlerIntaRNA1detailed( std::ostream & out
				, const InteractionEnergy & energy );

	/**
	 * destruction, enforces a flush on the output stream.
	 */
	virtual ~OutputHandlerIntaRNA1detailed();

	/**
	 * Write a given RNA-RNA interaction in simple text format to the output
	 * stream.
	 *
	 * @param interaction the interaction to output
	 */
	virtual
	void
	add( const Interaction & interaction );

	/**
	 * Handles a given RNA-RNA interaction range as a
	 * RNA-RNA interaction with two base pairs and writes it in simple
	 * text format to the output stream.
	 *
	 * @param range the interaction range to add
	 */
	virtual
	void
	add( const InteractionRange & range );

protected:

	//! the output stream to write the interaction text representation to
	std::ostream & out;

	//! the interaction energy handler used for the energy computations
	const InteractionEnergy & energy;


};

#endif /* OUTPUTHANDLERINTARNA1DETAILED_H_ */
