
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
class OutputHandlerIntaRNA1: public OutputHandler {
public:

	/**
	 * Construct a simple text output handler for interaction reporting
	 * similar to the detailed output of IntaRNA version 1.2.*
	 *
	 * @param out the stream to write to
	 * @param energy the interaction energy object used for computation
	 * @param detailedOutput whether or not to produce detailed (true) or
	 *        normal (false) IntaRNA v1* output
	 */
	OutputHandlerIntaRNA1( std::ostream & out
				, const InteractionEnergy & energy
				, const bool detailedOutput );

	/**
	 * destruction, enforces a flush on the output stream.
	 */
	virtual ~OutputHandlerIntaRNA1();

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

	/**
	 * Defines whether or not a leading separator is to be printed before any
	 * output for this output handler.
	 * @param yesNo whether or not a separator is to be printed
	 */
	virtual
	void
	addSeparator (const bool yesNo );

protected:

	//! whether or not to produce detailed (true) or normal (false) output
	const bool detailedOutput;

	//! the output stream to write the interaction text representation to
	std::ostream & out;

	//! the interaction energy handler used for the energy computations
	const InteractionEnergy & energy;

	//! whether or not the initial output was done
	bool initialOutputDone;

	//! whether or not to print a leading separator with the initial output
	bool printSeparator;


};

#endif /* OUTPUTHANDLERINTARNA1DETAILED_H_ */
