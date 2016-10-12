
#ifndef OUTPUTHANDLERTEXT_H_
#define OUTPUTHANDLERTEXT_H_

#include "OutputHandler.h"
#include "InteractionEnergy.h"

#include <iostream>

/**
 * Interaction output handler that writes all interactions directly to stream
 * in a simple text format.
 *
 * @author Martin Mann 2014
 */
class OutputHandlerText: public OutputHandler {
public:

	/**
	 * Construct a simple text output handler for interaction reporting.
	 *
	 * @param out the stream to write to
	 * @param flankingLength maximal number of nucleotides flanking the
	 *        interaction to be printed in the output
	 */
	OutputHandlerText( std::ostream & out
				, const InteractionEnergy & energy
				, const size_t flankingLength = 10 );

	/**
	 * destruction, enforces a flush on the output stream.
	 */
	virtual ~OutputHandlerText();

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

	//! number of flanking bases left/right of the interaction to print
	const size_t flankingLength;

};

#endif /* OUTPUTHANDLERTEXT_H_ */
