
#ifndef INTARNA_OUTPUTHANDLERTEXT_H_
#define INTARNA_OUTPUTHANDLERTEXT_H_

#include "IntaRNA/OutputConstraint.h"
#include "IntaRNA/OutputHandler.h"
#include "IntaRNA/InteractionEnergy.h"

#include <iostream>

namespace IntaRNA {

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
	 * @param outConstraint the output constraint applied to find the reported
	 *        interaction
	 * @param out the stream to write to
	 * @param energy the interaction energy object used for computation
	 * @param flankingLength maximal number of nucleotides flanking the
	 *        interaction to be printed in the output
	 * @param detailedOutput if (true) detailed output is provided; normal
	 *        reduced output otherwise
	 */
	OutputHandlerText( const OutputConstraint & outConstraint
				, std::ostream & out
				, const InteractionEnergy & energy
				, const size_t flankingLength = 10
				, const bool detailedOutput = false
				);

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

protected:

	//! counter of reported interactions
	using OutputHandler::reportedInteractions;

	//! aggregated overall partition function
	using OutputHandler::Z;

	//! the output stream to write the interaction text representation to
	std::ostream & out;

	//! the interaction energy handler used for the energy computations
	const InteractionEnergy & energy;

	//! number of flanking bases left/right of the interaction to print
	const size_t flankingLength;

	//! whether or not detailed output has to be provided
	const bool detailedOutput;

};


////////////////////////////////////////////////////////////////////////////

} // namespace

#endif /* OUTPUTHANDLERTEXT_H_ */
