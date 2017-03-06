
#ifndef OUTPUTHANDLERRANGEONLY_H_
#define OUTPUTHANDLERRANGEONLY_H_

#include "IntaRNA/OutputHandler.h"

namespace IntaRNA {

/**
 * Output handler that transforms a given interaction into an interaction
 * range and forwards all reports to a provided successive output handler.
 *
 * @author Martin Mann
 *
 */
class OutputHandlerRangeOnly : public OutputHandler {

protected:

	//! the output handler to be called for each interaction range
	OutputHandler& successiveOutput;

public:

	/**
	 * construction
	 * @param successiveOutput the output handler to be called for each
	 *   interaction range
	 */
	OutputHandlerRangeOnly( OutputHandler & successiveOutput );

	/**
	 * destruction
	 */
	virtual ~OutputHandlerRangeOnly();


	/**
	 * Converts the given RNA-RNA interaction into an interaction range
	 * and forwards it to the successive output handler
	 *
	 * @param interaction the interaction to transform into an interaction range
	 */
	virtual
	void
	add( const Interaction & interaction );

	/**
	 * Forwards the given RNA-RNA interaction range to the successive output
	 * handler
	 *
	 * @param range the interaction range to add
	 */
	virtual
	void
	add( const InteractionRange & range );


};

} // namespace

#endif /* OUTPUTHANDLERRANGEONLY_H_ */
