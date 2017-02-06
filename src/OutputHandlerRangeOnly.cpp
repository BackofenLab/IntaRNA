
#include "OutputHandlerRangeOnly.h"

OutputHandlerRangeOnly::OutputHandlerRangeOnly( OutputHandler & successiveOutput )
	:
	successiveOutput(successiveOutput)
{
}

OutputHandlerRangeOnly::~OutputHandlerRangeOnly() {
}

void
OutputHandlerRangeOnly::
add( const Interaction & interaction )
{
	// convert and forward
	successiveOutput.add( InteractionRange(interaction) );
}

void
OutputHandlerRangeOnly::
add( const InteractionRange & range )
{
	// forward
	successiveOutput.add( range );
}
