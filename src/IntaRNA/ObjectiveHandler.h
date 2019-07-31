#ifndef INTARNA_OBJECTIVEHANDLER_H_
#define INTARNA_OBJECTIVEHANDLER_H_

#include "IntaRNA/Interaction.h"
#include "IntaRNA/InteractionEnergy.h"

namespace IntaRNA
{

/**
 * Defines the optimization objective of Predictor instances
 */
class ObjectiveHandler
{
public:
	ObjectiveHandler();
	virtual ~ObjectiveHandler();


	static
	E_type
	getLcE( const Interaction & i, const InteractionEnergy & energy );

	static
	E_type
	getLcE( const size_t & curLength, const E_type & fullE, const InteractionEnergy & energy );

};

} /* namespace IntaRNA */

#endif /* INTARNA_OBJECTIVEHANDLER_H_ */
