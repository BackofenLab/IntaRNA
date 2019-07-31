

#include "IntaRNA/ObjectiveHandler.h"

namespace IntaRNA
{

ObjectiveHandler::ObjectiveHandler()
{
}

ObjectiveHandler::~ObjectiveHandler()
{
}

E_type
ObjectiveHandler::
getLcE( const Interaction & i, const InteractionEnergy & energy )
{
	const size_t curLength = std::max(
							( 1 + i.basePairs.rbegin()->first - i.basePairs.begin()->first )
							, ( 1 + i.basePairs.begin()->second - i.basePairs.rbegin()->second )
							);

	return getLcE( curLength, i.energy, energy);
}

E_type
ObjectiveHandler::
getLcE( const size_t & curLength, const E_type & fullE, const InteractionEnergy & energy )
{
//	const size_t maxLength = std::max( std::min(energy.getAccessibility1().getSequence().size(), energy.getAccessibility1().getMaxLength())
//							, std::min(energy.getAccessibility2().getSequence().size(), energy.getAccessibility2().getMaxLength()) );
//	return (E_type)( (double)(fullE) / std::log2(1+(double)curLength) );
//	return (E_type)( (double)(fullE) / std::log(5+(double)curLength) );
//	return (E_type)( (double)(fullE) / std::log(2.0*(double)curLength) );
//	return (fullE) / (E_type)curLength;
	return E_type((float)(fullE) / (float)curLength);
}



} /* namespace IntaRNA */
