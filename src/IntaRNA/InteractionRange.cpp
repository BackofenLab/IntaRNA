
#include "IntaRNA/Interaction.h"
#include "IntaRNA/InteractionRange.h"

namespace IntaRNA {

InteractionRange &
InteractionRange::
operator= ( const Interaction & interaction )
{
#if INTARNA_IN_DEBUG_MODE
	if (interaction.isEmpty())
		throw std::runtime_error("InteractionRange::=(interaction) is empty!");
	if (!interaction.isValid())
		throw std::runtime_error("InteractionRange::=(interaction) not valid!");
#endif
	// init data
	s1 = interaction.s1;
	s2 = interaction.s2;
	r1.from = interaction.basePairs.begin()->first;
	r1.to = interaction.basePairs.rbegin()->first;
	r2.from = interaction.basePairs.begin()->second;
	r2.to = interaction.basePairs.rbegin()->second;
	energy = interaction.energy;

#if INTARNA_IN_DEBUG_MODE
	if (!isSane())
		throw std::runtime_error("InteractionRange::=(interaction)="+toString(*this)+" not sane!");
#endif
	return *this;
}




} // namespace
