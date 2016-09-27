
#include "InteractionEnergy.h"


////////////////////////////////////////////////////////////////////////////

InteractionEnergy::InteractionEnergy( const Accessibility & accS1
				, const ReverseAccessibility & accS2
				, const size_t maxInternalLoopSize1
				, const size_t maxInternalLoopSize2
		)
  :
	accS1(accS1)
	, accS2(accS2)
	, maxInternalLoopSize1(maxInternalLoopSize1)
	, maxInternalLoopSize2(maxInternalLoopSize2)

{
}

////////////////////////////////////////////////////////////////////////////

InteractionEnergy::~InteractionEnergy()
{
}

////////////////////////////////////////////////////////////////////////////

bool
InteractionEnergy::
isAllowedLoopRegion( const RnaSequence& seq, const size_t i, const size_t j, const size_t maxInternalLoopSize )
{
	// ensure index and loop size validity
	return	   i < seq.size()
			&& j < seq.size()
			&& seq.asString().at(i) != 'N'
			&& seq.asString().at(j) != 'N'
			&& i <= j
			&& (j-i) < maxInternalLoopSize;

}

////////////////////////////////////////////////////////////////////////////

size_t
InteractionEnergy::
getLength1() const
{
	// sequence length access
	return accS1.getSequence().size();
}

////////////////////////////////////////////////////////////////////////////

size_t
InteractionEnergy::
getLength2() const
{
	// sequence length access
	return accS2.getSequence().size();
}

////////////////////////////////////////////////////////////////////////////

