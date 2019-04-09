
#include "IntaRNA/PredictorMfeEns.h"

#include <iostream>
#include <algorithm>

namespace IntaRNA {

////////////////////////////////////////////////////////////////////////////

PredictorMfeEns::PredictorMfeEns(
		const InteractionEnergy & energy
		, OutputHandler & output
		, PredictionTracker * predTracker
		)
	: PredictorMfe(energy,output,predTracker)
	, overallZhybrid(0)
{

}

////////////////////////////////////////////////////////////////////////////

PredictorMfeEns::~PredictorMfeEns()
{
}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns::
initZ( const OutputConstraint & outConstraint )
{
	// reinit overall partition function
	overallZhybrid = Z_type(0);

	// reinit mfe information
	initOptima( outConstraint );
}

////////////////////////////////////////////////////////////////////////////

Z_type
PredictorMfeEns::
getHybridZ() const
{
	return overallZhybrid;
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns::
updateZ( const size_t i1, const size_t j1
		, const size_t i2, const size_t j2
		, const Z_type partZ
		, const bool isHybridZ )
{
	// update overall hybridization partition function
	if (isHybridZ) {
		overallZhybrid += partZ;
	} else {
		// remove ED, dangling end contributions, etc. before adding
		overallZhybrid += ( partZ / energy.getBoltzmannWeight(energy.getE(i1,j1,i2,j2, E_type(0))) );
	}
}

////////////////////////////////////////////////////////////////////////////


} // namespace
