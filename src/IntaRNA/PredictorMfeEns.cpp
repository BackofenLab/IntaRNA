
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
	, overallZ(0)
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
	overallZ = Z_type(0);
}

////////////////////////////////////////////////////////////////////////////

Z_type
PredictorMfeEns::
getHybridZ() const
{
	return overallZ;
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
#if INTARNA_IN_DEBUG_MODE
		if ( (std::numeric_limits<Z_type>::max() - (partZ*energy.getBoltzmannWeight(energy.getE(i1,j1,i2,j2, E_type(0))))) <= overallZ) {
			LOG(WARNING) <<"PredictorMfeEns::updateZ() : partition function overflow! Recompile with larger partition function data type!";
		}
#endif
		overallZ += partZ*energy.getBoltzmannWeight(energy.getE(i1,j1,i2,j2, E_type(0)));
	} else {
#if INTARNA_IN_DEBUG_MODE
		if ( (std::numeric_limits<Z_type>::max() - partZ) <= overallZ) {
			LOG(WARNING) <<"PredictorMfeEns::updateZ() : partition function overflow! Recompile with larger partition function data type!";
		}
#endif
		// remove ED, dangling end contributions, etc. before adding
		overallZ += partZ;
	}
}

////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns::
reportOptima( const OutputConstraint & outConstraint )
{
	// store overall partition function
	output.incrementZ( getHybridZ() );
	// report optima
	PredictorMfe::reportOptima( outConstraint );
}

////////////////////////////////////////////////////////////////////////////


} // namespace
