
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
getOverallZ() const
{
	return overallZ;
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns::
reportZ() const
{
	if (predTracker != NULL) {
		predTracker->updateZ(this);
	}
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns::
checkKeyBoundaries( const size_t maxLength )
{
	// check if getMaxLength > sqrt3(size_t) -> error
	if (maxLength > cbrt(std::numeric_limits<size_t>::max())) {
		throw std::runtime_error("PredictorMfeEns::checkKeyBoundaries() : maxLength too big for key generation (out of bounds)");
	}
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns::
updateZ( const size_t i1, const size_t j1
		, const size_t i2, const size_t j2
		, const Z_type partZ
		, const bool isHybridZ )
{
	// check if something to be done
	if (Z_equal(partZ,0) || Z_isINF(overallZ))
		return;
	// update overall partition function
	if (isHybridZ) {
#if INTARNA_IN_DEBUG_MODE
		if ( (std::numeric_limits<Z_type>::max() - (partZ*energy.getBoltzmannWeight(energy.getE(i1,j1,i2,j2, E_type(0))))) <= overallZ) {
			LOG(WARNING) <<"PredictorMfeEns::updateZ() : partition function overflow! Recompile with larger partition function data type!";
		}
#endif
		// add ED penalties etc.
		overallZ += partZ*energy.getBoltzmannWeight(energy.getE(i1,j1,i2,j2, E_type(0)));
	} else {
#if INTARNA_IN_DEBUG_MODE
		if ( (std::numeric_limits<Z_type>::max() - partZ) <= overallZ) {
			LOG(WARNING) <<"PredictorMfeEns::updateZ() : partition function overflow! Recompile with larger partition function data type!";
		}
#endif
		// just increase
		overallZ += partZ;
	}

// TODO : was soll das hier?
	// store partial Z
	size_t maxLength = std::max(energy.getAccessibility1().getMaxLength(), energy.getAccessibility2().getMaxLength());
	size_t key = 0;
	key += i1;
	key += j1 * pow(maxLength, 1);
	key += i2 * pow(maxLength, 2);
	key += j2 * pow(maxLength, 3);
	if ( Z_partitions.find(key) == Z_partitions.end() ) {
		// create new entry
		ZPartition zPartition;
		zPartition.i1 = i1;
		zPartition.j1 = j1;
		zPartition.i2 = i2;
		zPartition.j2 = j2;
		zPartition.partZ = partZ;
		Z_partitions[key] = zPartition;
	} else {
		// update entry
		ZPartition & zPartition = Z_partitions[key];
		zPartition.partZ += partZ;
	}
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns::
reportOptima( const OutputConstraint & outConstraint )
{
	// store overall partition function
	output.incrementZ( getOverallZ() );
	// report optima
	PredictorMfe::reportOptima( outConstraint );
}

////////////////////////////////////////////////////////////////////////////


} // namespace
