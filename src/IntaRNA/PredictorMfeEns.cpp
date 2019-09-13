
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

Z_type
PredictorMfeEns::
getHybridZ( const size_t i1, const size_t j1
	 , const size_t i2, const size_t j2)
{
	size_t key = generateMapKey(i1, j1, i2, j2);
	if ( Z_partitions.find(key) == Z_partitions.end() ) {
		return Z_type(0);
	} else {
		return Z_partitions[key].partZ;
	}
}

////////////////////////////////////////////////////////////////////////////

Z_type
PredictorMfeEns::
getZ( const size_t i1, const size_t j1
	 , const size_t i2, const size_t j2)
{
	size_t key = generateMapKey(i1, j1, i2, j2);
	if ( Z_partitions.find(key) == Z_partitions.end() ) {
		return Z_type(0);
	} else {
		return Z_partitions[key].partZ * energy.getBoltzmannWeight(energy.getE(i1,j1,i2,j2, E_type(0)));
	}
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns::
reportZ( SeedHandler* seedHandler )
{
	if (predTracker != NULL) {
		predTracker->updateZ(this, seedHandler);
	}
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns::
checkKeyBoundaries( const size_t maxLength )
{
	// reset storage
	Z_partition.clear();
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns::
updateZ( const size_t i1, const size_t j1
		, const size_t i2, const size_t j2
		, const Z_type partZ )
{
	// check if something to be done
	if (Z_equal(partZ,0) || Z_isINF(Zall))
		return;
	// update overall partition function
#if INTARNA_IN_DEBUG_MODE
		if ( (std::numeric_limits<Z_type>::max() - (partZ*energy.getBoltzmannWeight(energy.getE(i1,j1,i2,j2, E_type(0))))) <= Zall) {
			LOG(WARNING) <<"PredictorMfeEns::updateZ() : partition function overflow! Recompile with larger partition function data type!";
		}
#endif
		// add ED penalties etc.
		Zall += partZ*energy.getBoltzmannWeight(energy.getE(i1,j1,i2,j2, E_type(0)));
	} else {
#if INTARNA_IN_DEBUG_MODE
		if ( (std::numeric_limits<Z_type>::max() - partZ) <= Zall) {
			LOG(WARNING) <<"PredictorMfeEns::updateZ() : partition function overflow! Recompile with larger partition function data type!";
		}
#endif
		// just increase
		Zall += partZ;
	}
#endif
	// add ED penalties etc.
	overallZ += partZ * energy.getBoltzmannWeight(energy.getE(i1,j1,i2,j2, E_type(0)));

	// store partial Z
	size_t key = generateMapKey(i1, j1, i2, j2);
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
		keyEntry->second += partZ;
	}

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns::
updateOptimaUsingZ()
{
	for (auto it = Z_partition.begin(); it != Z_partition.end(); ++it)
	{
		// if partition function is > 0
		if (Z_isNotINF(it->second) && it->second > 0) {
			PredictorMfe::updateOptima( it->first.i1, it->first.j1, it->first.i2, it->first.j2, energy.getE(it->second), true, false );
		}
	}
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns::
reportOptima()
{
	// update optima from Z information
	updateOptimaUsingZ();

	// call super-class function
	PredictorMfe::reportOptima();
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns::
traceBack( Interaction & interaction )
{
	// temporary access
	const OutputConstraint & outConstraint = output.getOutputConstraint();
	// check if something to trace
	if (interaction.basePairs.size() < 2) {
		return;
	}

#if INTARNA_IN_DEBUG_MODE
	// sanity checks
	if ( interaction.basePairs.size() != 2 ) {
		throw std::runtime_error("PredictorMfeEns::traceBack() : given interaction does not contain boundaries only");
	}
#endif

	// ensure sorting
	interaction.sort();

	// check for single base pair interaction
	if (interaction.basePairs.at(0).first == interaction.basePairs.at(1).first) {
		// delete second boundary (identical to first)
		interaction.basePairs.resize(1);
		// update done
		return;
	}

#if INTARNA_IN_DEBUG_MODE
	// sanity checks
	if ( ! interaction.isValid() ) {
		throw std::runtime_error("PredictorMfeEns2d::traceBack() : given interaction is not valid");
	}
#endif

}

////////////////////////////////////////////////////////////////////////////

size_t
PredictorMfeEns::
generateMapKey( const size_t i1, const size_t j1
					, const size_t i2, const size_t j2 ) const
{
	size_t maxLength = std::max(energy.getAccessibility1().getMaxLength(), energy.getAccessibility2().getMaxLength());
	size_t key = 0;
	key += i1;
	key += j1 * pow(maxLength, 1);
	key += i2 * pow(maxLength, 2);
	key += j2 * pow(maxLength, 3);
	return key;
}

////////////////////////////////////////////////////////////////////////////


} // namespace
