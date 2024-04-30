
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
initZ()
{
	// reset storage
	Z_partition.clear();
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
	if (Z_equal(partZ,0) || Z_isINF(Zall))
		return;
	// handle whether or not partZ includes ED values or not
	Z_type partZ_widthED = 0, partZ_noED = 0;
	if (isHybridZ) {
#if INTARNA_IN_DEBUG_MODE
		if ( (std::numeric_limits<Z_type>::max() - (partZ*energy.getBoltzmannWeight(energy.getE(i1,j1,i2,j2, E_type(0))))) <= Zall) {
			LOG(WARNING) <<"PredictorMfeEns::updateZ() : partition function overflow! Recompile with larger partition function data type!";
		}
#endif
		// add ED penalties etc.
		partZ_noED = partZ;
		partZ_widthED = partZ*energy.getBoltzmannWeight(energy.getE(i1,j1,i2,j2, E_type(0)));
	} else {
#if INTARNA_IN_DEBUG_MODE
		if ( (std::numeric_limits<Z_type>::max() - partZ) <= Zall) {
			LOG(WARNING) <<"PredictorMfeEns::updateZ() : partition function overflow! Recompile with larger partition function data type!";
		}
#endif
		// remove ED
		partZ_noED = partZ / energy.getBoltzmannWeight(energy.getE(i1,j1,i2,j2, E_type(0)));;
		partZ_widthED = partZ;
	}

	// increase overall partition function
	Zall += partZ_withED;

	// store partial Z (without ED)
	Interaction::Boundary key(i1,j1,i2,j2);
	auto keyEntry = Z_partition.find(key);
	if ( Z_partition.find(key) == Z_partition.end() ) {
		Z_partition[key] = partZ_noED;
	} else {
		// update entry
		keyEntry->second += partZ_noED;
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


} // namespace
