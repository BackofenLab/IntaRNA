
#include "IntaRNA/PredictorMfeEns.h"
#include "IntaRNA/PartitionArithmetic.h"

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
	, updateZisComplete(false)
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

bool
PredictorMfeEns::
addPartitionContribution( const size_t i1, const size_t j1
		, const size_t i2, const size_t j2
		, const Z_type partZ
		, const bool isHybridZ
		, Z_type & partZ_noED )
{
	if (exactContribution) {
		// The coefficient was admitted and computed once by the forward adapter.
		// Reconstructing it here could disagree with the reverse objective.
		if (!isHybridZ || partZ != exactContribution->first)
			throw std::logic_error("exact seeded updateZ override changed a boundary weight; override exactBoundaryWeight instead");
		partZ_noED=PartitionArithmetic::check(partZ);
		const Z_type weighted=PartitionArithmetic::multiply(partZ,exactContribution->second);
		if (weighted==0) return false;
		Zall=PartitionArithmetic::add(Zall,weighted);
		return true;
	}
	// check if something to be done
	if (Z_equal(partZ,0) || Z_isINF(Zall))
		return false;

	// Apply the same site filters used for MFE candidates before changing
	// either the global or boundary-specific partition.
	if (!isValidOutputSite(i1, j1, i2, j2)) {
		return false;
	}

	// handle whether or not partZ includes ED values or not
	Z_type partZ_withED = 0;
	if (isHybridZ) {
#if INTARNA_IN_DEBUG_MODE
		if ( (std::numeric_limits<Z_type>::max() - (partZ*energy.getBoltzmannWeight(energy.getE(i1,j1,i2,j2, E_type(0))))) <= Zall) {
			LOG(WARNING) <<"PredictorMfeEns::updateZ() : partition function overflow! Recompile with larger partition function data type!";
		}
#endif
		// add ED penalties etc.
		partZ_noED = partZ;
		partZ_withED = partZ*energy.getBoltzmannWeight(energy.getE(i1,j1,i2,j2, E_type(0)));
	} else {
#if INTARNA_IN_DEBUG_MODE
		if ( (std::numeric_limits<Z_type>::max() - partZ) <= Zall) {
			LOG(WARNING) <<"PredictorMfeEns::updateZ() : partition function overflow! Recompile with larger partition function data type!";
		}
#endif
		// remove ED
		partZ_noED = partZ / energy.getBoltzmannWeight(energy.getE(i1,j1,i2,j2, E_type(0)));;
		partZ_withED = partZ;
	}

	// increase overall partition function
	Zall += partZ_withED;
	return true;
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns::
updateZ( const size_t i1, const size_t j1
		, const size_t i2, const size_t j2
		, const Z_type partZ
		, const bool isHybridZ )
{
	Z_type partZ_noED = 0;
	if (!addPartitionContribution(i1, j1, i2, j2, partZ, isHybridZ, partZ_noED)) {
		return;
	}

	// Exact 2D prediction finalizes each boundary once. In this scoped mode the
	// partition can update the optima immediately instead of entering the map.
	// Trackers retain the map path to preserve deferred callback ordering.
	if (updateZisComplete && predTracker == NULL) {
		if (Z_isNotINF(partZ_noED) && partZ_noED > 0) {
			PredictorMfe::updateOptima(i1, j1, i2, j2,
					energy.getE(partZ_noED), true, false);
		}
		return;
	}

	// store partial Z (without ED)
	Interaction::Boundary key(i1,j1,i2,j2);
	auto [keyEntry, inserted] = Z_partition.try_emplace(key, partZ_noED);
	if (!inserted) {
		// update entry
		keyEntry->second += partZ_noED;
	}

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns::
updateCompleteZ( const size_t i1, const size_t j1
		, const size_t i2, const size_t j2
		, const Z_type partZ
		, const bool isHybridZ )
{
	class CompleteUpdateScope {
	public:
		explicit CompleteUpdateScope(bool & mode)
		 : mode(mode), previous(mode)
		{
			mode = true;
		}

		~CompleteUpdateScope()
		{
			mode = previous;
		}

	private:
		bool & mode;
		const bool previous;
	};

	CompleteUpdateScope completeUpdate(updateZisComplete);
	updateZ(i1, j1, i2, j2, partZ, isHybridZ);
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns::updateExactCompleteZ(size_t i1,size_t j1,size_t i2,size_t j2,
		Z_type hybrid,Z_type coefficient)
{
	const std::pair<Z_type,Z_type> contribution(hybrid,coefficient);
	struct Restore {
		const std::pair<Z_type,Z_type> * & slot;
		const std::pair<Z_type,Z_type> * previous;
		~Restore() { slot=previous; }
	} restore{exactContribution,exactContribution};
	exactContribution=&contribution;
	updateCompleteZ(i1,j1,i2,j2,hybrid,true);
}

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
