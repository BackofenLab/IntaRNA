
#include "IntaRNA/PredictorMfeEnsSeedOnly.h"

namespace IntaRNA {

//////////////////////////////////////////////////////////////////////////

PredictorMfeEnsSeedOnly::
PredictorMfeEnsSeedOnly(
		const InteractionEnergy & energy
		, OutputHandler & output
		, PredictionTracker * predTracker
		, SeedHandler * seedHandlerInstance )
 :
	PredictorMfeEns(energy,output,predTracker)
	, seedHandler(seedHandlerInstance)
	, seedLastPos1(0)
	, seedLastPos2(0)
{
}

//////////////////////////////////////////////////////////////////////////

PredictorMfeEnsSeedOnly::
~PredictorMfeEnsSeedOnly()
{
}

//////////////////////////////////////////////////////////////////////////

void
PredictorMfeEnsSeedOnly::
predict( const IndexRange & r1, const IndexRange & r2 )
{
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) <<"predicting seed interactions only with partition function ..."; }
	// measure timing
	TIMED_FUNC_IF(timerObj,VLOG_IS_ON(9));


#if INTARNA_IN_DEBUG_MODE
	// check indices
	if (!(r1.isAscending() && r2.isAscending()) )
		throw std::runtime_error("PredictorMfeEnsSeedOnly::predict("+toString(r1)+","+toString(r2)+") is not sane");
#endif

	// setup index offset
	energy.setOffset1(r1.from);
	energy.setOffset2(r2.from);
	seedHandler.setOffset1(r1.from);
	seedHandler.setOffset2(r2.from);

	seedLastPos1 = std::min( energy.size1()
			, (r1.to==RnaSequence::lastPos?energy.size1()-1:r1.to)-r1.from+1 ) -1;
	seedLastPos2 = std::min( energy.size2()
			, (r2.to==RnaSequence::lastPos?energy.size2()-1:r2.to)-r2.from+1 ) -1;

	// trigger empty interaction reporting
	initOptima();
	initZ();

	// compute seed interactions for whole range
	// and check if any seed possible
	if (seedHandler.fillSeed( 0, seedLastPos1, 0, seedLastPos2 ) == 0) {
		reportOptima();
		// stop computation
		return;
	}

	// iterate over all seeds within range
	size_t i1 = RnaSequence::lastPos, i2 = RnaSequence::lastPos;
	size_t j1 = i1, j2 = i2;
	while( seedHandler.updateToNextSeed( i1, i2, 0, seedLastPos1, 0, seedLastPos2) ) {
		// check of seed end is within range
		j1 = i1 + seedHandler.getSeedLength1(i1,i2) -1;
		if (j1 > seedLastPos1) continue;
		j2 = i2 + seedHandler.getSeedLength2(i1,i2) -1;
		if (j2 > seedLastPos2) continue;
		const E_type seedEhybrid = seedHandler.getSeedE(i1,i2) + energy.getE_init();
		// report seed energy (including initialization term)
		updateOptima( i1, j1, i2, j2, seedEhybrid, true, true );
	}

	// report mfe interaction
	reportOptima();

}


//////////////////////////////////////////////////////////////////////////

void
PredictorMfeEnsSeedOnly::
traceBack( Interaction & interaction )
{
	// check if something to trace
	if (interaction.basePairs.size() < 2) {
		return;
	}

#if INTARNA_IN_DEBUG_MODE
	// sanity checks
	if ( interaction.basePairs.size() != 2 ) {
		throw std::runtime_error("PredictorMfeEnsSeedOnly::traceBack() : given interaction does not contain boundaries only");
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
		throw std::runtime_error("PredictorMfeEnsSeedOnly::traceBack() : given interaction not valid");
	}
#endif

	// get indices in hybridE for boundary base pairs
	size_t	i1 = energy.getIndex1(interaction.basePairs.at(0)),
			j1 = energy.getIndex1(interaction.basePairs.at(1)),
			i2 = energy.getIndex2(interaction.basePairs.at(0)),
			j2 = energy.getIndex2(interaction.basePairs.at(1));

#if INTARNA_IN_DEBUG_MODE
	// check if seed interval
	if (!seedHandler.isSeedBound(i1,i2)
		|| (j1 != i1+seedHandler.getSeedLength1(i1,i2)-1)
		|| (j2 != i2+seedHandler.getSeedLength2(i1,i2)-1) )
	{
		// no seed possible, abort computation
		throw std::runtime_error("PredictorMfeEnsSeedOnly::traceBack() : given boundaries "+toString(interaction)+" are not only containing a seed");
	}
#endif


	// store seed information
	interaction.setSeedRange(
			interaction.basePairs.at(0),
			interaction.basePairs.at(1),
			energy.getE(i1,j1,i2,j2, seedHandler.getSeedE(i1,i2)+energy.getE_init()));
	// trace back seed base pairs
	seedHandler.traceBackSeed( interaction, i1, i2 );

	// sort final interaction (to make valid) (faster than calling sort())
	if (interaction.basePairs.size() > 2) {
		Interaction::PairingVec & bps = interaction.basePairs;
		// shift all added base pairs to the front
		for (size_t i=2; i<bps.size(); i++) {
			bps.at(i-1).first = bps.at(i).first;
			bps.at(i-1).second = bps.at(i).second;
		}
		// set last to j1-j2
		(*bps.rbegin()) = energy.getBasePair( j1, j2 );
	}

}

//////////////////////////////////////////////////////////////////////////




} // namespace
