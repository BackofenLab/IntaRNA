#include "IntaRNA/PredictorMfeEns2dHeuristicSeedExtension.h"

namespace IntaRNA {

//////////////////////////////////////////////////////////////////////////

PredictorMfeEns2dHeuristicSeedExtension::
PredictorMfeEns2dHeuristicSeedExtension(
		const InteractionEnergy & energy
		, OutputHandler & output
		, PredictionTracker * predTracker
		, SeedHandler * seedHandlerInstance )
 :
	PredictorMfeEns2dSeedExtension(energy,output,predTracker,seedHandlerInstance)
, E_right_opt(E_INF)
, j1opt(0)
, j2opt(0)
{
	assert( seedHandler.getConstraint().getBasePairs() > 1 );
}

//////////////////////////////////////////////////////////////////////////

PredictorMfeEns2dHeuristicSeedExtension::
~PredictorMfeEns2dHeuristicSeedExtension()
{
}

//////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2dHeuristicSeedExtension::
predict( const IndexRange & r1, const IndexRange & r2 )
{
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) <<"predicting ensemble mfe interactions with seed-extension heuristically in O(n^2) space and O(n^2) time..."; }
	// measure timing
	TIMED_FUNC_IF(timerObj,VLOG_IS_ON(9));

#if INTARNA_IN_DEBUG_MODE
	// check indices
	if (!(r1.isAscending() && r2.isAscending()) )
		throw std::runtime_error("PredictorMfe2dHeuristicSeedExtension::predict("+toString(r1)+","+toString(r2)+") is not sane");
#endif

	// setup index offset
	energy.setOffset1(r1.from);
	energy.setOffset2(r2.from);
	seedHandler.setOffset1(r1.from);
	seedHandler.setOffset2(r2.from);

	const size_t range_size1 = std::min( energy.size1()
			, (r1.to==RnaSequence::lastPos?energy.size1()-1:r1.to)-r1.from+1 );
	const size_t range_size2 = std::min( energy.size2()
			, (r2.to==RnaSequence::lastPos?energy.size2()-1:r2.to)-r2.from+1 );

	// compute seed interactions for whole range
	// and check if any seed possible
	if (seedHandler.fillSeed( 0, range_size1-1, 0, range_size2-1 ) == 0) {
		// trigger empty interaction reporting
		initOptima();
		reportOptima();
		// stop computation
		return;
	}

	// initialize mfe interaction for updates
	initOptima();
	// initialize overall partition function for updates
	initZ();

	size_t si1 = RnaSequence::lastPos, si2 = RnaSequence::lastPos;
	while( seedHandler.updateToNextSeed(si1,si2
			, 0, range_size1+1-seedHandler.getConstraint().getBasePairs()
			, 0, range_size2+1-seedHandler.getConstraint().getBasePairs()) )
	{
		const Z_type seedZ = energy.getBoltzmannWeight( seedHandler.getSeedE(si1, si2) );
		const size_t sl1 = seedHandler.getSeedLength1(si1, si2);
		const size_t sl2 = seedHandler.getSeedLength2(si1, si2);
		const size_t sj1 = si1+sl1-1;
		const size_t sj2 = si2+sl2-1;
		// check if seed fits into interaction range
		if (sj1 > range_size1 || sj2 > range_size2)
			continue;
		const size_t maxMatrixLen1 = energy.getAccessibility1().getMaxLength()-sl1+1;
		const size_t maxMatrixLen2 = energy.getAccessibility2().getMaxLength()-sl2+1;

		// init optimal right boundaries
		j1opt = sj1;
		j2opt = sj2;
		E_right_opt = E_INF;

		// ER
		hybridZ_right.resize( std::min(range_size1-sj1, maxMatrixLen1), std::min(range_size2-sj2, maxMatrixLen2) );
		fillHybridZ_right(sj1, sj2, si1, si2);

		// ensure there is a valid right-extension
		if (!output.getOutputConstraint().noGUend || (!E_isINF(E_right_opt) || !energy.isGU(sj1,sj2))) {
			// EL
			hybridZ_left.resize( std::min(si1+1, maxMatrixLen1), std::min(si2+1, maxMatrixLen2) );
			fillHybridZ_left(si1, si2);
		}

	} // si1 / si2

	// report mfe interaction
	reportOptima();

}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2dHeuristicSeedExtension::
fillHybridZ_right( const size_t sj1, const size_t sj2
			, const size_t si1, const size_t si2 )
{

	// static data
	const Z_type seedZ = energy.getBoltzmannWeight(seedHandler.getSeedE(si1, si2));
	const Z_type initZ = energy.getBoltzmannWeight(energy.getE_init());

	// compute right-extensions of current seed
	PredictorMfeEns2dSeedExtension::fillHybridZ_right(sj1,sj2);

	// update partition function information
	for (size_t r1=1; r1 < hybridZ_right.size1(); r1++ ) {
		for (size_t r2=1; r2 < hybridZ_right.size2(); r2++ ) {

			// referencing cell access
			Z_type & rightExtZ = hybridZ_right(r1,r2);

			// update overall partition function
			if (!Z_equal(rightExtZ,Z_type(0))) {
				// update optimal right extension if needed
				updateOptRightZ( si1,sj1+r1,si2,sj2+r2, energy.getE(seedZ * rightExtZ * initZ) );
				// update overall partition function information for true right-extensions of the current seed
				// seed only not covered due to min-val of r1,r2
				updateZ(si1, sj1+r1, si2, sj2+r2, seedZ * rightExtZ * initZ, true);
			}

		}
	}

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2dHeuristicSeedExtension::
fillHybridZ_left( const size_t si1, const size_t si2 )
{
	// static data
	const Z_type seedZ = energy.getBoltzmannWeight(seedHandler.getSeedE(si1, si2));
	const size_t sj1 = si1 + seedHandler.getSeedLength1(si1,si2) -1;
	const size_t sj2 = si2 + seedHandler.getSeedLength2(si1,si2) -1;
	const Z_type rightOptZ = hybridZ_right(j1opt-sj1,j2opt-sj2);

	// compute left-extensions of current seed
	PredictorMfeEns2dSeedExtension::fillHybridZ_left(si1,si2);

	// update partition function information
	for (size_t l1=0; l1 < hybridZ_left.size1(); l1++ ) {
		for (size_t l2=0; l2 < hybridZ_left.size2(); l2++ ) {

			// referencing cell access
			const Z_type & curZ = hybridZ_left(l1,l2);

			// update overall partition function information given the current seed
			if ( !Z_equal(curZ,0.0) ) {

				// Z( left + seed ); covers seed only
				updateZ(si1-l1, sj1, si2-l2, sj2, curZ * seedZ, true);

				// Z( left + seed + rightOpt ) and rightOpt true seed extension
				if (	l1 > 0 // true left seed extension
					&&	j1opt != sj1 // true right seed extension
					&&	!Z_equal(rightOptZ,0.0) // true right seed extension
					// check if interaction width is within boundaries
					&&	j1opt+1-(si1-l1) <= energy.getAccessibility1().getMaxLength()
					&& 	j2opt+1-(si2-l2) <= energy.getAccessibility2().getMaxLength()
					)
				{
					updateZ(si1-l1, j1opt, si2-l2, j2opt, curZ * seedZ * rightOptZ, true);
				}
			}
		} // i2
	} // i1

}

////////////////////////////////////////////////////////////////////////////


} // namespace
