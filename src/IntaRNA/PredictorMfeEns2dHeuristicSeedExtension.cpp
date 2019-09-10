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
	// temporary access
	const OutputConstraint & outConstraint = output.getOutputConstraint();
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
		const size_t maxMatrixLen1 = energy.getAccessibility1().getMaxLength()-sl1+1;
		const size_t maxMatrixLen2 = energy.getAccessibility2().getMaxLength()-sl2+1;
		// check if seed fits into interaction range
		if (sj1 > range_size1 || sj2 > range_size2)
			continue;

		// init optimal right boundaries
		j1opt = sj1;
		j2opt = sj2;
		E_right_opt = E_INF;

		// ER
		hybridZ_right.resize( std::min(range_size1-sj1, maxMatrixLen1), std::min(range_size2-sj2, maxMatrixLen2) );
		fillHybridZ_right(sj1, sj2, si1, si2);

		// EL
		hybridZ_left.resize( std::min(si1+1, maxMatrixLen1), std::min(si2+1, maxMatrixLen2) );
		fillHybridZ_left(si1, si2);

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
	// temporary access
	const OutputConstraint & outConstraint = output.getOutputConstraint();

	// global vars to avoid reallocation
	size_t j1,j2,k1,k2;

	// determine whether or not lonely base pairs are allowed or if we have to
	// ensure a stacking to the right of the left boundary (i1,i2)
	const size_t noLpShift = outConstraint.noLP ? 1 : 0;
	Z_type iStackZ = Z_type(1);

	// current minimal value
	const Z_type seedZ = energy.getBoltzmannWeight(seedHandler.getSeedE(si1, si2));
	const Z_type initZ = energy.getBoltzmannWeight(energy.getE_init());
	// iterate over all window starts j1 (seq1) and j2 (seq2)
	for (j1=sj1; j1-sj1 < hybridZ_right.size1(); j1++ ) {
		for (j2=sj2; j2-sj2 < hybridZ_right.size2(); j2++ ) {

			// referencing cell access
			Z_type & curZ = hybridZ_right(j1-sj1,j2-sj2);

			// init current cell (0 if just left (i1,i2) base pair)
			curZ = sj1==j1 && sj2==j2 ? energy.getBoltzmannWeight(0.0) : 0.0;

			// check if complementary
			// check if not only seed
			if( sj1<j1 && sj2<j2 && energy.areComplementary(j1,j2) ) {

				// left-stacking of j if no-LP
				if (outConstraint.noLP) {
					// skip if no stacking possible
					if (!energy.areComplementary(j1-noLpShift,j2-noLpShift)) {
						continue;
					}
					// get stacking energy to avoid recomputation in recursion below
					iStackZ = energy.getBoltzmannWeight(energy.getE_interLeft(j1-noLpShift,j1,j2-noLpShift,j2));
				}

				// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
				for (k1=j1-noLpShift; k1-- > sj1; ) {
					// ensure maximal loop length
					if (j1-noLpShift-k1 > energy.getMaxInternalLoopSize1()+1) break;
					for (k2=j2-noLpShift; k2-- > sj2; ) {
						// ensure maximal loop length
						if (j2-noLpShift-k2 > energy.getMaxInternalLoopSize2()+1) break;
						// check if (k1,k2) is valid left boundary
						if ( ! Z_equal(hybridZ_right(k1-sj1,k2-sj2), 0.0) ) {
							// store value
							curZ += (hybridZ_right(k1-sj1,k2-sj2)
									* energy.getBoltzmannWeight(energy.getE_interLeft(k1,j1-noLpShift,k2,j2-noLpShift))
									* iStackZ );
						}
					}
				}
				// update overall partition function
				if (Z_isNotINF(curZ) && !Z_equal(curZ,Z_type(0))) {
					// update optimal right extension if needed
					updateOptRightZ( si1,j1,si2,j2, energy.getE(seedZ * curZ * initZ) );
					// update overall partition function information given the current seed
					// seed only not covered due to enclosing check
					if (sj1 != j1) {
						updateZ(si1, j1, si2, j2, seedZ * curZ * initZ, true);
					}
				}
			}

		}
	}

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2dHeuristicSeedExtension::
fillHybridZ_left( const size_t si1, const size_t si2 )
{
	// temporary access
	const OutputConstraint & outConstraint = output.getOutputConstraint();

#if INTARNA_IN_DEBUG_MODE
	// check indices
	if (!energy.areComplementary(si1,si2) )
		throw std::runtime_error("PredictorMfeEns2dSeedExtension::fillHybridZ_left("+toString(si1)+","+toString(si2)+",..) are not complementary");
#endif

	// global vars to avoid reallocation
	size_t i1,i2,k1,k2;
	const Z_type seedZ = energy.getBoltzmannWeight(seedHandler.getSeedE(si1, si2));
	const size_t sj1 = si1 + seedHandler.getSeedLength1(si1,si2) -1;
	const size_t sj2 = si2 + seedHandler.getSeedLength2(si1,si2) -1;
	const Z_type initZ = energy.getBoltzmannWeight(energy.getE_init());
	const Z_type rightOptZ = hybridZ_right(j1opt-sj1,j2opt-sj2);

	// determine whether or not lonely base pairs are allowed or if we have to
	// ensure a stacking to the right of the left boundary (i1,i2)
	const size_t noLpShift = outConstraint.noLP ? 1 : 0;
	Z_type iStackZ = Z_type(0);

	// iterate over all window starts i1 (seq1) and i2 (seq2)
	for (i1=si1; si1-i1 < hybridZ_left.size1(); i1-- ) {
		for (i2=si2; si2-i2 < hybridZ_left.size2(); i2-- ) {

			// referencing cell access
			Z_type & curZ = hybridZ_left(si1-i1,si2-i2);

			// init current cell (0 if not just right-most (j1,j2) base pair)
			curZ = i1==si1 && i2==si2 ? initZ : 0.0;

			// check if complementary (use global sequence indexing)
			if( i1<si1 && i2<si2 && energy.areComplementary(i1,i2) ) {

				// left-stacking of j if no-LP
				if (outConstraint.noLP) {
					// skip if no stacking possible
					if (!energy.areComplementary(i1+noLpShift,i2+noLpShift)) {
						continue;
					}
					// get stacking energy to avoid recomputation in recursion below
					iStackZ = energy.getBoltzmannWeight(energy.getE_interLeft(i1,i1+noLpShift,i2,i2+noLpShift));
				}

				// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
				for (k1=i1+noLpShift; k1++ < si1; ) {
					// ensure maximal loop length
					if (k1-i1-noLpShift > energy.getMaxInternalLoopSize1()+1) break;
					for (k2=i2+noLpShift; k2++ < si2; ) {
						// ensure maximal loop length
						if (k2-i2-noLpShift > energy.getMaxInternalLoopSize2()+1) break;
						// check if (k1,k2) are valid left boundary
						if ( ! Z_equal(hybridZ_left(si1-k1,si2-k2), 0.0) ) {
							curZ += (iStackZ
									* energy.getBoltzmannWeight(energy.getE_interLeft(i1+noLpShift,k1,i2+noLpShift,k2))
									* hybridZ_left(si1-k1,si2-k2) );
						}
					} // k2
				} // k1
			}

			// correction for left seeds
			if (i1 < si1 && i2 < si2 && seedHandler.isSeedBound(i1, i2) ) {

				// check if seed is to be processed:
				bool subtractThisSeed =
						// check if left of anchor seed
									( i1+seedHandler.getSeedLength1(i1,i2)-1 <= si1
									&& i2+seedHandler.getSeedLength2(i1,i2)-1 <= si2 )
						// check if overlapping with anchor seed
								||	seedHandler.areLoopOverlapping(i1,i2,si1,si2);
				if (subtractThisSeed) {

					// iterate seeds in S region
					size_t sj1 = RnaSequence::lastPos, sj2 = RnaSequence::lastPos;
					size_t si1overlap = si1+1, si2overlap = si2+1;
					// find left-most loop-overlapping seed
					while( seedHandler.updateToNextSeed(sj1,sj2
							, i1, std::min(si1,i1+seedHandler.getSeedLength1(i1, i2)-2)
							, i2, std::min(si2,i2+seedHandler.getSeedLength2(i1, i2)-2)) )
					{
						// check if right of i1,i2 and overlapping
						if (sj1 > i1 && seedHandler.areLoopOverlapping(i1, i2, sj1, sj2)) {
							// update left-most loop-overlapping seed
							if (sj1 < si1overlap) {
								si1overlap = sj1;
								si2overlap = sj2;
							}
						}
					}

					// if we found an overlapping seed
					if (si1overlap <= si1) {
						// check if right side is non-empty
						if ( ! E_equal(hybridZ_left( si1-si1overlap, si2-si2overlap),0) ) {
							// compute Energy of loop S \ S'
							E_type nonOverlapE = getNonOverlappingEnergy(i1, i2, si1overlap, si2overlap);
							// subtract energy.getBoltzmannWeight( nonOverlapE ) * hybridZ_left( si1overlap, si2overlap bis anchor seed [==1 falls gleich])
							Z_type correctionTerm = energy.getBoltzmannWeight( nonOverlapE ) * hybridZ_left( si1-si1overlap, si2-si2overlap);
							curZ -= correctionTerm;
							// sanity insurance
							if (curZ < 0) {
								curZ = Z_type(0.0);
							}
						}
					} else {
						// get energy of seed
						E_type seedE = seedHandler.getSeedE(i1, i2);
						// if no S'
						// subtract seedZ * hybridZ_left(right end seed bis anchor seed)
						Z_type correctionTerm = energy.getBoltzmannWeight(seedE) * hybridZ_left( si1-(i1+seedHandler.getSeedLength1(i1, i2)-1), si2-(i2+seedHandler.getSeedLength2(i1, i2)-1));
						curZ -= correctionTerm;
						// sanity insurance
						if (curZ < 0) {
							curZ = Z_type(0.0);
						}
					}
				} // substractThisSeed
			} // curZ correction for left-seeds

			// update overall partition function information given the current seed
			if (Z_isNotINF(curZ)) {

				// Z( left + seed ); covers seed only
				updateZ(i1, sj1, i2, sj2, curZ * seedZ, true);

				// Z( left + seed + rightOpt ) and rightOpt true seed extension
				if (	i1 != si1 // true left seed extension
					&&	j1opt != sj1 // true right seed extension
					// check if interaction width is within boundaries
					&&	j1opt+1-i1 <= energy.getAccessibility1().getMaxLength()
					&& 	j2opt+1-i2 <= energy.getAccessibility2().getMaxLength())
				{
					updateZ(i1, j1opt, i2, j2opt, curZ * seedZ * rightOptZ, true);
				}
			}
		} // i2
	} // i1

}

////////////////////////////////////////////////////////////////////////////


} // namespace
