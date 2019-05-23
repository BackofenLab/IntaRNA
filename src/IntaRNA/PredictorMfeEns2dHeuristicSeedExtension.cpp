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

	checkKeyBoundaries(std::max(energy.getAccessibility1().getMaxLength(), energy.getAccessibility2().getMaxLength()));
}

//////////////////////////////////////////////////////////////////////////

PredictorMfeEns2dHeuristicSeedExtension::
~PredictorMfeEns2dHeuristicSeedExtension()
{
}

//////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2dHeuristicSeedExtension::
predict( const IndexRange & r1, const IndexRange & r2
		, const OutputConstraint & outConstraint )
{
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) <<"predicting mfe interactions with seed in O(n^2) space and O(n^4) time..."; }
	// measure timing
	TIMED_FUNC_IF(timerObj,VLOG_IS_ON(9));

	// suboptimal setup check
	if (outConstraint.reportMax>1 && outConstraint.reportOverlap != OutputConstraint::ReportOverlap::OVERLAP_BOTH) {
		throw std::runtime_error("PredictorMfeEns2dHeuristicSeedExtension : the enumeration of non-overlapping suboptimal interactions is not supported in this prediction mode");
	}

#if INTARNA_IN_DEBUG_MODE
	// check indices
	if (!(r1.isAscending() && r2.isAscending()) )
		throw std::runtime_error("PredictorMfe2d::predict("+toString(r1)+","+toString(r2)+") is not sane");
#endif

	// setup index offset
	energy.setOffset1(r1.from);
	energy.setOffset2(r2.from);
	seedHandler.setOffset1(r1.from);
	seedHandler.setOffset2(r2.from);

	const size_t interaction_size1 = std::min( energy.size1()
			, (r1.to==RnaSequence::lastPos?energy.size1()-1:r1.to)-r1.from+1 );
	const size_t interaction_size2 = std::min( energy.size2()
			, (r2.to==RnaSequence::lastPos?energy.size2()-1:r2.to)-r2.from+1 );

	// compute seed interactions for whole range
	// and check if any seed possible
	if (seedHandler.fillSeed( 0, interaction_size1-1, 0, interaction_size2-1 ) == 0) {
		// trigger empty interaction reporting
		initOptima(outConstraint);
		reportOptima(outConstraint);
		// stop computation
		return;
	}

	// initialize mfe interaction for updates
	initOptima( outConstraint );
	// initialize overall partition function for updates
	initZ( outConstraint );

	size_t si1 = RnaSequence::lastPos, si2 = RnaSequence::lastPos;
	while( seedHandler.updateToNextSeed(si1,si2
			, 0, interaction_size1+1-seedHandler.getConstraint().getBasePairs()
			, 0, interaction_size2+1-seedHandler.getConstraint().getBasePairs()) )
	{
		const Z_type seedZ = energy.getBoltzmannWeight( seedHandler.getSeedE(si1, si2) );
		const size_t sl1 = seedHandler.getSeedLength1(si1, si2);
		const size_t sl2 = seedHandler.getSeedLength2(si1, si2);
		const size_t sj1 = si1+sl1-1;
		const size_t sj2 = si2+sl2-1;
		const size_t maxMatrixLen1 = energy.getAccessibility1().getMaxLength()-sl1+1;
		const size_t maxMatrixLen2 = energy.getAccessibility2().getMaxLength()-sl2+1;
		// check if seed fits into interaction range
		if (sj1 > interaction_size1 || sj2 > interaction_size2)
			continue;

		// init optimal right boundaries
		j1opt = sj1;
		j2opt = sj2;
		E_right_opt = E_INF;

		// update mfe for seed only
		PredictorMfeEns::updateOptima( si1,sj1,si2,sj2, energy.getE_init() + seedHandler.getSeedE(si1, si2), true );

		// ER
		hybridZ_right.resize( std::min(interaction_size1-sj1, maxMatrixLen1), std::min(interaction_size2-sj2, maxMatrixLen2) );
		fillHybridZ_right(sj1, sj2, outConstraint, si1, si2);

		// EL for fixed right boundary jopt
		hybridZ_left.resize( std::min(si1+1, maxMatrixLen1-(j1opt-sj1)), std::min(si2+1, maxMatrixLen2-(j2opt-sj2)) );
		fillHybridZ_left(si1, si2, outConstraint);

		// updateZ for all boundary combinations
		for (size_t i1 = 0; i1<hybridZ_left.size1(); i1++) {
			// ensure max interaction length in seq 1
			for (size_t j1 = 0; j1 < j1opt-sj1+1 ; j1++) {
				// check interaction length
				if (sj1+j1-si1+i1 >= energy.getAccessibility1().getMaxLength()) continue;
				for (size_t i2 = 0; i2 < hybridZ_left.size2(); i2++) {
					// check complementarity of boundary
					if ( Z_equal(hybridZ_left(i1,i2), 0.0) ) continue;
					// ensure max interaction length in seq 2
					for (size_t j2 = 0; j2 < j2opt-sj2+1 ; j2++) {
						// check interaction length
						if (sj2+j2-si2+i2 >= energy.getAccessibility2().getMaxLength()) continue;
						// check complementarity of boundary
						if (Z_equal(hybridZ_right(j1,j2),0.0)) continue;
						// compute partition function given the current seed
						updateZ(si1-i1, sj1+j1, si2-i2, sj2+j2, hybridZ_left(i1,i2) * seedZ * hybridZ_right(j1,j2), true);
					} // dj2
				} // di2
			} // dj1
		} // di1

	} // si1 / si2

	// report mfe interaction
	reportOptima( outConstraint );

	LOG(DEBUG) <<"Overall Z = "<<getOverallZ();
	LOG(DEBUG) <<"Overall E = "<<E_2_Ekcal(energy.getE(getOverallZ()));

}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2dHeuristicSeedExtension::
fillHybridZ_right( const size_t i1, const size_t i2
			, const OutputConstraint & outConstraint
			, const size_t si1, const size_t si2 )
{

	// global vars to avoid reallocation
	size_t j1,j2,k1,k2;

	// current minimal value
	const Z_type seedZ = energy.getBoltzmannWeight(seedHandler.getSeedE(si1, si2));
	const Z_type initZ = energy.getBoltzmannWeight(energy.getE_init());
	// iterate over all window starts j1 (seq1) and j2 (seq2)
	for (j1=i1; j1-i1 < hybridZ_right.size1(); j1++ ) {
		for (j2=i2; j2-i2 < hybridZ_right.size2(); j2++ ) {

			// init current cell (0 if just left (i1,i2) base pair)
			hybridZ_right(j1-i1,j2-i2) = i1==j1 && i2==j2 ? energy.getBoltzmannWeight(0.0) : 0.0;

			// check if complementary
			if( i1<j1 && i2<j2 && energy.areComplementary(j1,j2) ) {

				// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
				for (k1=j1; k1-- > i1; ) {
					// ensure maximal loop length
					if (j1-k1 > energy.getMaxInternalLoopSize1()+1) break;
					for (k2=j2; k2-- > i2; ) {
						// ensure maximal loop length
						if (j2-k2 > energy.getMaxInternalLoopSize2()+1) break;
						// check if (k1,k2) is valid left boundary
						if ( ! Z_equal(hybridZ_right(k1-i1,k2-i2), 0.0) ) {
							// store value
							hybridZ_right(j1-i1,j2-i2) += energy.getBoltzmannWeight(energy.getE_interLeft(k1,j1,k2,j2)) * hybridZ_right(k1-i1,k2-i2);
						}
					}
				}
			}

			if (Z_isNotINF(hybridZ_right(j1-i1,j2-i2)) && !Z_equal(hybridZ_right(j1-i1,j2-i2),Z_type(0))) {
				// update optimum of right extension if needed
				updateOptRightZ( si1,j1,si2,j2, energy.getE(hybridZ_right(j1-i1,j2-i2) * initZ * seedZ) );
			}

		}
	}

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2dHeuristicSeedExtension::
fillHybridZ_left( const size_t j1, const size_t j2
			, const OutputConstraint & outConstraint )
{

#if INTARNA_IN_DEBUG_MODE
	// check indices
	if (!energy.areComplementary(j1,j2) )
		throw std::runtime_error("PredictorMfeEns2dSeedExtension::fillHybridZ_left("+toString(j1)+","+toString(j2)+",..) are not complementary");
#endif

	// global vars to avoid reallocation
	size_t i1,i2,k1,k2;
	const Z_type seedZ = energy.getBoltzmannWeight(seedHandler.getSeedE(j1, j2));
	const size_t sj1 = j1 + seedHandler.getSeedLength1(j1,j2) -1;
	const size_t sj2 = j2 + seedHandler.getSeedLength2(j1,j2) -1;

	// iterate over all window starts i1 (seq1) and i2 (seq2)
	for (i1=j1; j1-i1 < hybridZ_left.size1(); i1-- ) {
		for (i2=j2; j2-i2 < hybridZ_left.size2(); i2-- ) {

			// init current cell (0 if not just right-most (j1,j2) base pair)
			hybridZ_left(j1-i1,j2-i2) = i1==j1 && i2==j2 ? energy.getBoltzmannWeight(energy.getE_init()) : 0.0;

			// check if complementary (use global sequence indexing)
			if( i1<j1 && i2<j2 && energy.areComplementary(i1,i2) ) {

				// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
				for (k1=i1; k1++ < j1; ) {
					// ensure maximal loop length
					if (k1-i1 > energy.getMaxInternalLoopSize1()+1) break;
					for (k2=i2; k2++ < j2; ) {
						// ensure maximal loop length
						if (k2-i2 > energy.getMaxInternalLoopSize2()+1) break;
						// check if (k1,k2) are valid left boundary
						if ( ! Z_equal(hybridZ_left(j1-k1,j2-k2), 0.0) ) {
							hybridZ_left(j1-i1,j2-i2) += energy.getBoltzmannWeight(energy.getE_interLeft(i1,k1,i2,k2)) * hybridZ_left(j1-k1,j2-k2);
						}
					} // k2
				} // k1
			}

			// correction for left seeds
			if (i1 < j1 && i2 < j2 && seedHandler.isSeedBound(i1, i2) ) {

				// check if seed is to be processed:
				bool substractThisSeed =
						// check if left of anchor seed
									( i1+seedHandler.getSeedLength1(i1,i2)-1 < j1
									&& i2+seedHandler.getSeedLength2(i1,i2)-1 < j2 )
						// check if overlapping with anchor seed
								||	seedHandler.areLoopOverlapping(i1,i2,j1,j2);
				if (substractThisSeed) {

					// iterate seeds in S region
					size_t sj1 = RnaSequence::lastPos, sj2 = RnaSequence::lastPos;
					size_t si1overlap = j1+1, si2overlap = j2+1;
					// find left-most loop-overlapping seed
					while( seedHandler.updateToNextSeed(sj1,sj2
							, i1, std::min(j1,i1+seedHandler.getSeedLength1(i1, i2)-2)
							, i2, std::min(j2,i2+seedHandler.getSeedLength2(i1, i2)-2)) )
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
					if (si1overlap <= j1) {
						// check if right side is non-empty
						if ( ! E_equal(hybridZ_left( j1-si1overlap, j2-si2overlap),0) ) {
							// compute Energy of loop S \ S'
							E_type nonOverlapE = getNonOverlappingEnergy(i1, i2, si1overlap, si2overlap);
							// substract energy.getBoltzmannWeight( nonOverlapE ) * hybridZ_left( si1overlap, si2overlap bis anchor seed [==1 falls gleich])
							Z_type correctionTerm = energy.getBoltzmannWeight( nonOverlapE ) * hybridZ_left( j1-si1overlap, j2-si2overlap);
							hybridZ_left(j1-i1, j2-i2) -= correctionTerm;
							// sanity ensurence
							if (hybridZ_left(j1-i1, j2-i2) < 0) {
								hybridZ_left(j1-i1, j2-i2) = Z_type(0.0);
							}
						}
					} else {
						// get energy of seed
						E_type seedE = seedHandler.getSeedE(i1, i2);
						// if no S'
						// substract seedZ * hybridZ_left(right end seed bis anchor seed)
						Z_type correctionTerm = energy.getBoltzmannWeight(seedE) * hybridZ_left( j1-(i1+seedHandler.getSeedLength1(i1, i2)-1), j2-(i2+seedHandler.getSeedLength2(i1, i2)-1));
						hybridZ_left(j1-i1, j2-i2) -= correctionTerm;
						// sanity ensurence
						if (hybridZ_left(j1-i1, j2-i2) < 0) {
							hybridZ_left(j1-i1, j2-i2) = Z_type(0.0);
						}
					}
				} // substractThisSeed
			}

			// update mfe if needed
			if ( ! Z_equal(hybridZ_left(j1-i1,j2-i2), 0.0) && Z_isNotINF(hybridZ_left(j1-i1,j2-i2)) ) {
				const size_t sl1 = seedHandler.getSeedLength1(j1, j2)-1;
				const size_t sl2 = seedHandler.getSeedLength2(j1, j2)-1;
				const size_t sj1 = j1+sl1;
				const size_t sj2 = j2+sl2;
				PredictorMfeEns::updateOptima( i1,j1opt,i2,j2opt,
						(energy.getE(hybridZ_right(j1opt-sj1, j2opt-sj2))
								+ energy.getE(hybridZ_left(j1-i1,j2-i2))
								+ seedHandler.getSeedE(j1, j2))
								, true );
			}
		}
	}

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2dHeuristicSeedExtension::
getNextBest( Interaction & curBest )
{
	INTARNA_NOT_IMPLEMENTED("PredictorMfeEns2dHeuristicSeedExtension::getNextBest() not implemented yet");
}

//////////////////////////////////////////////////////////////////////////


} // namespace
