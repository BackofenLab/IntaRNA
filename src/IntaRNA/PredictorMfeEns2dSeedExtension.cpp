
#include "IntaRNA/PredictorMfeEns2dSeedExtension.h"

namespace IntaRNA {

//////////////////////////////////////////////////////////////////////////

PredictorMfeEns2dSeedExtension::
PredictorMfeEns2dSeedExtension(
		const InteractionEnergy & energy
		, OutputHandler & output
		, PredictionTracker * predTracker
		, SeedHandler * seedHandlerInstance )
 :
	PredictorMfeEns(energy,output,predTracker)
	, seedHandler(seedHandlerInstance)
	, hybridZ_right( 0,0 )
{
	assert( seedHandler.getConstraint().getBasePairs() > 1 );
}

//////////////////////////////////////////////////////////////////////////

PredictorMfeEns2dSeedExtension::
~PredictorMfeEns2dSeedExtension()
{
}

//////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2dSeedExtension::
predict( const IndexRange & r1, const IndexRange & r2 )
{
	// temporary access
	const OutputConstraint & outConstraint = output.getOutputConstraint();
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) <<"predicting mfe interactions with seed in O(n^2) space and O(n^4) time..."; }
	// measure timing
	TIMED_FUNC_IF(timerObj,VLOG_IS_ON(9));

	// suboptimal setup check
	if (outConstraint.reportMax>1 && outConstraint.reportOverlap != OutputConstraint::ReportOverlap::OVERLAP_BOTH) {
		throw std::runtime_error("PredictorMfeEns2dSeedExtension : the enumeration of non-overlapping suboptimal interactions is not supported in this prediction mode");
	}

	// no-LP setup check
	if (outConstraint.noLP) {
		INTARNA_NOT_IMPLEMENTED("PredictorMfeEns2dSeedExtension : prediction without lonely base pairs is not implemented yet");
	}

#if INTARNA_IN_DEBUG_MODE
	// check indices
	if (!(r1.isAscending() && r2.isAscending()) )
		throw std::runtime_error("PredictorMfeEns2dSeedExtension::predict("+toString(r1)+","+toString(r2)+") is not sane");
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

		// EL
		hybridZ_left.resize( std::min(si1+1, maxMatrixLen1), std::min(si2+1, maxMatrixLen2) );
		fillHybridZ_left(si1, si2);

		// ER
		hybridZ_right.resize( std::min(range_size1-sj1, maxMatrixLen1), std::min(range_size2-sj2, maxMatrixLen2) );
		fillHybridZ_right(sj1, sj2);

		// updateZ for all boundary combinations
		for (size_t l1 = 0; l1<hybridZ_left.size1(); l1++) {
			for (size_t l2 = 0; l2< hybridZ_left.size2(); l2++) {
				// check complementarity of boundary
				if ( Z_equal(hybridZ_left(l1,l2), 0.0) ) continue;
				// ensure max interaction length in seq 1
				const size_t r1max = std::min(hybridZ_right.size1(), energy.getAccessibility1().getMaxLength()+1-sl1-l1);
				for (size_t r1 = 0; r1 < r1max ; r1++) {
					// check interaction length
					assert(sj1+r1-si1+l1 < energy.getAccessibility1().getMaxLength());
					// ensure max interaction length in seq 2
					for (size_t r2 = 0; r2 < hybridZ_right.size2() ; r2++) {
						// check interaction length
						if (sj2+r2-si2+l2 >= energy.getAccessibility2().getMaxLength()) break;
						// check complementarity of boundary
						if (Z_equal(hybridZ_right(r1,r2),0.0)) continue;
						// compute partition function given the current seed
						updateZ(si1-l1, sj1+r1, si2-l2, sj2+r2, hybridZ_left(l1,l2) * seedZ * hybridZ_right(r1,r2), true);
					} // r2
				} // l2
			} // r1
		} // l1

	} // si1 / si2

	// report mfe interaction
	reportOptima();

}

//////////////////////////////////////////////////////////////////////////

E_type
PredictorMfeEns2dSeedExtension::
getNonOverlappingEnergy( const size_t si1, const size_t si2, const size_t si1p, const size_t si2p ) {

#if INTARNA_IN_DEBUG_MODE
	// check indices
	if( !seedHandler.isSeedBound(si1,si2) )
		throw std::runtime_error("PredictorMfeEns2dSeedExtension::getNonOverlappingEnergy( si "+toString(si1)+","+toString(si2)+",..) is no seed bound");
	if( !seedHandler.isSeedBound(si1p,si2p) )
		throw std::runtime_error("PredictorMfeEns2dSeedExtension::getNonOverlappingEnergy( sip "+toString(si1p)+","+toString(si2p)+",..) is no seed bound");
#endif

	// sanity check
	if( si1 == si1p || ! seedHandler.areLoopOverlapping(si1,si2,si1p,si2p) ) {
		return E_type(0);
	}

	// trace S
	Interaction interaction = Interaction(energy.getAccessibility1().getSequence(), energy.getAccessibility2().getAccessibilityOrigin().getSequence());
	interaction.basePairs.push_back( energy.getBasePair(si1, si2) );
	seedHandler.traceBackSeed( interaction, si1, si2 );

	E_type fullE = 0;
	size_t k1old = si1, k2old = si2;
	for (size_t i = 1; i < interaction.basePairs.size(); i++) {
		// get index of current base pair
		size_t k1 = energy.getIndex1(interaction.basePairs[i]);
		// check if overlap done
		if (k1 > si1p) break;
		size_t k2 = energy.getIndex2(interaction.basePairs[i]);
		// add hybridization energy
		fullE += energy.getE_interLeft(k1old,k1,k2old,k2);
		// store
		k1old = k1;
		k2old = k2;
	}
	return fullE;
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2dSeedExtension::
fillHybridZ_left( const size_t j1, const size_t j2 )
{
	// temporary access
	const OutputConstraint & outConstraint = output.getOutputConstraint();
#if INTARNA_IN_DEBUG_MODE
	// check indices
	if (!energy.areComplementary(j1,j2) )
		throw std::runtime_error("PredictorMfeEns2dSeedExtension::fillHybridZ_left("+toString(j1)+","+toString(j2)+",..) are not complementary");
#endif

	// global vars to avoid reallocation
	size_t i1,i2,k1,k2;

	// iterate over all window starts i1 (seq1) and i2 (seq2)
	for (i1=j1; j1-i1 < hybridZ_left.size1(); i1-- ) {
		for (i2=j2; j2-i2 < hybridZ_left.size2(); i2-- ) {

			// referencing cell access
			Z_type & curZ = hybridZ_left(j1-i1,j2-i2);

			// init current cell (0 if not just right-most (j1,j2) base pair)
			curZ = (i1==j1 && i2==j2) ? energy.getBoltzmannWeight(energy.getE_init()) : 0.0;

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
							curZ += energy.getBoltzmannWeight(energy.getE_interLeft(i1,k1,i2,k2))
									* hybridZ_left(j1-k1,j2-k2);
						}
					} // k2
				} // k1

				// correction for left seeds
				if ( i1<j1 && i2<j2 && seedHandler.isSeedBound(i1, i2) ) {

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
								Z_type correctionTerm = energy.getBoltzmannWeight( nonOverlapE )
														* hybridZ_left( j1-si1overlap, j2-si2overlap);
								curZ -= correctionTerm;
								// sanity ensurence
								if (curZ < 0) {
									curZ = Z_type(0.0);
								}
							}
						} else {
							// get energy of seed
							E_type seedE = seedHandler.getSeedE(i1, i2);
							// if no S'
							// substract seedZ * hybridZ_left(right end seed bis anchor seed)
							Z_type correctionTerm = energy.getBoltzmannWeight(seedE)
									* hybridZ_left( j1-(i1+seedHandler.getSeedLength1(i1, i2)-1)
												  , j2-(i2+seedHandler.getSeedLength2(i1, i2)-1));
							curZ -= correctionTerm;
							// sanity ensurence
							if (curZ < 0) {
								curZ = Z_type(0.0);
							}
						}
					} // substractThisSeed
				}

			} // complementary

		} // i2
	} // i1

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2dSeedExtension::
fillHybridZ_right( const size_t i1, const size_t i2 )
{
	// temporary access
	const OutputConstraint & outConstraint = output.getOutputConstraint();
#if INTARNA_IN_DEBUG_MODE
	// check indices
	if (!energy.areComplementary(i1,i2) )
		throw std::runtime_error("PredictorMfeEns2dSeedExtension::fillHybridZ_right("+toString(i1)+","+toString(i2)+",..) are not complementary");
#endif


	// global vars to avoid reallocation
	size_t j1,j2,k1,k2;

	// iterate over all window ends j1 (seq1) and j2 (seq2)
	for (j1=i1; j1-i1 < hybridZ_right.size1(); j1++ ) {
		for (j2=i2; j2-i2 < hybridZ_right.size2(); j2++ ) {

			// referencing cell access
			Z_type & curZ = hybridZ_right(j1-i1,j2-i2);

			// init partition function for current cell -> (i1,i2) are complementary per definition
			curZ = i1==j1 && i2==j2 ? energy.getBoltzmannWeight(0.0) : 0.0;

			// check if complementary free base pair
			if( i1<j1 && i2<j2 && energy.areComplementary(j1,j2) ) {

				// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
				for (k1=j1; k1-- > i1; ) {
					// ensure maximal loop length
					if (j1-k1 > energy.getMaxInternalLoopSize1()+1) break;
					for (k2=j2; k2-- > i2; ) {
						// ensure maximal loop length
						if (j2-k2 > energy.getMaxInternalLoopSize2()+1) break;
						// check if (k1,k2) are valid left boundary
						if ( ! Z_equal(hybridZ_right(k1-i1,k2-i2), 0.0) ) {
							// update partition function
							curZ += energy.getBoltzmannWeight(energy.getE_interLeft(k1,j1,k2,j2)) * hybridZ_right(k1-i1,k2-i2);
						}
					} // k2
				} // k1
			}
		}
	}

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2dSeedExtension::
traceBack( Interaction & interaction )
{

	// forward tracing
	PredictorMfeEns::traceBack(interaction);

	// add seeds in region
	seedHandler.addSeeds( interaction );

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2dSeedExtension::
getNextBest( Interaction & curBest )
{
	INTARNA_NOT_IMPLEMENTED("PredictorMfeEns2dSeedExtension::getNextBest() not implemented yet");
}

//////////////////////////////////////////////////////////////////////////


} // namespace
