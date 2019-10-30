
#include "IntaRNA/PredictorMfe2dHeuristicSeedExtension.h"

namespace IntaRNA {

//////////////////////////////////////////////////////////////////////////

PredictorMfe2dHeuristicSeedExtension::
PredictorMfe2dHeuristicSeedExtension(
		const InteractionEnergy & energy
		, OutputHandler & output
		, PredictionTracker * predTracker
		, SeedHandler * seedHandlerInstance )
 :
	PredictorMfe2dSeedExtension(energy,output,predTracker,seedHandlerInstance)
	, E_right_opt(E_INF)
	, j1opt(0)
	, j2opt(0)
{
	assert( seedHandler.getConstraint().getBasePairs() > 1 );
}

//////////////////////////////////////////////////////////////////////////

PredictorMfe2dHeuristicSeedExtension::
~PredictorMfe2dHeuristicSeedExtension()
{
}

//////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dHeuristicSeedExtension::
predict( const IndexRange & r1, const IndexRange & r2 )
{
	// temporary access
	const OutputConstraint & outConstraint = output.getOutputConstraint();
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) <<"predicting mfe interactions with seed-extension heuristically in O(n^2) space and O(n^2) time..."; }
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

	// initialize mfe interaction for updates
	initOptima();

	// compute seed interactions for whole range
	// and check if any seed possible
	if (seedHandler.fillSeed( 0, range_size1-1, 0, range_size2-1 ) == 0) {
		// trigger empty interaction reporting
		reportOptima();
		// stop computation
		return;
	}

	size_t si1 = RnaSequence::lastPos, si2 = RnaSequence::lastPos;
	while( seedHandler.updateToNextSeed(si1,si2
			, 0, range_size1+1-seedHandler.getConstraint().getBasePairs()
			, 0, range_size2+1-seedHandler.getConstraint().getBasePairs()) )
	{
		E_type seedE = seedHandler.getSeedE(si1, si2);
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

		// ER : update for seed + ER
		hybridE_right.resize( std::min(range_size1-sj1, maxMatrixLen1), std::min(range_size2-sj2, maxMatrixLen2) );
		fillHybridE_right(sj1, sj2, si1, si2);

		// ensure there is a valid right-extension
		if (!output.getOutputConstraint().noGUend || (!E_isINF(E_right_opt) || !energy.isGU(sj1,sj2))) {
			// EL
			hybridE_left.resize( std::min(si1+1, maxMatrixLen1), std::min(si2+1, maxMatrixLen2) );
			fillHybridE_left(si1, si2);
		}

	} // si1 / si2

	// report mfe interaction
	reportOptima();
}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dHeuristicSeedExtension::
fillHybridE_right( const size_t sj1, const size_t sj2
			, const size_t si1, const size_t si2 )
{
	// temporary access
	const OutputConstraint & outConstraint = output.getOutputConstraint();

	// global vars to avoid reallocation
	size_t j1,j2,k1,k2;

	// seed energy
	const E_type seedE = seedHandler.getSeedE(si1, si2);

	// determine whether or not lonely base pairs are allowed or if we have to
	// ensure a stacking to the right of the left boundary (i1,i2)
	const size_t noLpShift = outConstraint.noLP ? 1 : 0;
	E_type iStackE = E_type(0);

	// current minimal value
	E_type curMinE = E_INF;
	// iterate over all window starts j1 (seq1) and j2 (seq2)
	for (j1=sj1; j1-sj1 < hybridE_right.size1(); j1++) {
		// screen for right boundaries in seq2
		for (j2=sj2; j2-sj2 < hybridE_right.size2(); j2++) {

			// referencing cell access
			E_type & curMinE = hybridE_right(j1-sj1,j2-sj2);

			// init current cell (0 if just left (i1,i2) base pair)
			curMinE = (sj1==j1 && sj2==j2) ? 0 : E_INF;

			// check if complementary
			if( sj1<j1 && sj2<j2 && energy.areComplementary(j1,j2) ) {

				// left-stacking of j if no-LP
				if (outConstraint.noLP) {
					// skip if no stacking possible
					if (!energy.areComplementary(j1-noLpShift,j2-noLpShift)) {
						continue;
					}
					// get stacking energy to avoid recomputation in recursion below
					iStackE = energy.getE_interLeft(j1-noLpShift,j1,j2-noLpShift,j2);
					// check just stacked seed extension
					if (j1-noLpShift==sj1 && j2-noLpShift==sj2) {
						curMinE = std::min( curMinE, iStackE + hybridE_right(0,0) );
					}
				}

				// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
				for (k1=j1-noLpShift; k1-- > sj1; ) {
					// ensure maximal loop length
					if (j1-noLpShift-k1 > energy.getMaxInternalLoopSize1()+1) break;
				for (k2=j2-noLpShift; k2-- > sj2; ) {
					// ensure maximal loop length
					if (j2-noLpShift-k2 > energy.getMaxInternalLoopSize2()+1) break;
					// check if (k1,k2) are valid left boundary
					if ( E_isNotINF( hybridE_right(k1-sj1,k2-sj2) ) ) {
						curMinE = std::min( curMinE,
								(hybridE_right(k1-sj1,k2-sj2) // left part
								+ energy.getE_interLeft(k1,j1-noLpShift,k2,j2-noLpShift) // loop
								+ iStackE) // right stack if no-LP
								);
					}
				}
				}

				// update mfe if needed
				if (E_isNotINF(curMinE)) {
					// seedE + rightE + E_init (not covered by rightE)
					updateOptimalRightExt( si1,j1,si2,j2, seedE + curMinE + energy.getE_init(), true );
					// update mfe for seed+rightExt
					updateOptima( si1,j1,si2,j2, seedE + curMinE + energy.getE_init(),true,true);
				}
			}
		}
	}

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dHeuristicSeedExtension::
fillHybridE_left( const size_t si1, const size_t si2 )
{
	// temporary access
	const OutputConstraint & outConstraint = output.getOutputConstraint();

#if INTARNA_IN_DEBUG_MODE
	// check indices
	if (!energy.areComplementary(si1,si2) )
		throw std::runtime_error("PredictorMfe2dSeedExtension::fillHybridE_left("+toString(si1)+","+toString(si2)+",..) are not complementary");
#endif

	// global vars to avoid reallocation
	size_t i1,i2,k1,k2;

	const E_type seedE = seedHandler.getSeedE(si1, si2);
	const size_t sl1 = seedHandler.getSeedLength1(si1, si2);
	const size_t sl2 = seedHandler.getSeedLength2(si1, si2);
	const size_t sj1 = si1+sl1-1;
	const size_t sj2 = si2+sl2-1;
	const E_type rightOptE = hybridE_right(j1opt-sj1, j2opt-sj2);

	// determine whether or not lonely base pairs are allowed or if we have to
	// ensure a stacking to the right of the left boundary (i1,i2)
	const size_t noLpShift = outConstraint.noLP ? 1 : 0;
	E_type iStackE = E_type(0);

	// iterate over all window starts
	for (size_t l1=0; l1 < hybridE_left.size1(); l1++) {
		for (size_t l2=0; l2 < hybridE_left.size2(); l2++) {
			i1 = si1-l1;
			i2 = si2-l2;

			// referencing cell access
			E_type & curMinE = hybridE_left(l1,l2);
			// init cell
			curMinE = (i1==si1 && i2==si2) ? energy.getE_init() : E_INF;
			// check if complementary
			if( i1<si1 && i2<si2 && energy.areComplementary(i1,i2) ) {

				// left-stacking of j if no-LP
				if (outConstraint.noLP) {
					// skip if no stacking possible
					if (!energy.areComplementary(i1+noLpShift,i2+noLpShift)) {
						continue;
					}
					// get stacking energy to avoid recomputation in recursion below
					iStackE = energy.getE_interLeft(i1,i1+noLpShift,i2,i2+noLpShift);
					// check just stacked
					curMinE = std::min( curMinE, iStackE + hybridE_left(l1-noLpShift,l2-noLpShift));
				}

				// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
				for (k1=i1+noLpShift; k1++ < si1; ) {
					// ensure maximal loop length
					if (k1-i1-noLpShift > energy.getMaxInternalLoopSize1()+1) break;
				for (k2=i2+noLpShift; k2++ < si2; ) {
					// ensure maximal loop length
					if (k2-i2-noLpShift > energy.getMaxInternalLoopSize2()+1) break;
					// check if (k1,k2) are valid left boundary
					if ( E_isNotINF( hybridE_left(si1-k1,si2-k2) ) ) {
						curMinE = std::min( curMinE,
								(iStackE // i stacking if no-LP
										+ energy.getE_interLeft(i1+noLpShift,k1,i2+noLpShift,k2) // loop
										+ hybridE_left(si1-k1,si2-k2) ) // right part up to seed
								);
					}
				} // k2
			  } // k1
			}
			// update mfe if needed
			if ( E_isNotINF( curMinE ) ) {
				// leftE (incl. E_init) + seedE
				updateOptima( i1,sj1,i2,sj2, curMinE + seedE, true, i1==si1 && i2==si2 );

				// check if right opt != right seed boundary
				// check for max interaction length
				if (j1opt != sj1
					&& (j1opt+1-i1)<=energy.getAccessibility1().getMaxLength()
					&& (j2opt+1-i2)<=energy.getAccessibility2().getMaxLength() )
				{
					// leftE (incl. E_init) + seedE + rightOptE
					// no Z update for left extensions to avoid duplicated handling (underestimates Zall)
					updateOptima( i1,j1opt,i2,j2opt, curMinE + seedE + rightOptE, true, false );
				}
			}
		}
	}

}

////////////////////////////////////////////////////////////////////////////


} // namespace
