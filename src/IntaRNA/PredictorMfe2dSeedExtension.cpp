
#include "IntaRNA/PredictorMfe2dSeedExtension.h"

namespace IntaRNA {

//////////////////////////////////////////////////////////////////////////

PredictorMfe2dSeedExtension::
PredictorMfe2dSeedExtension(
		const InteractionEnergy & energy
		, OutputHandler & output
		, PredictionTracker * predTracker
		, SeedHandler * seedHandlerInstance )
 :
	PredictorMfe2d(energy,output,predTracker)
	, seedHandler(seedHandlerInstance)
	, hybridE_left( 0,0 )
	, hybridE_right( 0,0 )
{
	assert( seedHandler.getConstraint().getBasePairs() > 1 );
}

//////////////////////////////////////////////////////////////////////////

PredictorMfe2dSeedExtension::
~PredictorMfe2dSeedExtension()
{
}

//////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dSeedExtension::
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

	size_t si1 = RnaSequence::lastPos, si2 = RnaSequence::lastPos;
	while( seedHandler.updateToNextSeed(si1,si2
			, 0, range_size1+1-seedHandler.getConstraint().getBasePairs()
			, 0, range_size2+1-seedHandler.getConstraint().getBasePairs()) )
	{
		// get energy and boundaries of seed
		E_type seedE = seedHandler.getSeedE(si1, si2);
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
		hybridE_left.resize( std::min(si1+1, maxMatrixLen1), std::min(si2+1, maxMatrixLen2) );
		fillHybridE_left(si1, si2);

		// ER
		hybridE_right.resize( std::min(range_size1-sj1, maxMatrixLen1), std::min(range_size2-sj2, maxMatrixLen2) );
		fillHybridE_right(sj1, sj2);

		// update Optimum for all boundary combinations
		for (int i1 = 0; i1 < hybridE_left.size1(); i1++) {
			for (int i2 = 0; i2 < hybridE_left.size2(); i2++) {
				if (E_isINF(hybridE_left(i1,i2))) continue;
				const size_t j1max = std::min(maxMatrixLen1-i1, hybridE_right.size1());
				// ensure max interaction length in seq 1
				for (int j1 = 0; j1 < j1max; j1++) {
					assert(sj1+j1-si1+i1 < energy.getAccessibility1().getMaxLength());
					const size_t j2max = std::min(maxMatrixLen2-i2, hybridE_right.size2());
					// ensure max interaction length in seq 2
					for (int j2 = 0; j2 < j2max; j2++) {
						assert(sj2+j2-si2+i2 < energy.getAccessibility2().getMaxLength());
						if (E_isINF(hybridE_right(j1,j2))) continue;
						// no Z update for left extensions to avoid duplicated handling (underestimates Zall)
						updateOptima( si1-i1, sj1+j1, si2-i2, sj2+j2, seedE + hybridE_left(i1,i2) + hybridE_right(j1,j2), true, i1==0 && i2==0 );
					} // j2
				} // i2
			} // j1
		} // i1

	} // si1 / si2

	// report mfe interaction
	reportOptima();

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dSeedExtension::
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

	// determine whether or not lonely base pairs are allowed or if we have to
	// ensure a stacking to the right of the left boundary (i1,i2)
	const size_t noLpShift = outConstraint.noLP ? 1 : 0;
	E_type iStackE = E_type(0);

	// iterate over all window starts j1 (seq1) and j2 (seq2)
	for (size_t l1=0; l1 < hybridE_left.size1(); l1++) {
		for (size_t l2=0; l2 < hybridE_left.size2(); l2++) {
			i1 = si1-l1;
			i2 = si2-l2;

			// referencing cell access
			E_type & curE = hybridE_left(l1,l2);
			// init current cell (e_init if just left (i1,i2) base pair; assuming seed is internally stacked on the left end if noLP)
			curE = (i1==si1 && i2==si2) ? energy.getE_init() : E_INF;
			// check if complementary
			if( i1<si1 && i2<si2 && energy.areComplementary(i1,i2) ) {

				// right-stacking of i if no-LP
				if (outConstraint.noLP) {
					// skip if no stacking possible
					if (!energy.areComplementary(i1+noLpShift,i2+noLpShift)) {
						continue;
					}
					// get stacking energy to avoid recomputation in recursion below
					iStackE = energy.getE_interLeft(i1,i1+noLpShift,i2,i2+noLpShift);
					// check just stacked
					curE = std::min( curE, iStackE + hybridE_left(l1-noLpShift,l2-noLpShift));
				}

				// check all combinations of decompositions into (i1,i2)..(k1,k2)-(si1,si2)
				for (k1=i1+noLpShift; k1++ < si1; ) {
					// ensure maximal loop length
					if (k1-i1-noLpShift > energy.getMaxInternalLoopSize1()+1) break;
					for (k2=i2+noLpShift; k2++ < si2; ) {
						// ensure maximal loop length
						if (k2-i2-noLpShift > energy.getMaxInternalLoopSize2()+1) break;
						// check if (k1,k2) are valid left boundary
						if ( E_isNotINF( hybridE_left(si1-k1,si2-k2) ) ) {
							curE = std::min( curE,
									(iStackE
											+ energy.getE_interLeft(i1+noLpShift,k1,i2+noLpShift,k2)
											+ hybridE_left(si1-k1,si2-k2) )
							);
						}
					} // k2
				} // k1
			} // i is valid right end
		} // i2
	} // i1

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dSeedExtension::
fillHybridE_right( const size_t sj1, const size_t sj2 )
{
	// temporary access
	const OutputConstraint & outConstraint = output.getOutputConstraint();
#if INTARNA_IN_DEBUG_MODE
	// check indices
	if (!energy.areComplementary(sj1,sj2) )
		throw std::runtime_error("PredictorMfe2dSeedExtension::fillHybridE_right("+toString(sj1)+","+toString(sj2)+",..) are not complementary");
#endif

	// global vars to avoid reallocation
	size_t j1,j2,k1,k2;

	// determine whether or not lonely base pairs are allowed or if we have to
	// ensure a stacking to the right of the left boundary (i1,i2)
	const size_t noLpShift = outConstraint.noLP ? 1 : 0;
	E_type iStackE = E_type(0);

	// iterate over all window starts j1 (seq1) and j2 (seq2)
	for (j1=sj1; j1-sj1 < hybridE_right.size1(); j1++) {
		// screen for right boundaries in seq2
		for (j2=sj2; j2-sj2 < hybridE_right.size2(); j2++) {

			// referencing access to cell
			E_type & curE = hybridE_right(j1-sj1,j2-sj2);
			// init current cell (0 if just left (i1,i2) base pair)
			curE = (sj1==j1 && sj2==j2) ? 0 : E_INF;

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
						curE = std::min( curE, iStackE + hybridE_right(0,0) );
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
						curE = std::min( curE,
								(hybridE_right(k1-sj1,k2-sj2)
										+ energy.getE_interLeft(k1,j1-noLpShift,k2,j2-noLpShift)
										+ iStackE )
								);
					}
				} // k2
			  } // k1
			} // j is valid right bound
		} // j2
	} // j1

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dSeedExtension::
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
		throw std::runtime_error("PredictorMfe2dSeedExtension::traceBack() : given interaction does not contain boundaries only");
	}
#endif

	// check for single interaction
	if (interaction.basePairs.at(0).first == interaction.basePairs.at(1).first) {
		// delete second boundary (identical to first)
		interaction.basePairs.resize(1);
		// update done
		return;
	}

#if INTARNA_IN_DEBUG_MODE
	// sanity checks
	if ( ! interaction.isValid() ) {
		throw std::runtime_error("PredictorMfe2dSeedExtension::traceBack() : given interaction not valid");
	}
#endif

	// ensure sorting
	interaction.sort();
	// get indices in hybridE for boundary base pairs
	size_t	i1 = energy.getIndex1(interaction.basePairs.at(0)),
			j1 = energy.getIndex1(interaction.basePairs.at(1)),
			i2 = energy.getIndex2(interaction.basePairs.at(0)),
			j2 = energy.getIndex2(interaction.basePairs.at(1))
			;

#if INTARNA_IN_DEBUG_MODE
	// check if intervals are larger enough to contain a seed
	if (std::min(j1-i1,j2-i2)+1 < seedHandler.getConstraint().getBasePairs()) {
		// no seed possible, abort computation
		throw std::runtime_error("PredictorMfe2dSeedExtension::traceBack() : given boundaries "+toString(interaction)+" can not hold a seed of "+toString(seedHandler.getConstraint().getBasePairs())+" base pairs");
	}
#endif

	E_type fullE = interaction.energy;

	// determine whether or not lonely base pairs are allowed or if we have to
	// ensure a stacking to the right of the left boundary (i1,i2)
	const size_t noLpShift = outConstraint.noLP ? 1 : 0;
	E_type iStackE = E_type(0);

	size_t si1 = RnaSequence::lastPos, si2 = RnaSequence::lastPos;
	while( seedHandler.updateToNextSeed(si1,si2
			, i1,j1+1-seedHandler.getConstraint().getBasePairs()
			, i2,j2+1-seedHandler.getConstraint().getBasePairs() ) )
	{
			E_type seedE = seedHandler.getSeedE(si1, si2);
			const size_t sl1 = seedHandler.getSeedLength1(si1, si2);
			const size_t sl2 = seedHandler.getSeedLength2(si1, si2);
			const size_t sj1 = si1+sl1-1;
			const size_t sj2 = si2+sl2-1;
			const size_t maxMatrixLen1 = energy.getAccessibility1().getMaxLength()-sl1+1;
			const size_t maxMatrixLen2 = energy.getAccessibility2().getMaxLength()-sl2+1;

			hybridE_left.resize( std::min(si1+1, maxMatrixLen1), std::min(si2+1, maxMatrixLen2) );
			fillHybridE_left( si1, si2 );
			hybridE_right.resize( std::min(j1-sj1+1, maxMatrixLen1), std::min(j2-sj2+1, maxMatrixLen2) );
			fillHybridE_right( sj1, sj2 );

			if ( E_equal( fullE,
					(energy.getE(i1, j1, i2, j2, seedE + hybridE_left( si1-i1, si2-i2 ) + hybridE_right( j1-sj1, j2-sj2 )))))
			{
				// found seed -> traceback
				// the currently traced value for i1-si1, i2-si2
				E_type curE = hybridE_left(si1-i1, si2-i2);

				// trace back left
				while( i1+noLpShift < si1 ) {

					// right-stacking of i if no-LP
					if (outConstraint.noLP) {
						// assure stacking is possible
						assert(energy.areComplementary(i1+noLpShift,i2+noLpShift));
						// check if just stacking loop
						if ( si1-i1==2 && si2-i2==2 ) {
							interaction.basePairs.push_back( energy.getBasePair(i1+noLpShift,i2+noLpShift) );
							break;
						}
						// get stacking energy to avoid recomputation in recursion below
						iStackE = energy.getE_interLeft(i1,i1+noLpShift,i2,i2+noLpShift);
					}

					// check if just internal loop
					if ( E_equal( curE, (iStackE + energy.getE_interLeft(i1+noLpShift,si1,i2+noLpShift,si2) + hybridE_left(0,0)) ) )
					{
						// store noLP base pair
						if ( outConstraint.noLP ) {
							interaction.basePairs.push_back( energy.getBasePair(i1+noLpShift,i2+noLpShift) );
						}
						break;
					}
					// check all interval splits
					if ( (si1-i1-noLpShift) > 1 && (si2-i2-noLpShift) > 1) {

						bool traceNotFound = true;
						// temp variables
						size_t k1=i1+noLpShift,k2=i2+noLpShift;
						if (outConstraint.noLP && E_equal( curE, (iStackE + hybridE_left(si1-k1,si2-k2)) )) {
							traceNotFound = false;
							interaction.basePairs.push_back( energy.getBasePair(k1,k2) );
							// trace right part of stack
							i1=k1;
							i2=k2;
							curE = hybridE_left(si1-i1,si2-i2);
						}
						// check all combinations of decompositions into (i1,i2)..(k1,k2)-(si1,si2)
						const size_t k1max = std::min(si1-1,i1+noLpShift+energy.getMaxInternalLoopSize1()+1);
						const size_t k2max = std::min(si2-1,i2+noLpShift+energy.getMaxInternalLoopSize2()+1);
						for (k1++; traceNotFound && k1<=k1max; k1++) {
						for (k2++; traceNotFound && k2<=k2max; k2++) {
							// check if (k1,k2) are valid left boundary
							if ( E_isNotINF( hybridE_left(si1-k1,si2-k2) ) ) {
								if ( E_equal( curE,
										(iStackE + energy.getE_interLeft(i1+noLpShift,k1,i2+noLpShift,k2) + hybridE_left(si1-k1,si2-k2)) ) )
								{
									// stop searching
									traceNotFound = false;
									// store no-LP base pair right of i
									if (outConstraint.noLP) {
										interaction.basePairs.push_back( energy.getBasePair(i1+noLpShift,i2+noLpShift) );
									}
									// store splitting base pair
									interaction.basePairs.push_back( energy.getBasePair(k1,k2) );
									// trace right part of split
									i1=k1;
									i2=k2;
									curE = hybridE_left(si1-i1,si2-i2);
								}
							}
						}
						}
						if (traceNotFound) {
							LOG(ERROR) <<"left-not-found:  i "<<i1<<","<<i2<<"  si "<<si1<<","<<si2;
							throw std::runtime_error("trace not found");
						}
					}

				} // traceback left

				// trace seed
				if (si1 > i1 && si2 > i2) {
					interaction.basePairs.push_back( energy.getBasePair(si1,si2) );
				}
				seedHandler.traceBackSeed( interaction, si1, si2 );
				if (sj1 < j1 && sj2 < j2) {
					interaction.basePairs.push_back( energy.getBasePair(sj1,sj2) );
				}

				// the currently traced value for sj1-j1, sj2-j2
				curE = hybridE_right(j1-sj1,j2-sj2);

				// trace back right
				while( j1-noLpShift > sj1 ) {

					// left-stacking of j if no-LP
					if (outConstraint.noLP) {
						// assure stacking is possible
						assert (energy.areComplementary(j1-noLpShift,j2-noLpShift));
						// handle stacking only
						if (j1-sj1==2 && j2-sj2==2) {
							interaction.basePairs.push_back( energy.getBasePair(j1-noLpShift,j2-noLpShift) );
							break;
						}
						// get stacking energy to avoid recomputation in recursion below
						iStackE = energy.getE_interLeft(j1-noLpShift,j1,j2-noLpShift,j2);
					}

					// check if just internal loop
					if ( E_equal( curE, (hybridE_right(0,0) + energy.getE_interLeft(sj1,j1-noLpShift,sj2,j2-noLpShift) + iStackE) ) )
					{
						// store no-LP base pair left of j
						if (outConstraint.noLP) {
							interaction.basePairs.push_back( energy.getBasePair(j1-noLpShift,j2-noLpShift) );
						}
						break;
					}
					// check all interval splits
					if ( (j1-noLpShift-sj1) > 1 && (j2-noLpShift-sj2) > 1) {
						// temp variables
						size_t k1,k2;
						bool traceNotFound = true;
						// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
						for (k1=j1-noLpShift-1; traceNotFound && k1 > sj1 && (j1-noLpShift-k1 <= energy.getMaxInternalLoopSize1()+1); k1--) {
						for (k2=j2-noLpShift-1; traceNotFound && k2 > sj2 && (j2-noLpShift-k2 <= energy.getMaxInternalLoopSize2()+1); k2--) {
							// check if (k1,k2) are valid left boundary
							if ( E_isNotINF( hybridE_right(k1-sj1,k2-sj2) ) ) {
								if ( E_equal( curE,
										(iStackE + energy.getE_interLeft(k1,j1-noLpShift,k2,j2-noLpShift) + hybridE_right(k1-sj1,k2-sj2)) ) )
								{
									// stop searching
									traceNotFound = false;
									// store no-LP base pair left of j
									if (outConstraint.noLP) {
										interaction.basePairs.push_back( energy.getBasePair(j1-noLpShift,j2-noLpShift) );
									}
									// store splitting base pair
									interaction.basePairs.push_back( energy.getBasePair(k1,k2) );
									// trace right part of split
									j1=k1;
									j2=k2;
									curE = hybridE_right(j1-sj1,j2-sj2);
								}
							}
						}
						}
						if (traceNotFound) {
							LOG(ERROR) <<"right-not-found:  sji "<<sj1<<","<<sj2<<"  j "<<j1<<","<<j2;
							throw std::runtime_error("trace not found");
						}
					}
				}  // traceback right

				interaction.sort();
				seedHandler.addSeeds( interaction );

				// stop searching for seeds
				return;
			} // is valid seed for this interaction

	} // si1 / si2

}

////////////////////////////////////////////////////////////////////////////


} // namespace
