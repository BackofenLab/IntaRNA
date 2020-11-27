
#include "IntaRNA/PredictorMfe2dHeuristicSeed.h"

#include <stdexcept>

namespace IntaRNA {

////////////////////////////////////////////////////////////////////////////

PredictorMfe2dHeuristicSeed::
PredictorMfe2dHeuristicSeed( const InteractionEnergy & energy
		, OutputHandler & output
		, PredictionTracker * predTracker
		, SeedHandler * seedHandlerInstance
		)
 : PredictorMfe2dHeuristic(energy,output,predTracker)
	, seedHandler( seedHandlerInstance )
{
}


////////////////////////////////////////////////////////////////////////////

PredictorMfe2dHeuristicSeed::
~PredictorMfe2dHeuristicSeed()
{
	// clean up
}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dHeuristicSeed::
predict( const IndexRange & r1
		, const IndexRange & r2 )
{
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) <<"predicting mfe interactions with seed heuristically in O(n^2) space and time..."; }
	// measure timing
	TIMED_FUNC_IF(timerObj,VLOG_IS_ON(9));

#if INTARNA_IN_DEBUG_MODE
	// check indices
	if (!(r1.isAscending() && r2.isAscending()) )
		throw std::runtime_error("PredictorMfe2dHeuristicSeed::predict("+toString(r1)+","+toString(r2)+") is not sane");
#endif

	// set index offset
	energy.setOffset1(r1.from);
	energy.setOffset2(r2.from);
	seedHandler.setOffset1(r1.from);
	seedHandler.setOffset2(r2.from);

	const size_t hybridEsize1 = std::min( energy.size1()
			, (r1.to==RnaSequence::lastPos?energy.size1()-1:r1.to)-r1.from+1 );
	const size_t hybridEsize2 = std::min( energy.size2()
			, (r2.to==RnaSequence::lastPos?energy.size2()-1:r2.to)-r2.from+1 );

	// compute seed interactions for whole range
	// and check if any seed possible
	if (seedHandler.fillSeed( 0, hybridEsize1-1, 0, hybridEsize2-1 ) == 0) {
		// trigger empty interaction reporting
		initOptima();
		reportOptima();
		// stop computation
		return;
	}

	// resize matrices
	hybridE.resize( hybridEsize1, hybridEsize2 );
	hybridE_seed.resize( hybridE.size1(), hybridE.size2() );

	// reinit mfe for later updates with final information
	initOptima();

	// fill matrices and update optima
	fillHybridE();

	// report mfe interaction
	reportOptima();

}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dHeuristicSeed::
fillHybridE()
{
	// temporary access
	const OutputConstraint & outConstraint = output.getOutputConstraint();
	// compute entries
	// current minimal value
	E_type curE = E_INF, curEtotal = E_INF, curCellEtotal = E_INF;
	E_type curEseedtotal = E_INF, curCellSeedEtotal = E_INF;
	E_type curEloop = E_INF;
	size_t i1,i2,w1,w2;

	// determine whether or not lonely base pairs are allowed or if we have to
	// ensure a stacking to the right of the left boundary (i1,i2)
	const size_t noLpShift = outConstraint.noLP ? 1 : 0;
	E_type iStackE = E_type(0);

	BestInteractionE * curCell = NULL, *curCellSeed = NULL;
	const BestInteractionE * rightExt = NULL;
	// iterate (decreasingly) over all left interaction starts
	for (i1=hybridE.size1(); i1-- > 0;) {
		for (i2=hybridE.size2(); i2-- > 0;) {
			// direct cell access
			curCell = &(hybridE(i1,i2));
			curCellSeed = &(hybridE_seed(i1,i2));
			// init as invalid boundary
			*curCell = BestInteractionE(E_INF, RnaSequence::lastPos, RnaSequence::lastPos);
			*curCellSeed = BestInteractionE(E_INF, RnaSequence::lastPos, RnaSequence::lastPos);

			// check if positions can form interaction
			if (	energy.isAccessible1(i1)
					&& energy.isAccessible2(i2)
					&& energy.areComplementary(i1,i2) )
			{

				// no lp allowed
				if (noLpShift != 0) {
					// check if right-side stacking of (i1,i2) is possible
					if ( i1+noLpShift < energy.size1()
						&& i2+noLpShift < energy.size2()
						&& energy.isAccessible1(i1+noLpShift)
						&& energy.isAccessible2(i2+noLpShift)
						&& energy.areComplementary(i1+noLpShift,i2+noLpShift))
					{
						// get stacking term to avoid recomputation
						iStackE = energy.getE_interLeft(i1,i1+noLpShift,i2,i2+noLpShift);
					} else {
						// skip further processing, since no stacking possible
						iStackE = E_INF;
					}
				}

				// set to interaction initiation with according boundary
				// if valid iStackE value
				// if valid right boundary
				if (E_isNotINF(iStackE)
						&& (!outConstraint.noGUend || !energy.isGU(i1+noLpShift,i2+noLpShift)))
				{
					*curCell = BestInteractionE(iStackE+energy.getE_init(), i1+noLpShift, i2+noLpShift);
					// current best total energy value (covers to far E_init only)
					curCellEtotal = energy.getE(i1,curCell->j1,i2,curCell->j2,curCell->val);
				}

				// no base case with seed so far
				curEseedtotal = E_INF;
				curCellSeedEtotal = E_INF;

				///////////////////////////////////////////////////////////////////
				// check if seed is starting here
				///////////////////////////////////////////////////////////////////

				// check if seed is possible for this left boundary
				if ( seedHandler.isSeedBound(i1,i2) ) {
					// get right extension
					size_t sj1 = i1 + seedHandler.getSeedLength1(i1,i2)-1;
					size_t sj2 = i2 + seedHandler.getSeedLength2(i1,i2)-1;
					E_type seedE = seedHandler.getSeedE(i1,i2);
					if (sj1 < hybridE.size1() && sj2 < hybridE.size2()) {

						// get energy of seed only explicitly
						curE = seedE + energy.getE_init();
						// check if this combination yields better energy
						curEseedtotal = energy.getE(i1,sj1,i2,sj2,curE);
						// update Zall for seed only (if otherwise stacking enforced)
						updateZall( i1,sj1,i2,sj2, curEseedtotal, false );

						// for noLP : check for explicit interior loop after seed
						// assumption: seed fulfills noLP
						if (outConstraint.noLP) {
							// update optimum
							if ( curEseedtotal < curCellSeedEtotal )
							{
								// update current best to seed-only
								*curCellSeed = BestInteractionE( curE, sj1, sj2 );
								// store total energy
								curCellSeedEtotal = curEseedtotal;
							}
							// iterate over all loop sizes w1 (seq1) and w2 (seq2) (minus 1)
							for (w1=1; w1-1 <= energy.getMaxInternalLoopSize1() && sj1+w1<hybridE.size1(); w1++) {
							for (w2=1; w2-1 <= energy.getMaxInternalLoopSize2() && sj2+w2<hybridE.size2(); w2++) {

								// ensure there is at least one unpaired position in loop right after the seed
								if (w1+w2==2) { continue; }

								// direct cell access (const)
								rightExt = &(hybridE(sj1+w1,sj2+w2));

								// check if right side can pair
								// check if interaction length is within boundary
								if ( E_isNotINF(rightExt->val)
									&& (rightExt->j1 +1 -i1) <= energy.getAccessibility1().getMaxLength()
									&& (rightExt->j2 +1 -i2) <= energy.getAccessibility2().getMaxLength() )
								{
									// compute energy for this loop sizes
									curEloop = energy.getE_interLeft(sj1,sj1+w1,sj2,sj2+w2);
									curE = seedE + curEloop + rightExt->val;
									// check if this combination yields better energy
									curEseedtotal = energy.getE(i1,rightExt->j1,i2,rightExt->j2,curE);
									// update Zall
									updateZall( i1,rightExt->j1,i2,rightExt->j2, curEseedtotal, false );
									if ( curEseedtotal < curCellSeedEtotal )
									{
										// update current best for this left boundary
										// copy right boundary
										*curCellSeed = *rightExt;
										// set new energy
										curCellSeed->val = curE;
										// store total energy to avoid recomputation
										curCellSeedEtotal = curEseedtotal;
									}
								}

							} } // w1 w2
						}  // noLP

						// check direct true-right extension of seed
						rightExt = &(hybridE(sj1,sj2));
						if (E_isNotINF(rightExt->val)) {
							// get energy of seed interaction with best right extension
							curE = seedE + rightExt->val;
							// check if this combination yields better energy
							curEseedtotal = energy.getE(i1,rightExt->j1,i2,rightExt->j2,curE);
							// update Zall (avoid duplicated consideration of seed-only)
							if (sj1 < rightExt->j1) {
								updateZall( i1,rightExt->j1,i2,rightExt->j2, curEseedtotal, false );
							}
							// update optimum
							if ( curEseedtotal < curCellSeedEtotal )
							{
								// update current best for this left boundary
								// copy right boundary
								*curCellSeed = *rightExt;
								// set new energy
								curCellSeed->val = curE;
								// store total energy
								curCellSeedEtotal = curEseedtotal;
							}
						} // direct right extension

					}
				}

				// if !noLP or stacking bp possible
				if (E_isNotINF(iStackE))  {
					// iterate over all loop sizes w1 (seq1) and w2 (seq2) (minus 1)
					for (w1=1; w1-1 <= energy.getMaxInternalLoopSize1() && i1+noLpShift+w1<hybridE.size1(); w1++) {
					for (w2=1; w2-1 <= energy.getMaxInternalLoopSize2() && i2+noLpShift+w2<hybridE.size2(); w2++) {

						// reset loop energy
						curEloop = E_INF;

						//////////////////////////////////////////////////////////
						// update hybridE without seed constraint
						//////////////////////////////////////////////////////////

						// direct cell access (const)
						rightExt = &(hybridE(i1+noLpShift+w1,i2+noLpShift+w2));

						// check if right side can pair
						// check if interaction length is within boundary
						if ( E_isNotINF(rightExt->val)
							&& (rightExt->j1 +1 -i1) <= energy.getAccessibility1().getMaxLength()
							&& (rightExt->j2 +1 -i2) <= energy.getAccessibility2().getMaxLength() )
						{
							// compute energy for this loop sizes
							curEloop = energy.getE_interLeft(i1+noLpShift,i1+noLpShift+w1,i2+noLpShift,i2+noLpShift+w2);
							curE = iStackE + curEloop + rightExt->val;
							// check if this combination yields better energy
							curEtotal = energy.getE(i1,rightExt->j1,i2,rightExt->j2,curE);
							if ( curEtotal < curCellEtotal )
							{
								// update current best for this left boundary
								// copy right boundary
								*curCell = *rightExt;
								// set new energy
								curCell->val = curE;
								// store total energy to avoid recomputation
								curCellEtotal = curEtotal;
							}
						}

						//////////////////////////////////////////////////////////
						// update hybridE_seed including seed constraint
						//////////////////////////////////////////////////////////

						// direct cell access to right side end of loop (seed has to be to the right of it)
						rightExt = &(hybridE_seed(i1+noLpShift+w1,i2+noLpShift+w2));
						// check if right side of loop can pair
						// check if interaction length is within boundary
						if (E_isNotINF(rightExt->val)
							&& (rightExt->j1 +1 -i1) <= energy.getAccessibility1().getMaxLength()
							&& (rightExt->j2 +1 -i2) <= energy.getAccessibility2().getMaxLength() )
						{
							// compute loop energy if not already known
							if (E_isINF(curEloop)) {
								curEloop = energy.getE_interLeft(i1+noLpShift,i1+noLpShift+w1,i2+noLpShift,i2+noLpShift+w2);
							}
							// compute energy for this loop sizes
							curE = iStackE + curEloop + rightExt->val;
							// check if this combination yields better energy
							curEseedtotal = energy.getE(i1,rightExt->j1,i2,rightExt->j2,curE);
							if ( curEseedtotal < curCellSeedEtotal )
							{
								// update current best for this left boundary
								// copy right boundary
								*curCellSeed = *rightExt;
								// set new energy
								curCellSeed->val = curE;
								// store overall energy
								curCellSeedEtotal = curEseedtotal;
							}
							// avoid Z update; otherwise double-counting of interactions
							//   --> that way, underestimation of Z
						}

					} // w2
					} // w1
				}

			// update mfe
			updateOptima(i1,curCellSeed->j1,i2,curCellSeed->j2,curCellSeedEtotal,false, false);

			} // valid base pair

		} // i2
	} // i1

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dHeuristicSeed::
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
		throw std::runtime_error("PredictorMfe2dHeuristicSeed::traceBack() : given interaction does not contain boundaries only");
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
		throw std::runtime_error("PredictorMfe2dHeuristicSeed::traceBack() : given interaction not valid");
	}
#endif

	// ensure sorting
	interaction.sort();
	// get indices in hybridE for boundary base pairs
	size_t	i1 = energy.getIndex1(interaction.basePairs.at(0)),
			i2 = energy.getIndex2(interaction.basePairs.at(0));
	const size_t j1 = energy.getIndex1(interaction.basePairs.at(1));
	const size_t j2 = energy.getIndex2(interaction.basePairs.at(1));

	// determine whether or not lonely base pairs are allowed or if we have to
	// ensure a stacking to the right of the left boundary (i1,i2)
	const size_t noLpShift = outConstraint.noLP ? 1 : 0;
	E_type iStackE = E_type(0);

	// the currently traced value for i1-j1, i2-j2
	E_type curE = hybridE_seed(i1,i2).val;
	assert( hybridE_seed(i1,i2).j1 == j1 );
	assert( hybridE_seed(i1,i2).j2 == j2 );
	assert( i1 <= j1 );
	assert( i2 <= j2 );
	assert( j1 < hybridE_seed.size1() );
	assert( j2 < hybridE_seed.size2() );

	// trace back
	// temp variables
	size_t k1,k2;
	// do until only right boundary is left over
	while( (j1-i1) > 1 ) {

		curE = hybridE_seed(i1,i2).val;

		if (noLpShift != 0) {
			// get stacking term to avoid recomputation
			iStackE = energy.getE_interLeft(i1,i1+noLpShift,i2,i2+noLpShift);
		}
		
		const BestInteractionE * curCell = NULL;
		bool traceNotFound = true;
		// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
		for (k1=std::min(j1,i1+energy.getMaxInternalLoopSize1()+1+noLpShift); traceNotFound && k1>i1+noLpShift; k1--) {
		for (k2=std::min(j2,i2+energy.getMaxInternalLoopSize2()+1+noLpShift); traceNotFound && k2>i2+noLpShift; k2--) {
			// temp access to current cell
			curCell = &(hybridE_seed(k1,k2));
			// check if right boundary is equal (part of the heuristic)
			if ( curCell->j1 == j1 && curCell->j2 == j2 &&
					// and energy is the source of curE
					E_equal( curE, (iStackE + energy.getE_interLeft(i1+noLpShift,k1,i2+noLpShift,k2) + curCell->val ) ) )
			{
				// stop searching
				traceNotFound = false;
				if ( noLpShift != 0 ) {
					interaction.basePairs.push_back( energy.getBasePair(i1+noLpShift,i2+noLpShift) );
				}
				// store splitting base pair if not last one of interaction range
				if ( k1 < j1 ) {
					interaction.basePairs.push_back( energy.getBasePair(k1,k2) );
				}
				// trace right part of split
				i1=k1;
				i2=k2;
				curE = curCell->val;
			}
		}
		}
		// has to be interaction with seed on the left starting at (i1,i2)..seed..(k1,k2)..rest..(j1,j2)
		if (traceNotFound) {
			assert( seedHandler.isSeedBound(i1,i2));
			k1 = i1+seedHandler.getSeedLength1(i1,i2)-1; assert(k1<hybridE.size1());
			k2 = i2+seedHandler.getSeedLength2(i1,i2)-1; assert(k2<hybridE.size2());
			E_type seedE = seedHandler.getSeedE(i1,i2);
			// store seed information
			interaction.setSeedRange(
							energy.getBasePair(i1,i2),
							energy.getBasePair(k1,k2),
							energy.getE(i1,k1,i2,k2,seedE)+energy.getE_init());
			// traceback seed base pairs (excludes right most = (k1,k2))
			seedHandler.traceBackSeed( interaction, i1, i2 );
			// check if seed only
			if (k1 == j1 && k2 == j2 && E_equal( curE, seedE + energy.getE_init()) ) {
				i1 = k1;
				i2 = k2;
				break;
			}
			// trace remaining base pairs
			if (outConstraint.noLP) {
				for (size_t l1=std::min(j1-1,k1+energy.getMaxInternalLoopSize1()+1); traceNotFound && l1>k1; l1--) {
				for (size_t l2=std::min(j2-1,k2+energy.getMaxInternalLoopSize2()+1); traceNotFound && l2>k2; l2--) {
					// skip invalid boundaries
					curCell = &(hybridE(l1,l2));
					if (E_isINF(curCell->val)) {
						continue;
					}
					// check boundary and energy
					if ( j1 == curCell->j1 && j2 == curCell->j2
						&& E_equal( curE, (seedE + energy.getE_interLeft(k1,l1,k2,l2) + (curCell->val)) ) )
					{
						// store seend end
						interaction.basePairs.push_back( energy.getBasePair(k1,k2) );
						// update position to point after seed interaction
						i1 = l1;
						i2 = l2;
						traceNotFound = false;
					}
				}} // l1 l2
			}
			if (traceNotFound){
				curCell = &(hybridE(k1,k2));
				// sanity check
				assert( E_equal( curE, (seedE+(curCell->val)) )
						 && j1 == curCell->j1
						 && j2 == curCell->j2	);
				// update cell
				i1 = k1;
				i2 = k2;
				traceNotFound = false;
			}
			// traceback remaining right interaction via hybridE
			if (i1<j1) {
				Interaction bpsRight(*(interaction.s1), *(interaction.s2) );
				bpsRight.basePairs.push_back( energy.getBasePair(i1,i2) );
				bpsRight.basePairs.push_back( energy.getBasePair(j1,j2) );
				PredictorMfe2dHeuristic::traceBack( bpsRight );
				// copy remaining base pairs
				Interaction::PairingVec & bps = bpsRight.basePairs;
				// copy all base pairs excluding the right most
				for (size_t i=0; i+1<bps.size(); i++) {
					interaction.basePairs.push_back( bps.at(i) );
				}
			}
			// stop search since all trace back done
			i1 = j1;
			i2 = j2;
			break;
		}

		assert( !traceNotFound );
	}
#if INTARNA_IN_DEBUG_MODE
	if ( (j2-i2) > 1 ) {
		throw std::runtime_error("PredictorMfe2dHeuristicSeed::traceBack() : trace leaves ji<j2 : "+toString(i2)+"<"+toString(j2));
	}
#endif

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

	// add all seeds that are subsets of the interaction
	seedHandler.addSeeds( interaction );

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dHeuristicSeed::
getNextBest( Interaction & curBest )
{

	const E_type curBestE = curBest.energy;

	// identify cell with next best non-overlapping interaction site
	// iterate (decreasingly) over all left interaction starts
	size_t i1,i2;
	BestInteractionE * curBestCell = NULL;
	E_type curBestCellE = E_INF;
	Interaction::BasePair curBestCellStart;
	BestInteractionE * curCell = NULL;
	E_type curCellE = E_INF;
	IndexRange r1,r2;
	for (i1=hybridE_seed.size1(); i1-- > 0;) {
		// ensure interaction site start is not covered
		if (reportedInteractions.first.covers(i1)) {
			continue;
		}
		for (i2=hybridE_seed.size2(); i2-- > 0;) {
			// ensure interaction site start is not covered
			if (reportedInteractions.second.covers(i2)) {
				continue;
			}
			// direct cell access
			curCell = &(hybridE_seed(i1,i2));
			// check if left side can pair
			if (E_isINF(curCell->val))
			{
				continue;
			}
			// get overall energy of the interaction
			curCellE = energy.getE(i1,curCell->j1,i2,curCell->j2,curCell->val);
			// or energy is too low to be considered
			// or energy is higher than current best found so far
			if (curCellE < curBestE || curCellE >= curBestCellE )
			{
				continue;
			}
			// ensure site is not overlapping
			r1.from = i1;
			r1.to = curCell->j1;
			if ( reportedInteractions.first.overlaps( r1 )) {
				continue;
			}
			r2.from = i2;
			r2.to = curCell->j2;
			if ( reportedInteractions.second.overlaps( r2 )) {
				continue;
			}
			//// FOUND THE NEXT BETTER SOLUTION
			// overwrite current best found so far
			curBestCell = curCell;
			curBestCellE = curCellE;
			curBestCellStart.first = i1;
			curBestCellStart.second = i2;

		} // i2
	} // i1

	// overwrite curBest
	curBest.basePairs.resize(2);
	curBest.energy = curBestCellE;
	if (E_isNotINF(curBestCellE)) {
		curBest.basePairs[0] = energy.getBasePair( curBestCellStart.first, curBestCellStart.second );
		curBest.basePairs[1] = energy.getBasePair( curBestCell->j1, curBestCell->j2 );
	}

}

////////////////////////////////////////////////////////////////////////////

} // namespace

