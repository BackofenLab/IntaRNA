
#include "IntaRNA/PredictorMfe2dHelixHeuristicSeed.h"

#include <stdexcept>

namespace IntaRNA {

////////////////////////////////////////////////////////////////////////////

PredictorMfe2dHelixHeuristicSeed::
PredictorMfe2dHelixHeuristicSeed( const InteractionEnergy & energy
		, OutputHandler & output
		, PredictionTracker * predTracker
		, const HelixConstraint & helixConstraint
		, SeedHandler * seedHandlerInstance )

	: PredictorMfe2dHelixHeuristic(energy,output,predTracker,helixConstraint)
		, seedHandler(seedHandlerInstance)

{
	helixHandler.setSeedHandler(seedHandler);
}


////////////////////////////////////////////////////////////////////////////

PredictorMfe2dHelixHeuristicSeed::
~PredictorMfe2dHelixHeuristicSeed()
{
	// clean up
}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dHelixHeuristicSeed::
predict( const IndexRange & r1
		, const IndexRange & r2
		, const OutputConstraint & outConstraint )
{
#if INTARNA_MULITHREADING
#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) <<"predicting mfe interactions with seed based on helices heuristically in O(n^2) space and time..."; }
	// measure timing
	TIMED_FUNC_IF(timerObj,VLOG_IS_ON(9));

#if INTARNA_IN_DEBUG_MODE
	// check indices
	if (!(r1.isAscending() && r2.isAscending()) )
		throw std::runtime_error("PredictorMfe2dHelixHeuristicSeed::predict("+toString(r1)+","+toString(r2)+") is not sane");
#endif

	// set index offset
	energy.setOffset1(r1.from);
	energy.setOffset2(r2.from);
	helixHandler.setOffset1(r1.from);
	helixHandler.setOffset2(r2.from);
	seedHandler.setOffset1(r1.from);
	seedHandler.setOffset2(r2.from);

	const size_t hybridEsize1 = std::min( energy.size1()
			, (r1.to==RnaSequence::lastPos?energy.size1()-1:r1.to)-r1.from+1 );
	const size_t hybridEsize2 = std::min( energy.size2()
			, (r2.to==RnaSequence::lastPos?energy.size2()-1:r2.to)-r2.from+1 );

	// resize matrix
	hybridE.resize( hybridEsize1, hybridEsize2 );
	hybridE_seed.resize( hybridE.size1(), hybridE.size2() );

	// Fill seed / helix and helixSeed Matrices, if one is empty trigger empty interaction reporting
	if ((seedHandler.fillSeed(0, hybridEsize1-1, 0, hybridEsize2-1) == 0)
		|| (helixHandler.fillHelix( 0, hybridEsize1-1, 0, hybridEsize2-1) == 0)
		|| (helixHandler.fillHelixSeed( 0, hybridEsize1-1, 0, hybridEsize2-1) == 0)) {
		// trigger empty interaction reporting
		initOptima(outConstraint);
		reportOptima(outConstraint);
		// stop computation
		return;
	}

	// temp vars
	size_t i1,i2,h1,h2,w1,w2;

	// init hybridE matrix
	bool isValidCell = true;
	for (i1=0; i1<hybridE.size1(); i1++) {
	for (i2=0; i2<hybridE.size2(); i2++) {

		// check if positions can form interaction
		if (	energy.isAccessible1(i1)
				&& energy.isAccessible2(i2)
				&& energy.areComplementary(i1,i2) )
		{
			// set to interaction initiation with according boundary
			hybridE(i1,i2) = BestInteraction(energy.getE_init(), i1, i2);
		} else {
			// set to infinity, ie not used
			hybridE(i1,i2) = BestInteraction(E_INF, RnaSequence::lastPos, RnaSequence::lastPos);
		}
		// init seed data
		hybridE_seed(i1,i2) = BestInteraction(E_INF, RnaSequence::lastPos, RnaSequence::lastPos);

	} // i2
	} // i1

	// init mfe without seed condition
	OutputConstraint tmpOutConstraint(1, outConstraint.reportOverlap, outConstraint.maxE, outConstraint.deltaE);
	initOptima( tmpOutConstraint );

	// compute hybridization energies WITHOUT seed condition
	// sets also -energy -hybridE
	// -> no hybrid update since updateOptima overwritten
	PredictorMfe2dHelixHeuristic::fillHybridE();

	// check result of predictions without seed if any interaction possible
	// if not no seed-containing interaction is possible neither
	if (this->mfeInteractions.begin()->energy > tmpOutConstraint.maxE || E_equal(this->mfeInteractions.begin()->energy,tmpOutConstraint.maxE)) {
		// stop computation since no favorable interaction found
		reportOptima(tmpOutConstraint);
		return;
	}

	// init mfe for later updates
	initOptima( outConstraint );

	// compute entries
	// current minimal value
	E_type curE = E_INF, curEtotal = E_INF, curCellEtotal = E_INF;
	BestInteraction * curCell = NULL;
	const BestInteraction * rightExt = NULL;

	// iterate (decreasingly) over all left interaction starts
	for (i1=hybridE_seed.size1(); i1-- > 0;) {
	for (i2=hybridE_seed.size2(); i2-- > 0;) {

		// check if left side can pair
		if (E_isINF(hybridE(i1,i2).E)) {
			continue;
		}
		// direct cell access
		curCell = &(hybridE_seed(i1,i2));
		// reset temporary variables
		curEtotal = E_INF;
		curCellEtotal = E_INF;

		// check if helix containing a seed is possible for this left boundary
		if ( E_isNotINF( helixHandler.getHelixSeedE(i1,i2) ) ) {
			// helixHandlerSeed Lengths
			h1 = helixHandler.getHelixSeedLength1(i1,i2)-1; assert(i1+h1 < hybridE_seed.size1());
			h2 = helixHandler.getHelixSeedLength2(i1,i2)-1; assert(i2+h2 < hybridE_seed.size2());

			///////////////////////////////////////////////////////////////////
			// Case: Helix_seed + E_init
			///////////////////////////////////////////////////////////////////

			curE = helixHandler.getHelixSeedE(i1,i2) + energy.getE_init();
			// check if this combination yields better energy
			curEtotal = energy.getE(i1,i1+h1, i2, i2+h2, curE);
			if ( curEtotal < curCellEtotal )
			{

				// update current best for this left boundary
				// set right boundary
				curCell->j1 = i1+h1;
				curCell->j2 = i2+h2;
				// set new energy
				curCell->E = curE;
				// store total energy to avoid recomputation
				curCellEtotal = curEtotal;
			}

			///////////////////////////////////////////////////////////////////
			// Case: Helix_seed + interior loop + hybridE
			///////////////////////////////////////////////////////////////////

			for (w1=1; w1-1 <= energy.getMaxInternalLoopSize1() && i1+h1+w1<hybridE.size1(); w1++) {
			for (w2=1; w2-1 <= energy.getMaxInternalLoopSize2() && i2+h2+w2<hybridE.size2(); w2++) {

				// skip too small interior loops or stackings
				// ensure bulges and interior loops exceed allowed ones within helices
				if ( w1+w2-2 <= helixHandler.getConstraint().getMaxIL()) {
					continue;
				}

				// direct cell access (const)
				rightExt = &(hybridE(i1+h1+w1, i2+h2+w2));
				// check if right side can pair
				if (E_isINF(rightExt->E)) {
					continue;
				}
				// check if interaction length is within boundary
				if ( (rightExt->j1+1-i1) > energy.getAccessibility1().getMaxLength()
					 || (rightExt->j2+1-i2) > energy.getAccessibility2().getMaxLength() )
				{
					continue;
				}
				// compute energy for this loop sizes
				curE = helixHandler.getHelixSeedE(i1,i2) + energy.getE_interLeft(i1+h1,i1+h1+w1,i2+h2,i2+h2+w2) + rightExt->E;

				// check if this combination yields better energy
				curEtotal = energy.getE(i1,rightExt->j1,i2,rightExt->j2,curE);
				if ( !E_equal(curEtotal, curCellEtotal) && curEtotal < curCellEtotal )
				{
					// update current best for this left boundary
					// copy right boundary
					*curCell = *rightExt;
					// set new energy
					curCell->E = curE;
					// store total energy to avoid recomputation
					curCellEtotal = curEtotal;
				}

			} // w2
			} // w1
		} // helixSeed


		///////////////////////////////////////////////////////////////////
		// Case: Helix + interior loop + hybridE_seed
		///////////////////////////////////////////////////////////////////

		// check if helix containing no seed is possible for this left boundary
		if ( E_isNotINF(helixHandler.getHelixE(i1,i2) ) ) {
			// helixHandler lengths
			h1 = helixHandler.getHelixLength1(i1,i2)-1; assert(i1+h1<hybridE.size1());
			h2 = helixHandler.getHelixLength2(i1,i2)-1; assert(i2+h2<hybridE.size2());

			// iterate over all loop size w1 (seq1) and w2 (seq2)
			for (w1=1; w1-1 <= energy.getMaxInternalLoopSize1() && i1+h1+w1<hybridE_seed.size1(); w1++) {
			for (w2=1; w2-1 <= energy.getMaxInternalLoopSize2() && i2+h2+w2<hybridE_seed.size2(); w2++) {

				// skip too small interior loops or stackings
				// ensure bulges and interior loops exceed allowed ones within helices
				if ( w1+w2-2 <= helixHandler.getConstraint().getMaxIL()) {
					continue;
				}

				// direct cell access (const)
				rightExt = &(hybridE_seed(i1+h1+w1,i2+h2+w2));
				// check if right side can pair
				if (E_isINF(rightExt->E)) {
					continue;
				}
				// check if interaction length is within boundary
				if ( (rightExt->j1 +1 -i1) > energy.getAccessibility1().getMaxLength()
					 || (rightExt->j2 +1 -i2) > energy.getAccessibility2().getMaxLength() )
				{
					continue;
				}
				// compute energy for this loop sizes
				curE = helixHandler.getHelixE(i1,i2) + energy.getE_interLeft(i1+h1,i1+h1+w1,i2+h2,i2+h2+w2) + rightExt->E;

				// check if this combination yields better energy
				curEtotal = energy.getE(i1,rightExt->j1,i2,rightExt->j2,curE);
				if ( curEtotal < curCellEtotal )
				{
					// update current best for this left boundary
					// copy right boundary
					*curCell = *rightExt;
					// set new energy
					curCell->E = curE;
					// store total energy to avoid recomputation
					curCellEtotal = curEtotal;
				}

			} // w2
			} // w1

		} // helix

		// update mfe if needed (call superclass update routine)
		PredictorMfe2dHelixHeuristic::updateOptima( i1,curCell->j1, i2,curCell->j2, curCellEtotal, false );

	} // i2
	} // i1


	// report mfe interaction
	reportOptima( outConstraint );

}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dHelixHeuristicSeed::
traceBack( Interaction & interaction, const OutputConstraint & outConstraint )
{
	// check if something to trace
	if (interaction.basePairs.size() < 2) {
		return;
	}

#if INTARNA_IN_DEBUG_MODE
	// sanity checks
	if ( ! interaction.isValid() ) {
		throw std::runtime_error("PredictorMfe2dHelixHeuristicSeed::traceBack() : given interaction not valid");
	}
	if ( interaction.basePairs.size() != 2 ) {
		throw std::runtime_error("PredictorMfe2dHelixHeuristicSeed::traceBack() : given interaction does not contain boundaries only");
	}
#endif

	// check for single interaction
	if (interaction.basePairs.at(0).first == interaction.basePairs.at(1).first) {
		// delete second boundary (identical to first)
		interaction.basePairs.resize(1);
		// update done
		return;
	}

	// ensure sorting
	interaction.sort();
	// get indices in hybridE for boundary base pairs
	size_t	i1 = energy.getIndex1(interaction.basePairs.at(0)),
			i2 = energy.getIndex2(interaction.basePairs.at(0));
	const size_t j1 = energy.getIndex1(interaction.basePairs.at(1));
	const size_t j2 = energy.getIndex2(interaction.basePairs.at(1));

	// the currently traced value for i1-j1, i2-j2
	E_type curE = hybridE_seed(i1,i2).E;
	assert( hybridE_seed(i1,i2).j1 == j1 );
	assert( hybridE_seed(i1,i2).j2 == j2 );
	assert( i1 <= j1 );
	assert( i2 <= j2 );
	assert( j1 < hybridE_seed.size1() );
	assert( j2 < hybridE_seed.size2() );

	// trace back
	// temp variables
	size_t h1,h2,k1,k2;
	// do until only right boundary is left over
	while( (j1-i1) > 1 ) {
		const BestInteraction * curCell = NULL;
		bool traceNotFound = true;

		// Assure that atleast one case is possible
		assert(E_isNotINF(helixHandler.getHelixE(i1,i2)) || E_isNotINF(helixHandler.getHelixSeedE(i1,i2)));

		// helix + il + hybridE_seed
		if (E_isNotINF(helixHandler.getHelixE(i1,i2))) {
			h1 = helixHandler.getHelixLength1(i1,i2)-1; assert(i1+h1<hybridE.size1());
			h2 = helixHandler.getHelixLength2(i1,i2)-1; assert(i2+h2<hybridE.size2());

			// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
			for (size_t w1=1; traceNotFound && w1-1 <= energy.getMaxInternalLoopSize1() && i1+h1+w1<hybridE_seed.size1(); w1++) {
			for (size_t w2=1; traceNotFound && w2-1 <= energy.getMaxInternalLoopSize2() && i2+h2+w2<hybridE_seed.size2(); w2++) {

				// skip too small interior loops or stackings
				// ensure bulges and interior loops exceed allowed ones within helices
				if ( w1+w2-2 <= helixHandler.getConstraint().getMaxIL()) {
					continue;
				}

				k1 = i1+h1+w1;
				k2 = i2+h2+w2;

				// temp access to current cell
				curCell = &(hybridE_seed(k1,k2));
				// check if right boundary is equal (part of the heuristic)
				if ( curCell->j1 == j1 && curCell->j2 == j2 &&
					 // and energy is the source of curE
					 E_equal( curE, (helixHandler.getHelixE(i1,i2) + energy.getE_interLeft(i1+h1,k1,i2+h2,k2) + curCell->E ) ) )
				{
					// stop searching
					traceNotFound = false;
					// store helix base pairs
					helixHandler.traceBackHelix( interaction, i1, i2 );

					// stor last base pair of helix
					interaction.basePairs.push_back( energy.getBasePair(i1+h1,i2+h2));
					// store splitting base pair
					interaction.basePairs.push_back( energy.getBasePair(k1,k2) );

					// trace right part of split
					i1=k1;
					i2=k2;
					curE = curCell->E;
				}
			} // w1
			} // w2

		}

		// seed + il + hybridE
		if (E_isNotINF(helixHandler.getHelixSeedE(i1,i2))) {
			h1 = helixHandler.getHelixSeedLength1(i1,i2)-1; assert(i1+h1<hybridE_seed.size1());
			h2 = helixHandler.getHelixSeedLength2(i1,i2)-1; assert(i2+h2<hybridE_seed.size2());
			// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
			for (size_t w1=1; traceNotFound && w1-1 <= energy.getMaxInternalLoopSize1() && i1+h1+w1<hybridE.size1(); w1++) {
			for (size_t w2=1; traceNotFound && w2-1 <= energy.getMaxInternalLoopSize2() && i2+h2+w2<hybridE.size2(); w2++) {

				// skip too small interior loops or stackings
				// ensure bulges and interior loops exceed allowed ones within helices
				if ( w1+w2-2 <= helixHandler.getConstraint().getMaxIL()) {
					continue;
				}

				k1 = i1+h1+w1;
				k2 = i2+h2+w2;

				// temp access to current cell
				curCell = &(hybridE(k1,k2));
				// check if right boundary is equal (part of the heuristic)
				if ( curCell->j1 == j1 && curCell->j2 == j2 &&
					 // and energy is the source of curE
					 E_equal( curE, (helixHandler.getHelixSeedE(i1,i2) + energy.getE_interLeft(i1+h1,k1,i2+h2,k2) + curCell->E ) ) )
				{
					// store helix base pairs
					helixHandler.traceBackHelixSeed( interaction, i1, i2 );
					interaction.basePairs.push_back(energy.getBasePair(i1+h1, i2+h2));
					i1 = k1;
					i2 = k2;

					// traceback remaining right interaction via hybridE
					if (i1 < j1) {
						Interaction bpsRight(*(interaction.s1), *(interaction.s2));
						bpsRight.basePairs.push_back(energy.getBasePair(i1, i2));
						bpsRight.basePairs.push_back(energy.getBasePair(j1, j2));
						PredictorMfe2dHelixHeuristic::traceBack(bpsRight, outConstraint);
						// copy remaining base pairs
						Interaction::PairingVec &bps = bpsRight.basePairs;
						// copy all base pairs excluding the right most
						for (size_t i = 0; i + 1 < bps.size(); i++) {
							interaction.basePairs.push_back(bps.at(i));
						}
					}
					// stop search since all trace back done
					traceNotFound = false;
					i1=j1;
					i2=j2;
				}
			} // w1
			} // w2

			// seed + E_init()
			if ( traceNotFound && E_equal(curE, helixHandler.getHelixSeedE(i1,i2) + energy.getE_init())  ) {
				// stop searching
				traceNotFound = false;
				// traceback helix base pairs ( excluding right most = (k1,k2))
				helixHandler.traceBackHelixSeed(interaction, i1, i2);
				// trace right part of split
				i1=i1+h1;
				i2=i2+h2;
				curE=0.0;
			}
		}
		assert(!traceNotFound);

	}
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

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dHelixHeuristicSeed::
getNextBest( Interaction & curBest )
{

	const E_type curBestE = curBest.energy;

	// TODO replace index iteration with something based on ranges from reportedInteractions

	// identify cell with next best non-overlapping interaction site
	// iterate (decreasingly) over all left interaction starts
	size_t i1,i2;
	BestInteraction * curBestCell = NULL;
	E_type curBestCellE = E_INF;
	Interaction::BasePair curBestCellStart;
	BestInteraction * curCell = NULL;
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
			if (E_isINF(curCell->E))
			{
				continue;
			}
			// get overall energy of the interaction
			curCellE = energy.getE(i1,curCell->j1,i2,curCell->j2,curCell->E);
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

