
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
		, const IndexRange & r2
		, const OutputConstraint & outConstraint )
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
		initOptima(outConstraint);
		reportOptima(outConstraint);
		// stop computation
		return;
	}

	// resize matrix
	hybridE.resize( hybridEsize1, hybridEsize2 );
	hybridE_seed.resize( hybridE.size1(), hybridE.size2() );

	// temp vars
	size_t i1,i2,w1,w2;

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
	// -> no tracker update since updateOptima overwritten
	PredictorMfe2dHeuristic::fillHybridE();

	// check if any interaction possible
	// if not no seed-containing interaction is possible neither
	if (!(this->mfeInteractions.begin()->energy < tmpOutConstraint.maxE)) {
		// stop computation since no favorable interaction found
		reportOptima(tmpOutConstraint);
		return;
	}

	// reinit mfe for later updates with final information
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
		// current best total energy value
		// NOTE: by setting to E_INF instead of getE(curCell->E) we ignore the
		// single intermolecular bp case to avoid the effect of extremely low
		// EDs of single positions
		curCellEtotal = E_INF;

		///////////////////////////////////////////////////////////////////
		// check all extensions of interactions CONTAINING a seed already
		///////////////////////////////////////////////////////////////////

		// TODO PARALLELIZE THIS DOUBLE LOOP ?!
		// iterate over all loop sizes w1 (seq1) and w2 (seq2)
		for (w1=1; w1-1 <= energy.getMaxInternalLoopSize1() && i1+w1<hybridE_seed.size1(); w1++) {
		for (w2=1; w2-1 <= energy.getMaxInternalLoopSize2() && i2+w2<hybridE_seed.size2(); w2++) {
			// direct cell access to right side end of loop (seed has to be to the right of it)
			rightExt = &(hybridE_seed(i1+w1,i2+w2));
			// check if right side of loop can pair
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
			curE = energy.getE_interLeft(i1,i1+w1,i2,i2+w2) + rightExt->E;
			// check if this combination yields better energy
			curEtotal = energy.getE(i1,rightExt->j1,i2,rightExt->j2,curE);
			if ( curEtotal < curCellEtotal )
			{
				// update current best for this left boundary
				// copy right boundary
				*curCell = *rightExt;
				// set new energy
				curCell->E = curE;
				// store overall energy
				curCellEtotal = curEtotal;
			}
		} // w2
		} // w1

		///////////////////////////////////////////////////////////////////
		// check if seed is starting here
		///////////////////////////////////////////////////////////////////

		// check if seed is possible for this left boundary
		if ( E_isNotINF( seedHandler.getSeedE(i1,i2) ) ) {
			// get right extension
			w1 = seedHandler.getSeedLength1(i1,i2)-1; assert(i1+w1 < hybridE.size1());
			w2 = seedHandler.getSeedLength2(i1,i2)-1; assert(i2+w2 < hybridE.size2());
			rightExt = &(hybridE(i1+w1,i2+w2));
			// get energy of seed interaction with best right extension
			curE = seedHandler.getSeedE(i1,i2) + rightExt->E;
			// check if this combination yields better energy
			curEtotal = energy.getE(i1,rightExt->j1,i2,rightExt->j2,curE);
			if ( curEtotal < curCellEtotal )
			{
				// update current best for this left boundary
				// copy right boundary
				*curCell = *rightExt;
				// set new energy
				curCell->E = curE;
				// store total energy
				curCellEtotal = curEtotal;
			}
		}

		// update mfe if needed (call superclass update routine)
		PredictorMfe2dHeuristic::updateOptima( i1,curCell->j1, i2,curCell->j2, curCellEtotal, false );

	} // i2
	} // i1


	// report mfe interaction
	reportOptima( outConstraint );

}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dHeuristicSeed::
traceBack( Interaction & interaction, const OutputConstraint & outConstraint  )
{
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
	size_t k1,k2;
	// do until only right boundary is left over
	while( (j1-i1) > 1 ) {
		const BestInteraction * curCell = NULL;
		bool traceNotFound = true;
		// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
		for (k1=std::min(j1,i1+energy.getMaxInternalLoopSize1()+1); traceNotFound && k1>i1; k1--) {
		for (k2=std::min(j2,i2+energy.getMaxInternalLoopSize2()+1); traceNotFound && k2>i2; k2--) {
			// temp access to current cell
			curCell = &(hybridE_seed(k1,k2));
			// check if right boundary is equal (part of the heuristic)
			if ( curCell->j1 == j1 && curCell->j2 == j2 &&
					// and energy is the source of curE
					E_equal( curE, (energy.getE_interLeft(i1,k1,i2,k2) + curCell->E ) ) )
			{
				// stop searching
				traceNotFound = false;
				// store splitting base pair
				if (k1 < j1) {
					interaction.basePairs.push_back( energy.getBasePair(k1,k2) );
				}
				// trace right part of split
				i1=k1;
				i2=k2;
				curE = curCell->E;
			}
		}
		}
		// has to be interaction with seed on the left starting at (i1,i2)..seed..(k1,k2)..rest..(j1,j2)
		if (traceNotFound) {
			assert(E_isNotINF(seedHandler.getSeedE(i1,i2)));
			k1 = i1+seedHandler.getSeedLength1(i1,i2)-1; assert(k1<hybridE.size1());
			k2 = i2+seedHandler.getSeedLength2(i1,i2)-1; assert(k2<hybridE.size2());
			curCell = &(hybridE(k1,k2)); assert( E_equal( curE, (seedHandler.getSeedE(i1,i2)+curCell->E) ));
			// store seed information
			interaction.setSeedRange(
							energy.getBasePair(i1,i2),
							energy.getBasePair(k1,k2),
							energy.getE(i1,k1,i2,k2,seedHandler.getSeedE(i1,i2))+energy.getE_init());
			// traceback seed base pairs (excludes right most = (k1,k2))
			seedHandler.traceBackSeed( interaction, i1, i2 );
			// update position to point after seed interaction
			i1 = k1;
			i2 = k2;
			// traceback remaining right interaction via hybridE
			if (i1<j1) {
				Interaction bpsRight(*(interaction.s1), *(interaction.s2) );
				bpsRight.basePairs.push_back( energy.getBasePair(i1,i2) );
				bpsRight.basePairs.push_back( energy.getBasePair(j1,j2) );
				PredictorMfe2dHeuristic::traceBack( bpsRight, outConstraint );
				// copy remaining base pairs
				Interaction::PairingVec & bps = bpsRight.basePairs;
				// copy all base pairs excluding the right most
				for (size_t i=0; i+1<bps.size(); i++) {
					interaction.basePairs.push_back( bps.at(i) );
				}
			}
			traceNotFound = false;
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

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dHeuristicSeed::
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

