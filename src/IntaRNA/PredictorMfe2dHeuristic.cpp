
#include "IntaRNA/PredictorMfe2dHeuristic.h"

#include <stdexcept>

namespace IntaRNA {

////////////////////////////////////////////////////////////////////////////

PredictorMfe2dHeuristic::
PredictorMfe2dHeuristic(
		const InteractionEnergy & energy
		, OutputHandler & output
		, PredictionTracker * predTracker )
 : PredictorMfe(energy,output,predTracker)
	, hybridE( 0,0 )
{
}


////////////////////////////////////////////////////////////////////////////

PredictorMfe2dHeuristic::
~PredictorMfe2dHeuristic()
{
	// clean up
}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dHeuristic::
predict( const IndexRange & r1
		, const IndexRange & r2
		)
{
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) <<"predicting mfe interactions heuristically in O(n^2) space and time..."; }
	// measure timing
	TIMED_FUNC_IF(timerObj,VLOG_IS_ON(9));

#if INTARNA_IN_DEBUG_MODE
	// check indices
	if (!(r1.isAscending() && r2.isAscending()) )
		throw std::runtime_error("PredictorMfe2dHeuristic::predict("+toString(r1)+","+toString(r2)+") is not sane");
#endif


	// set index offset
	energy.setOffset1(r1.from);
	energy.setOffset2(r2.from);

	// resize matrix
	hybridE.resize( std::min( energy.size1()
						, (r1.to==RnaSequence::lastPos?energy.size1()-1:r1.to)-r1.from+1 )
				, std::min( energy.size2()
						, (r2.to==RnaSequence::lastPos?energy.size2()-1:r2.to)-r2.from+1 ) );

	// init mfe for later updates
	initOptima();

	// compute table and update mfeInteraction
	fillHybridE();

	// trace back and output handler update
	reportOptima();

}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dHeuristic::
fillHybridE()
{
	// temporary access
	const OutputConstraint & outConstraint = output.getOutputConstraint();
	// compute entries
	// current minimal value
	E_type curE = E_INF, curEtotal = E_INF, curCellEtotal = E_INF;
	size_t i1,i2,w1,w2;

	// determine whether or not lonely base pairs are allowed or if we have to
	// ensure a stacking to the right of the left boundary (i1,i2)
	const size_t noLpShift = outConstraint.noLP ? 1 : 0;
	E_type iStackE = E_type(0);

	BestInteractionE * curCell = NULL;
	const BestInteractionE * rightExt = NULL;
	// iterate (decreasingly) over all left interaction starts
	for (i1=hybridE.size1(); i1-- > 0;) {
		for (i2=hybridE.size2(); i2-- > 0;) {
			// direct cell access
			curCell = &(hybridE(i1,i2));

			// init as invalid boundary
			*curCell = BestInteractionE(E_INF, RnaSequence::lastPos, RnaSequence::lastPos);

			// check if positions can form interaction
			if (	energy.isAccessible1(i1)
					&& energy.isAccessible2(i2)
					&& energy.areComplementary(i1,i2) )
			{
				// no LP allowed
				if (outConstraint.noLP) {
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
						continue;
					}
				}

				// set to interaction initiation with according boundary
				// if valid right boundary
				if (!outConstraint.noGUend || !energy.isGU(i1+noLpShift,i2+noLpShift))
				{
					*curCell = BestInteractionE(iStackE + energy.getE_init(), i1+noLpShift, i2+noLpShift);
					// current best total energy value (covers to far E_init only)
					curCellEtotal = energy.getE(i1,curCell->j1,i2,curCell->j2,curCell->val);
					// update Zall
					updateZall( i1,curCell->j1,i2,curCell->j2, curCellEtotal, false );
				}

				// iterate over all loop sizes w1 (seq1) and w2 (seq2) (minus 1)
				for (w1=1; w1-1 <= energy.getMaxInternalLoopSize1() && i1+w1+noLpShift<hybridE.size1(); w1++) {
				for (w2=1; w2-1 <= energy.getMaxInternalLoopSize2() && i2+w2+noLpShift<hybridE.size2(); w2++) {
					// direct cell access (const)
					rightExt = &(hybridE(i1+noLpShift+w1,i2+noLpShift+w2));
					// check if right side can pair
					if (E_isINF(rightExt->val)) {
						continue;
					}
					// check if interaction length is within boundary
					if ( (rightExt->j1 +1 -i1) > energy.getAccessibility1().getMaxLength()
						|| (rightExt->j2 +1 -i2) > energy.getAccessibility2().getMaxLength() )
					{
						continue;
					}
					// compute energy for this loop sizes
					curE = iStackE + energy.getE_interLeft(i1+noLpShift,i1+noLpShift+w1,i2+noLpShift,i2+noLpShift+w2) + rightExt->val;
					// check if this combination yields better energy
					curEtotal = energy.getE(i1,rightExt->j1,i2,rightExt->j2,curE);
					// update Zall
					updateZall( i1,rightExt->j1,i2,rightExt->j2, curEtotal, false );
					// update best extension
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

				} // w2
				} // w1

				// update mfe if needed
				updateOptima( i1,curCell->j1, i2,curCell->j2, curCellEtotal, false, false );

			} // valid base pair

		} // i2
	} // i1

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dHeuristic::
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
		throw std::runtime_error("PredictorMfe2dHeuristic::traceBack() : given interaction does not contain boundaries only");
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
		throw std::runtime_error("PredictorMfe2dHeuristic::traceBack() : given interaction not valid : "+toString(interaction));
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
	E_type curE = hybridE(i1,i2).val;
	assert( hybridE(i1,i2).j1 == j1 );
	assert( hybridE(i1,i2).j2 == j2 );
	assert( i1 <= j1 );
	assert( i2 <= j2 );
	assert( j1 < hybridE.size1() );
	assert( j2 < hybridE.size2() );

	// trace back
	// temp variables
	size_t k1,k2;
	// check all possible splits of the interval (i1,i2)-(j1,j2)
	// only reasonable, if there is an enclosed position k1 between i1-j1
	while( (j1-i1) > 1 ) {

		curE = hybridE(i1,i2).val;

		if (outConstraint.noLP) {

			// get stacking term to avoid recomputation
			iStackE = energy.getE_interLeft(i1,i1+1,i2,i2+1);
			assert(E_isNotINF(iStackE));

			// check just stacking
			if ( hybridE(i1+1,i2+1).j1 == j1
					&& hybridE(i1+1,i2+1).j2 == j2
					&& E_equal( curE, ( iStackE + hybridE(i1+1,i2+1).val ) ) )
			{
				// trace right part of split
				i1++;
				i2++;
				curE = hybridE(i1,i2).val;
				// store splitting base pair if not last one of interaction range
				if (i1 != j1) {
					interaction.basePairs.push_back( energy.getBasePair(i1,i2) );
				}
				// trace update i1,i2
				continue;
			}
		}

		const BestInteractionE * curCell = NULL;
		bool traceNotFound = true;
		// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
		for (k1=std::min(j1,i1+energy.getMaxInternalLoopSize1()+1+noLpShift); traceNotFound && k1>i1+noLpShift; k1--) {
		for (k2=std::min(j2,i2+energy.getMaxInternalLoopSize2()+1+noLpShift); traceNotFound && k2>i2+noLpShift; k2--) {
			// temp access to current cell
			curCell = &(hybridE(k1,k2));
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
		assert( !traceNotFound );
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
PredictorMfe2dHeuristic::
getNextBest( Interaction & curBest )
{

	// get original
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
	for (i1=hybridE.size1(); i1-- > 0;) {
		// ensure interaction site start is not covered
		if (reportedInteractions.first.covers(i1)) {
			continue;
		}
		for (i2=hybridE.size2(); i2-- > 0;) {
			// ensure interaction site start is not covered
			if (reportedInteractions.second.covers(i2)) {
				continue;
			}
			// direct cell access
			curCell = &(hybridE(i1,i2));
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
	curBest.energy = curBestCellE;
	curBest.basePairs.resize(2);
	if (E_isNotINF(curBestCellE)) {
		curBest.basePairs[0] = energy.getBasePair( curBestCellStart.first, curBestCellStart.second );
		curBest.basePairs[1] = energy.getBasePair( curBestCell->j1, curBestCell->j2 );
	}
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dHeuristic::
updateMfe4leftEnd(const size_t i1, const size_t j1
					, const size_t i2, const size_t j2
					, const Interaction & curInteraction )
{
	// do nothing since getNextBest() is based on local data structure
}

////////////////////////////////////////////////////////////////////////////

} // namespace

