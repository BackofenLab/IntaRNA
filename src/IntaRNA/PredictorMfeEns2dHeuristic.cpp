
#include "IntaRNA/PredictorMfeEns2dHeuristic.h"

#include <stdexcept>

namespace IntaRNA {

////////////////////////////////////////////////////////////////////////////

PredictorMfeEns2dHeuristic::
PredictorMfeEns2dHeuristic(
		const InteractionEnergy & energy
		, OutputHandler & output
		, PredictionTracker * predTracker )
 : PredictorMfeEns2d(energy,output,predTracker)
{
	checkKeyBoundaries(std::max(energy.getAccessibility1().getMaxLength(), energy.getAccessibility2().getMaxLength()));
}


////////////////////////////////////////////////////////////////////////////

PredictorMfeEns2dHeuristic::
~PredictorMfeEns2dHeuristic()
{
	// clean up
}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2dHeuristic::
predict( const IndexRange & r1
		, const IndexRange & r2
		, const OutputConstraint & outConstraint
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
		throw std::runtime_error("PredictorMfeEns2dHeuristic::predict("+toString(r1)+","+toString(r2)+") is not sane");
#endif


	// set index offset
	energy.setOffset1(r1.from);
	energy.setOffset2(r2.from);

	// resize matrix
	hybridZ.resize( std::min( energy.size1()
						, (r1.to==RnaSequence::lastPos?energy.size1()-1:r1.to)-r1.from+1 )
				, std::min( energy.size2()
						, (r2.to==RnaSequence::lastPos?energy.size2()-1:r2.to)-r2.from+1 ) );

	// init mfe for later updates
	initOptima( outConstraint );
	// initialize overall partition function for updates
	initZ( outConstraint );

	// compute table and update mfeInteraction
	fillHybridZ( outConstraint );

	// update ensemble mfe
	for (std::unordered_map<size_t, ZPartition >::const_iterator it = Z_partitions.begin(); it != Z_partitions.end(); ++it)
	{
		// if partition function is > 0
		if (Z_isNotINF(it->second.partZ) && it->second.partZ > 0) {
			//LOG(DEBUG) << "partZ: " << it->second.partZ;
			PredictorMfe::updateOptima( it->second.i1, it->second.j1, it->second.i2, it->second.j2, energy.getE(it->second.partZ), true );
		}
	}

	LOG(DEBUG) <<"Overall Z = "<< getOverallZ();
	LOG(DEBUG) <<"Overall E = "<<E_2_Ekcal(energy.getE(getOverallZ()));

	// trace back and output handler update
	reportOptima( outConstraint );

}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2dHeuristic::
fillHybridZ( const OutputConstraint & outConstraint )
{
	// compute entries
	// current minimal value
	Z_type curE = Z_INF;
	E_type curEtotal = E_INF, curCellEtotal = E_INF;
	size_t i1,i2,w1,w2;

	// determine whether or not lonely base pairs are allowed or if we have to
	// ensure a stacking to the right of the left boundary (i1,i2)
	const size_t noLpShift = outConstraint.noLP ? 1 : 0;
	Z_type iStackE = Z_type(0);

	BestInteraction * curCell = NULL;
	const BestInteraction * rightExt = NULL;
	// iterate (decreasingly) over all left interaction starts
	for (i1=hybridZ.size1(); i1-- > 0;) {
		for (i2=hybridZ.size2(); i2-- > 0;) {
			// direct cell access
			curCell = &(hybridZ(i1,i2));

			// init as invalid boundary
			*curCell = BestInteraction(0.0, RnaSequence::lastPos, RnaSequence::lastPos);

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
						iStackE = energy.getBoltzmannWeight(energy.getE_interLeft(i1,i1+noLpShift,i2,i2+noLpShift));
					} else {
						// skip further processing, since no stacking possible
						continue;
					}
				}

				// set to interaction initiation with according boundary
				*curCell = BestInteraction(iStackE + energy.getBoltzmannWeight(energy.getE_init()), i1+noLpShift, i2+noLpShift);

				// current best total energy value (covers to far E_init only)
				curCellEtotal = energy.getE(curCell->Z);

				// iterate over all loop sizes w1 (seq1) and w2 (seq2) (minus 1)
				for (w1=1; w1-1 <= energy.getMaxInternalLoopSize1() && i1+w1+noLpShift<hybridZ.size1(); w1++) {
				for (w2=1; w2-1 <= energy.getMaxInternalLoopSize2() && i2+w2+noLpShift<hybridZ.size2(); w2++) {
					// direct cell access (const)
					rightExt = &(hybridZ(i1+noLpShift+w1,i2+noLpShift+w2));
					// check if right side can pair
					if (Z_equal(rightExt->Z, 0.0)) {
						continue;
					}
					// check if interaction length is within boundary
					if ( (rightExt->j1 +1 -i1) > energy.getAccessibility1().getMaxLength()
						|| (rightExt->j2 +1 -i2) > energy.getAccessibility2().getMaxLength() )
					{
						continue;
					}

					// compute energy for this loop sizes
					curE = iStackE + energy.getBoltzmannWeight(energy.getE_interLeft(i1+noLpShift,i1+noLpShift+w1,i2+noLpShift,i2+noLpShift+w2)) * rightExt->Z;
					// check if this combination yields better energy
					curEtotal = energy.getE(curE);

					if ( curEtotal < curCellEtotal )
					{
						// update current best for this left boundary
						// copy right boundary
						*curCell = *rightExt;
						// set new energy
						curCell->Z += curE;
						// store total energy to avoid recomputation
						curCellEtotal = curEtotal;
					}

				} // w2
				} // w1

				// update mfe if needed
				updateZ( i1,curCell->j1, i2,curCell->j2, curCell->Z, true );

			} // valid base pair

		} // i2
	} // i1

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2dHeuristic::
getNextBest( Interaction & curBest )
{

	// get original
	const Z_type curBestE = curBest.energy;

	// TODO replace index iteration with something based on ranges from reportedInteractions

	// identify cell with next best non-overlapping interaction site
	// iterate (decreasingly) over all left interaction starts
	size_t i1,i2;
	BestInteraction * curBestCell = NULL;
	Z_type curBestCellE = Z_INF;
	Interaction::BasePair curBestCellStart;
	BestInteraction * curCell = NULL;
	Z_type curCellE = Z_INF;
	IndexRange r1,r2;
	for (i1=hybridZ.size1(); i1-- > 0;) {
		// ensure interaction site start is not covered
		if (reportedInteractions.first.covers(i1)) {
			continue;
		}
		for (i2=hybridZ.size2(); i2-- > 0;) {
			// ensure interaction site start is not covered
			if (reportedInteractions.second.covers(i2)) {
				continue;
			}
			// direct cell access
			curCell = &(hybridZ(i1,i2));
			// check if left side can pair
			if (Z_isINF(curCell->Z))
			{
				continue;
			}
			// get overall energy of the interaction
			curCellE = energy.getE(curCell->Z);
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

} // namespace
