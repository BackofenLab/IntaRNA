
#include "PredictorMfe2dHeuristic.h"

#include <stdexcept>


////////////////////////////////////////////////////////////////////////////

PredictorMfe2dHeuristic::
PredictorMfe2dHeuristic( const InteractionEnergy & energy, OutputHandler & output )
 : PredictorMfe(energy,output)
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
		, const size_t reportMax
		, const bool reportNonOverlapping
		)
{

	VLOG(2) <<"predicting mfe interactions heuristically in O(n^2) space and time...";
	// measure timing
	TIMED_FUNC_IF(timerObj,VLOG_IS_ON(9));

	// suboptimal setup check
	if (reportMax>1 && reportNonOverlapping) {
		NOTIMPLEMENTED("non-overlapping not implemented in this mode");
	}

#if IN_DEBUG_MODE
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

	// temp vars
	size_t i1,i2,w1,w2;

	// init matrix
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
//		// alternative init obsolete due to default constructor used in matrix cell creation
//		} else {
//			// set to infinity, ie not used
//			hybridE(i1,i2) = BestInteraction(E_INF, RnaSequence::lastPos, RnaSequence::lastPos);
		}

	} // i2
	} // i1

	// init mfe for later updates
	initOptima( reportMax, reportNonOverlapping );

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
	// compute entries
	// current minimal value
	E_type curE = E_INF;
	size_t i1,i2,w1,w2;
	BestInteraction * curCell = NULL;
	const BestInteraction * rightExt = NULL;
	// iterate (decreasingly) over all left interaction starts
	for (i1=hybridE.size1(); i1-- > 0;) {
	for (i2=hybridE.size2(); i2-- > 0;) {
		// direct cell access
		curCell = &(hybridE(i1,i2));
		// check if left side can pair
		if (E_isINF(hybridE(i1,i2).E)) {
			continue;
		}

		// TODO PARALLELIZE THIS DOUBLE LOOP ?!
		// iterate over all loop sizes w1 (seq1) and w2 (seq2) (minus 1)
		for (w1=1; w1-1 <= energy.getMaxInternalLoopSize1() && i1+w1<hybridE.size1(); w1++) {
		for (w2=1; w2-1 <= energy.getMaxInternalLoopSize2() && i2+w2<hybridE.size2(); w2++) {
			// direct cell access (const)
			rightExt = &(hybridE(i1+w1,i2+w2));
			// check if right side can pair
			if (E_isINF(rightExt->E)) {
				continue;
			}
			// compute energy for this loop sizes
			curE = energy.getE_interLeft(i1,i1+w1,i2,i2+w2) + rightExt->E;
			// check if this combination yields better energy
			if (energy.getE(i1,rightExt->j1,i2,rightExt->j2,curE)
					< energy.getE(i1,curCell->j1,i2,curCell->j2,curCell->E))
			{
				// update current best for this left boundary
				// copy right boundary
				*curCell = *rightExt;
				// set new energy
				curCell->E = curE;
			}

		} // w2
		} // w1

		// update mfe if needed
		updateOptima( i1,curCell->j1, i2,curCell->j2, curCell->E );

	} // i2
	} // i1

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dHeuristic::
traceBack( Interaction & interaction )
{
	// check if something to trace
	if (interaction.basePairs.size() < 2) {
		return;
	}

#if IN_DEBUG_MODE
	// sanity checks
	if ( ! interaction.isValid() ) {
		throw std::runtime_error("PredictorMfe2dHeuristic::traceBack() : given interaction not valid");
	}
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

	// ensure sorting
	interaction.sort();
	// get indices in hybridE for boundary base pairs
	size_t	i1 = energy.getIndex1(interaction.basePairs.at(0)),
			i2 = energy.getIndex2(interaction.basePairs.at(0));
	const size_t j1 = energy.getIndex1(interaction.basePairs.at(1));
	const size_t j2 = energy.getIndex2(interaction.basePairs.at(1));


	// the currently traced value for i1-j1, i2-j2
	E_type curE = hybridE(i1,i2).E;
	assert( hybridE(i1,i2).j1 == j1 );
	assert( hybridE(i1,i2).j2 == j2 );
	assert( i1 <= j1 );
	assert( i2 <= j2 );
	assert( j1 < hybridE.size1() );
	assert( j2 < hybridE.size2() );

	// trace back
	// temp variables
	size_t k1,k2;
	// do until only right boundary is left over
	while( (j1-i1) > 1 && (j2-i2) > 1 ) {
		const BestInteraction * curCell = NULL;
		bool traceNotFound = true;
		// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
		for (k1=std::min(j1,i1+energy.getMaxInternalLoopSize1()+1); traceNotFound && k1>i1; k1--) {
		for (k2=std::min(j2,i2+energy.getMaxInternalLoopSize2()+1); traceNotFound && k2>i2; k2--) {
			// temp access to current cell
			curCell = &(hybridE(k1,k2));
			// check if right boundary is equal (part of the heuristic)
			if ( curCell->j1 == j1 && curCell->j2 == j2 &&
					// and energy is the source of curE
					E_equal( curE, (energy.getE_interLeft(i1,k1,i2,k2) + curCell->E ) ) )
			{
				// stop searching
				traceNotFound = false;
				// store splitting base pair
				interaction.basePairs.push_back( energy.getBasePair(k1,k2) );
				// trace right part of split
				i1=k1;
				i2=k2;
				curE = curCell->E;
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

