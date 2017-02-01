
#include "PredictorMfe2d.h"

#include <stdexcept>

////////////////////////////////////////////////////////////////////////////

PredictorMfe2d::
PredictorMfe2d( const InteractionEnergy & energy, OutputHandler & output )
 : PredictorMfe(energy,output)
	, hybridE_pq( 0,0 )
	, hybridErange( energy.getAccessibility1().getSequence()
			, energy.getAccessibility2().getAccessibilityOrigin().getSequence() )
{
}


////////////////////////////////////////////////////////////////////////////

PredictorMfe2d::
~PredictorMfe2d()
{
	// clean up
}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2d::
predict( const IndexRange & r1
		, const IndexRange & r2
		, const OutputConstraint & outConstraint )
{
#ifdef INTARNA_MULITHREADING
	#pragma omp critical(intarna_logOutput)
#endif
	{ VLOG(2) <<"predicting mfe interactions in O(n^2) space..."; }
	// measure timing
	TIMED_FUNC_IF(timerObj,VLOG_IS_ON(9));

	// suboptimal setup check
	if (outConstraint.reportMax>1 && outConstraint.reportOverlap != OutputConstraint::ReportOverlap::OVERLAP_BOTH) {
		throw std::runtime_error("PredictorMfe2d : the enumeration of non-overlapping suboptimal interactions is not supported in this prediction mode");
	}

#if IN_DEBUG_MODE
	// check indices
	if (!(r1.isAscending() && r2.isAscending()) )
		throw std::runtime_error("PredictorMfe2d::predict("+toString(r1)+","+toString(r2)+") is not sane");
#endif


	// set index offset
	energy.setOffset1(r1.from);
	energy.setOffset2(r2.from);

	// resize matrix
	hybridE_pq.resize( std::min( energy.size1()
						, (r1.to==RnaSequence::lastPos?energy.size1()-1:r1.to)-r1.from+1 )
				, std::min( energy.size2()
						, (r2.to==RnaSequence::lastPos?energy.size2()-1:r2.to)-r2.from+1 ) );

	// initialize mfe interaction for updates
	initOptima( outConstraint );

	// for all right ends j1
	for (size_t j1 = hybridE_pq.size1(); j1-- > 0; ) {
		// check if j1 is accessible
		if (!energy.isAccessible1(j1))
			continue;
		// iterate over all right ends j2
		for (size_t j2 = hybridE_pq.size2(); j2-- > 0; ) {
			// check if j2 is accessible
			if (!energy.isAccessible2(j2))
				continue;
			// check if base pair (j1,j2) possible
			if (!energy.areComplementary( j1, j2 ))
				continue;

			// fill matrix and store best interaction
			fillHybridE( j1, j2 );

		}
	}

	// report mfe interaction
	reportOptima( outConstraint );
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2d::
initHybridE( const size_t j1, const size_t j2, const size_t i1init, const size_t i2init )
{
#if IN_DEBUG_MODE
	if (i1init > j1)
		throw std::runtime_error("PredictorMfe2d::initHybridE() : i1init > j1 : "+toString(i1init)+" > "+toString(j1));
	if (i2init > j2)
		throw std::runtime_error("PredictorMfe2d::initHybridE() : i2init > j2 : "+toString(i2init)+" > "+toString(j2));
#endif

	// global vars to avoid reallocation
	size_t i1,i2,w1,w2,k1,k2;

	// to mark as to be computed
	const E_type E_MAX = std::numeric_limits<E_type>::max();
	// to test whether computation is reasonable
	const E_type minInitDangleEndEnergy = minInitEnergy + 2.0*minDangleEnergy + 2.0*minEndEnergy;

	hybridErange.r1.from = std::max(i1init,j1-std::min(j1,energy.getAccessibility1().getMaxLength()+1));
	hybridErange.r1.to = j1;
	hybridErange.r2.from = std::max(i2init,j2-std::min(j2,energy.getAccessibility2().getMaxLength()+1));
	hybridErange.r2.to = j2;

	// temporary variable to reduce lookups
	E_type curED1 = 0.0;

	// init used part with E_INF - 1 and E_INF if ED test fails
	// iterating decreasing window size seq1 (only sane windows)
	for (i1=hybridErange.r1.from; i1<=j1; i1++ ) {
		w1 = j1-i1+1;
		curED1 = energy.getED1(i1,j1);
		// iterating decreasing window size seq2 (only sane windows)
		for (i2=hybridErange.r2.from; i2<=j2; i2++) {
			w2 = j2-i2+1;
			// check if all larger windows needing this site are already set to INF
			bool largerWindowsINF = i1==hybridErange.r1.from && i2==hybridErange.r2.from;
			// check all larger windows w1 + i2p..j2 (that might need this window for computation)
			for (size_t i2p=hybridErange.r2.from; largerWindowsINF && i2p>i2; i2p++) {
				// check if larger window is E_INF
				largerWindowsINF = E_isINF(hybridE_pq(i1,i2p));
			}
			// check all larger windows w2 + w1p (that might need this window for computation)
			for (size_t i1p=hybridErange.r1.from; largerWindowsINF && i1p>i1; i1p++) {
				// check if larger window is E_INF
				largerWindowsINF = E_isINF(hybridE_pq(i1,i2));
			}

			// if it holds for all w'>=w: ED1(i1+w1')+ED2(i2+w2') > -1*(min(w1',w2')*EmaxStacking + Einit + 2*Edangle + 2*Eend)
			// ie. the ED values exceed the max possible energy gain of an interaction
			if( largerWindowsINF
				&& ( -1.0*(std::min(w1,w2)*minStackingEnergy + minInitDangleEndEnergy)
					< (curED1 + energy.getED2(i2,j2)))
				)
			{
				// mark as NOT to be computed
				hybridE_pq(i1,i2) = E_INF;
				continue;
			}

			// mark as to be computed (has to be < E_INF)
			hybridE_pq(i1,i2) = E_MAX;
		}
	}

}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2d::
fillHybridE( const size_t j1, const size_t j2, const size_t i1init, const size_t i2init )
{

	// init for right interaction end (j1,j2)
	initHybridE( j1, j2, i1init, i2init );
	// sanity check of initialization
	assert(hybridErange.r1.to == j1);
	assert(hybridErange.r2.to == j2);

	// global vars to avoid reallocation
	size_t i1,i2,w1,w2,k1,k2;

	//////////  FIRST ROUND : COMPUTE HYBRIDIZATION ENERGIES ONLY  ////////////

	// current minimal value
	E_type curMinE = E_INF;
	// iterate over all window starts i1 (seq1) and i2 (seq2)
	// TODO PARALLELIZE THIS DOUBLE LOOP ?!
	for (i1=hybridErange.r1.to+1; i1-- > hybridErange.r1.from; ) {
		w1 = j1-i1+1;
		// w1 width check obsolete due to hybridErange setup
		// get maximal w2 for this w1
		const size_t maxW2 = getMaxInteractionWidth( w1, energy.getMaxInternalLoopSize1());
		// screen for left boundaries in seq2
		for (i2=hybridErange.r2.to+1; i2-- > hybridErange.r2.from; ) {
			w2 = j2-i2+1;
			// w2 width check obsolete due to hybridErange setup
			// check if widths' combination possible
			if ( w2 > maxW2 || w1 > getMaxInteractionWidth( w2, energy.getMaxInternalLoopSize2()) )
			{
				// combination not possible
				hybridE_pq(i1,i2) = E_INF;
				continue;
			}
			// check if left boundary (i1,i2) is complementary
			if (!energy.areComplementary(i1,i2)) {
				// interaction not possible
				hybridE_pq(i1,i2) = E_INF;
				continue;
			}

			// check if this cell is to be computed (!=E_INF)
			if( E_isNotINF( hybridE_pq(i1,i2) ) ) {

				// compute entry

				// either interaction initiation
				if ( w1==1 && w2==1 )  {
					curMinE = energy.getE_init();
				} else {
				// or only internal loop energy (nothing between i and j)
					curMinE = energy.getE_interLeft(i1,j1,i2,j2)
							+ hybridE_pq(j1,j2);
				}

				// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
				if (w1 > 2 && w2 > 2) {
					for (k1=std::min(j1-1,i1+energy.getMaxInternalLoopSize1()+1); k1>i1; k1--) {
					for (k2=std::min(j2-1,i2+energy.getMaxInternalLoopSize2()+1); k2>i2; k2--) {
						// check if (k1,k2) are valid left boundary
						if ( E_isNotINF( hybridE_pq(k1,k2) ) ) {
							curMinE = std::min( curMinE,
									(energy.getE_interLeft(i1,k1,i2,k2)
											+ hybridE_pq(k1,k2) )
									);
						}
					}
					}
				}
				// store value
				hybridE_pq(i1,i2) = curMinE;
				// update mfe if needed
				updateOptima( i1,j1,i2,j2, hybridE_pq(i1,i2), true );
				continue;
			}
		}
	}

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2d::
traceBack( Interaction & interaction )
{
	// check if something to trace
	if (interaction.basePairs.size() < 2) {
		return;
	}

#if IN_DEBUG_MODE
	// sanity checks
	if ( ! interaction.isValid() ) {
		throw std::runtime_error("PredictorMfe2d::traceBack() : given interaction not valid");
	}
	if ( interaction.basePairs.size() != 2 ) {
		throw std::runtime_error("PredictorMfe2d::traceBack() : given interaction does not contain boundaries only");
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
			j1 = energy.getIndex1(interaction.basePairs.at(1)),
			i2 = energy.getIndex2(interaction.basePairs.at(0)),
			j2 = energy.getIndex2(interaction.basePairs.at(1));


	// refill submatrix of mfe interaction
	fillHybridE( j1, j2, i1, i2 );

	// the currently traced value for i1-j1, i2-j2
	E_type curE = hybridE_pq(i1,i2);

	// trace back
	while( i1 != j1 ) {
		// check if just internal loop
		if ( E_equal( curE, (energy.getE_interLeft(i1,j1,i2,j2) )
				+ hybridE_pq(j1,j2)) )
		{
			break;
		}
		// check all interval splits
		if ( (j1-i1) > 1 && (j2-i2) > 1) {
			// temp variables
			size_t k1,k2;
			bool traceNotFound = true;
			// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
			for (k1=std::min(j1-1,i1+energy.getMaxInternalLoopSize1()+1); traceNotFound && k1>i1; k1--) {
			for (k2=std::min(j2-1,i2+energy.getMaxInternalLoopSize2()+1); traceNotFound && k2>i2; k2--) {
				// check if (k1,k2) are valid left boundary
				if ( E_isNotINF( hybridE_pq(k1,k2) ) ) {
					if ( E_equal( curE,
							(energy.getE_interLeft(i1,k1,i2,k2)
							+ hybridE_pq(k1,k2)) ) )
					{
						// stop searching
						traceNotFound = false;
						// store splitting base pair
						interaction.basePairs.push_back( energy.getBasePair(k1,k2) );
						// trace right part of split
						i1=k1;
						i2=k2;
						curE = hybridE_pq(i1,i2);
					}
				}
			}
			}
		}
	// do until only right boundary is left over
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
		(*bps.rbegin()) = energy.getBasePair(j1,j2);
	}

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2d::
getNextBest( Interaction & curBest )
{
	curBest.energy = E_INF;
	curBest.basePairs.clear();
}

////////////////////////////////////////////////////////////////////////////

