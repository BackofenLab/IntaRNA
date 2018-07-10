
#include "IntaRNA/PredictorMfe2d.h"

#include <stdexcept>

namespace IntaRNA {

////////////////////////////////////////////////////////////////////////////

PredictorMfe2d::
PredictorMfe2d(
		const InteractionEnergy & energy
		, OutputHandler & output
		, PredictionTracker * predTracker )
 : PredictorMfe(energy,output,predTracker)
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
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) <<"predicting mfe interactions in O(n^2) space..."; }
	// measure timing
	TIMED_FUNC_IF(timerObj,VLOG_IS_ON(9));

	// suboptimal setup check
	if (outConstraint.reportMax>1 && outConstraint.reportOverlap != OutputConstraint::ReportOverlap::OVERLAP_BOTH) {
		throw std::runtime_error("PredictorMfe2d : the enumeration of non-overlapping suboptimal interactions is not supported in this prediction mode");
	}

#if INTARNA_IN_DEBUG_MODE
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
			fillHybridE( j1, j2, outConstraint, 0, 0 );

		}
	}

	// report mfe interaction
	reportOptima( outConstraint );
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2d::
initHybridE( const size_t j1, const size_t j2
			, const OutputConstraint & outConstraint
			, const size_t i1init, const size_t i2init
			)
{
#if INTARNA_IN_DEBUG_MODE
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

	for (i1=hybridErange.r1.from; i1<=j1; i1++ ) {
		for (i2=hybridErange.r2.from; i2<=j2; i2++) {
			// check if complementary, i.e. to be computed
			if( energy.areComplementary(i1,i2) )
			{
				// mark as to be computed (has to be < E_INF)
				hybridE_pq(i1,i2) = E_MAX;
			} else {
				// mark as NOT to be computed
				hybridE_pq(i1,i2) = E_INF;
			}
		}
	}

}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2d::
fillHybridE( const size_t j1, const size_t j2
			, const OutputConstraint & outConstraint
			, const size_t i1init, const size_t i2init )
{

	// init for right interaction end (j1,j2)
	initHybridE( j1, j2, outConstraint, i1init, i2init );
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
		// screen for left boundaries in seq2
		for (i2=hybridErange.r2.to+1; i2-- > hybridErange.r2.from; ) {
			// w2 width check obsolete due to hybridErange setup
			w2 = j2-i2+1;
			curMinE = E_INF;

			// check if this cell is to be computed (!=E_INF)
			if( E_isNotINF( hybridE_pq(i1,i2) ) ) {

				// compute entry

				// either interaction initiation
				if ( i1==j1 && i2==j2 )  {
					curMinE = energy.getE_init();
				} else { // or more complex stuff
					// test only internal loop energy (nothing between i and j)
					// will be E_INF if loop is too large
					curMinE = energy.getE_interLeft(i1,j1,i2,j2)
							+ hybridE_pq(j1,j2);

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
traceBack( Interaction & interaction, const OutputConstraint & outConstraint  )
{
	// check if something to trace
	if (interaction.basePairs.size() < 2) {
		return;
	}

#if INTARNA_IN_DEBUG_MODE
	// sanity checks
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

#if INTARNA_IN_DEBUG_MODE
	// sanity checks
	if ( ! interaction.isValid() ) {
		throw std::runtime_error("PredictorMfe2d::traceBack() : given interaction is not valid");
	}
#endif

	// ensure sorting
	interaction.sort();
	// get indices in hybridE for boundary base pairs
	size_t	i1 = energy.getIndex1(interaction.basePairs.at(0)),
			j1 = energy.getIndex1(interaction.basePairs.at(1)),
			i2 = energy.getIndex2(interaction.basePairs.at(0)),
			j2 = energy.getIndex2(interaction.basePairs.at(1));

	// refill submatrix of mfe interaction
	fillHybridE( j1, j2, outConstraint, i1, i2 );

	// the currently traced value for i1-j1, i2-j2
	E_type curE = hybridE_pq(i1,i2);

	// trace back
	while( i1 != j1 ) {

		// check if just internal loop
		if ( E_equal( curE, (energy.getE_interLeft(i1,j1,i2,j2) + hybridE_pq(j1,j2)) ) )
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
							(energy.getE_interLeft(i1,k1,i2,k2) + hybridE_pq(k1,k2)) ) )
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
	// do until only right boundary base pair is left over
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


} // namespace
