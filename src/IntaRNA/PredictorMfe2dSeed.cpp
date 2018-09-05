
#include "IntaRNA/PredictorMfe2dSeed.h"

namespace IntaRNA {

//////////////////////////////////////////////////////////////////////////

PredictorMfe2dSeed::
PredictorMfe2dSeed(
		const InteractionEnergy & energy
		, OutputHandler & output
		, PredictionTracker * predTracker
		, SeedHandler * seedHandlerInstance )
 :
	PredictorMfe2d(energy,output,predTracker)
	, seedHandler(seedHandlerInstance)
	, hybridE_pq_seed()
{
	assert( seedHandler.getConstraint().getBasePairs() > 1 );
}

//////////////////////////////////////////////////////////////////////////

PredictorMfe2dSeed::
~PredictorMfe2dSeed()
{
}

//////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dSeed::
predict( const IndexRange & r1, const IndexRange & r2
		, const OutputConstraint & outConstraint )
{
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) <<"predicting mfe interactions with seed in O(n^2) space and O(n^4) time..."; }
	// measure timing
	TIMED_FUNC_IF(timerObj,VLOG_IS_ON(9));

	// suboptimal setup check
	if (outConstraint.reportMax>1 && outConstraint.reportOverlap != OutputConstraint::ReportOverlap::OVERLAP_BOTH) {
		throw std::runtime_error("PredictorMfe2dSeed : the enumeration of non-overlapping suboptimal interactions is not supported in this prediction mode");
	}

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

	const size_t hybridE_pqsize1 = std::min( energy.size1()
			, (r1.to==RnaSequence::lastPos?energy.size1()-1:r1.to)-r1.from+1 );
	const size_t hybridE_pqsize2 = std::min( energy.size2()
			, (r2.to==RnaSequence::lastPos?energy.size2()-1:r2.to)-r2.from+1 );


	// compute seed interactions for whole range
	// and check if any seed possible
	if (seedHandler.fillSeed( 0, hybridE_pqsize1-1, 0, hybridE_pqsize2-1 ) == 0) {
		// trigger empty interaction reporting
		initOptima(outConstraint);
		reportOptima(outConstraint);
		// stop computation
		return;
	}

	// resize matrix
	hybridE_pq.resize( hybridE_pqsize1, hybridE_pqsize2 );
	hybridE_pq_seed.resize( hybridE_pqsize1, hybridE_pqsize2 );

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

			// compute hybridE_pq_seed and update mfe via PredictorMfe2d::updateOptima()
			fillHybridE_seed( j1, j2, 0, 0, outConstraint );
		}
	}

	// report mfe interaction
	reportOptima( outConstraint );

}

//////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dSeed::
fillHybridE_seed( const size_t j1, const size_t j2, const size_t i1min, const size_t i2min
		, const OutputConstraint & outConstraint )
{

	// compute hybridE_pq
	fillHybridE( j1, j2, outConstraint, i1min, i2min );

	assert(i1min <= j1);
	assert(i2min <= j2);
	assert(hybridErange.r1.from <= i1min);
	assert(hybridErange.r2.from <= i2min);
	assert(j1==hybridErange.r1.to);
	assert(j2==hybridErange.r2.to);
	assert(j1<hybridE_pq.size1());
	assert(j2<hybridE_pq.size2());

	// check if it is possible to have a seed ending on the right at (j1,j2)
	if (std::min(j1-i1min,j2-i2min)+1 < seedHandler.getConstraint().getBasePairs()) {
		// no seed possible, abort table computation
		return;
	}

	// global vars to avoid reallocation
	size_t i1,i2,k1,k2,j1c,j2c;

	// get i1/i2 index boundaries for computation
	const IndexRange i1range( std::max(hybridErange.r1.from,i1min), j1+1-seedHandler.getConstraint().getBasePairs() );
	const IndexRange i2range( std::max(hybridErange.r2.from,i2min), j2+1-seedHandler.getConstraint().getBasePairs() );

	//////////  COMPUTE HYBRIDIZATION ENERGIES (WITH SEED)  ////////////

	// current minimal value
	E_type curMinE = E_INF;
	// iterate over all window starts i1 (seq1) and i2 (seq2)
	// TODO PARALLELIZE THIS DOUBLE LOOP ?!
	for (i1=1+i1range.to; i1-- > i1range.from; ) {
		// screen for left boundaries i2 in seq2
		for (i2=1+i2range.to; i2-- > i2range.from; ) {

			// compute entry
			curMinE = E_INF;

			// check if this cell is to be computed (!=E_INF)
			if( E_isNotINF( hybridE_pq(i1,i2) ) ) {

				// base case = incorporate mfe seed starting at (i1,i2)
				//             + interaction on right side up to (p,q)
				if ( E_isNotINF( seedHandler.getSeedE(i1,i2) ) ) {
					// decode right mfe boundary
					k1 = i1+seedHandler.getSeedLength1(i1,i2)-1;
					k2 = i2+seedHandler.getSeedLength2(i1,i2)-1;
					// compute overall energy of seed+upToPQ
					if ( k1 <= j1 && k2 <= j2 && E_isNotINF(hybridE_pq(k1,k2))) {
						curMinE = seedHandler.getSeedE(i1,i2) + hybridE_pq(k1,k2);
					}
				}

				// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
				// where k1..j1 contains a seed
				for (k1=std::min(i1range.to,i1+energy.getMaxInternalLoopSize1()+1); k1>i1; k1--) {
				for (k2=std::min(i2range.to,i2+energy.getMaxInternalLoopSize2()+1); k2>i2; k2--) {
					// check if (k1,k2) are valid left boundaries including a seed
					if ( E_isNotINF( hybridE_pq_seed(k1,k2) ) ) {
						curMinE = std::min( curMinE,
								(energy.getE_interLeft(i1,k1,i2,k2)
										+ hybridE_pq_seed(k1,k2) )
							);
					}
				}
				}
				// update mfe if needed (call super class)
				if (E_isNotINF(curMinE)) {
					PredictorMfe2d::updateOptima( i1,j1,i2,j2, curMinE, true );
				}
			}

			// store value
			hybridE_pq_seed(i1,i2) = curMinE;
		}
	}

}

//////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dSeed::
traceBack( Interaction & interaction, const OutputConstraint & outConstraint  )
{
	// check if something to trace
	if (interaction.basePairs.size() < 2) {
		return;
	}

#if INTARNA_IN_DEBUG_MODE
	// sanity checks
	if ( interaction.basePairs.size() != 2 ) {
		throw std::runtime_error("PredictorMfe2dSeed::traceBack() : given interaction does not contain boundaries only");
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
		throw std::runtime_error("PredictorMfe2dSeed::traceBack() : given interaction not valid");
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
		throw std::runtime_error("PredictorMfe2dSeed::traceBack() : given boundaries "+toString(interaction)+" can not hold a seed of "+toString(seedHandler.getConstraint().getBasePairs())+" base pairs");
	}
#endif

	// temp variables
	size_t k1,k2;


	// refill submatrices of mfe interaction
	fillHybridE_seed( j1, j2, i1, i2, outConstraint );

	// the currently traced value for i1-j1, i2-j2
	E_type curE = hybridE_pq_seed(i1,i2);

	// trace back
	bool seedNotTraced = true;
	while( i1 != j1 ) {

		// check if we still have to find the seed
		if (seedNotTraced) {

			// check base case == seed only
			if ( E_isNotINF( seedHandler.getSeedE(i1,i2) ) ) {

				// right boundary of seed
				k1 = i1 + seedHandler.getSeedLength1(i1,i2) -1;
				k2 = i2 + seedHandler.getSeedLength2(i1,i2) -1;

				// check if correct trace
				if ( E_equal( curE, seedHandler.getSeedE(i1,i2) + hybridE_pq(k1,k2) ) ) {
					// store seed information
					interaction.setSeedRange(
									energy.getBasePair(i1,i2),
									energy.getBasePair(k1,k2),
									energy.getE(i1,k1,i2,k2,seedHandler.getSeedE(i1,i2))+energy.getE_init());
					// trace back seed base pairs
					seedHandler.traceBackSeed( interaction, i1, i2 );
					// continue after seed
					i1 = k1;
					i2 = k2;
					curE = hybridE_pq(k1,k2);
					seedNotTraced = false;
					continue;
				}
			}
			// check all interval splits
			if ( (j1-i1) > 1 && (j2-i2) > 1) {
				// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
				// where k1..j1 contains a seed
				bool traceNotFound = true;
				for (k1=std::min(j1-seedHandler.getConstraint().getBasePairs()+1,i1+energy.getMaxInternalLoopSize1()+1); traceNotFound && k1>i1; k1--) {
				for (k2=std::min(j2-seedHandler.getConstraint().getBasePairs()+1,i2+energy.getMaxInternalLoopSize2()+1); traceNotFound && k2>i2; k2--) {
					// check if (k1,k2) are valid left boundaries including a seed
					if ( E_isNotINF( hybridE_pq_seed(k1,k2) ) ) {
						// check if correct split
						if (E_equal ( curE,
								(energy.getE_interLeft(i1,k1,i2,k2)
										+ hybridE_pq_seed(k1,k2) )
								) )
						{
							// update trace back boundary
							i1=k1;
							i2=k2;
							curE= hybridE_pq_seed(k1,k2);
							// stop search splits
							traceNotFound = false;
							// store splitting base pair
							interaction.basePairs.push_back( energy.getBasePair(k1,k2) );
						}
					}
				} // k2
				} // k1
				assert(!traceNotFound);
			}
		}
		// seed was already traced, do "normal" interaction trace
		else {
			// create temporary data structure to be filed
			Interaction rightSide( *interaction.s1, *interaction.s2 );
			rightSide.basePairs.push_back( energy.getBasePair(i1,i2) );
			rightSide.basePairs.push_back( energy.getBasePair(j1,j2) );
			// call traceback of super class
			PredictorMfe2d::traceBack( rightSide, outConstraint );
			// copy base pairs (excluding last)
			for (size_t i=0; i+1<rightSide.basePairs.size(); i++) {
				interaction.basePairs.push_back( rightSide.basePairs.at(i) );
			}
			i1 = j1;
			i2 = j2;
			// stop traceback
			break;
		}
	// do until only right boundary is left over (already part of interaction)
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
PredictorMfe2dSeed::
getNextBest( Interaction & curBest )
{
	throw std::runtime_error("PredictorMfe2dSeed::getNextBest() : This prediction mode does not support non-overlapping suboptimal interaction enumeration.");
}

//////////////////////////////////////////////////////////////////////////




} // namespace
