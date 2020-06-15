
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
predict( const IndexRange & r1, const IndexRange & r2  )
{
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) <<"predicting mfe interactions with seed in O(n^2) space and O(n^4) time..."; }
	// measure timing
	TIMED_FUNC_IF(timerObj,VLOG_IS_ON(9));

#if INTARNA_IN_DEBUG_MODE
	// check indices
	if (!(r1.isAscending() && r2.isAscending()) )
		throw std::runtime_error("PredictorMfe2dSeed::predict("+toString(r1)+","+toString(r2)+") is not sane");
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
		initOptima();
		reportOptima();
		// stop computation
		return;
	}

	// resize matrix
	hybridE_pq.resize( hybridE_pqsize1, hybridE_pqsize2 );
	hybridE_pq_seed.resize( hybridE_pqsize1, hybridE_pqsize2 );

	// initialize mfe interaction for updates
	initOptima();

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

			// compute both hybridE_pq and hybridE_pq_seed and update mfe
			fillHybridE( j1, j2, 0, 0, true );
		}
	}

	// report mfe interaction
	reportOptima();

}

//////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dSeed::
fillHybridE( const size_t j1, const size_t j2
		, const size_t i1min, const size_t i2min
		, const bool callUpdateOptima)
{
	// temporary access
	const OutputConstraint & outConstraint = output.getOutputConstraint();

#if INTARNA_IN_DEBUG_MODE
	if (i1min > j1)
		throw std::runtime_error("PredictorMfe2dSeed::fillHybridE() : i1init > j1 : "+toString(i1min)+" > "+toString(j1));
	if (i2min > j2)
		throw std::runtime_error("PredictorMfe2dSeed::fillHybridE() : i2init > j2 : "+toString(i2min)+" > "+toString(j2));
	assert(j1<hybridE_pq.size1());
	assert(j2<hybridE_pq.size2());
#endif

	// get minimal start indices heeding max interaction length
	const size_t i1start = std::max(i1min,j1-std::min(j1,energy.getAccessibility1().getMaxLength()+1));
	const size_t i2start = std::max(i2min,j2-std::min(j2,energy.getAccessibility2().getMaxLength()+1));


	// check if it is possible to have a seed within the range
	if (std::min(j1-i1min,j2-i2min)+1 < seedHandler.getConstraint().getBasePairs()) {
		// no seed possible, abort table computation
		return;
	}

	// global vars to avoid reallocation
	size_t i1,i2,k1,k2,j1c,j2c, w1,w2;

	// determine whether or not lonely base pairs are allowed or if we have to
	// ensure a stacking to the right of the left boundary (i1,i2)
	const size_t noLpShift = outConstraint.noLP ? 1 : 0;
	E_type iStackE = E_type(0);

	//////////  COMPUTE HYBRIDIZATION ENERGIES (WITHOUT AND WITH SEED)  ////////////

	// current minimal value
	E_type curEloop = E_INF;
	// iterate over all window starts i1 (seq1) and i2 (seq2)
	for (i1=1+j1; i1-- > i1start; ) {
		// w1 = interaction width in seq1
		w1 = j1-i1+1;
		// screen for left boundaries i2 in seq2
		for (i2=1+j2; i2-- > i2start; ) {

			// init: mark as invalid boundary
			hybridE_pq(i1,i2) = E_INF;
			hybridE_pq_seed(i1,i2) = E_INF;

			// check if this cell is to be computed (!=E_INF)
			if( energy.isAccessible1(i1)
				&& energy.isAccessible2(i2)
				&& energy.areComplementary(i1,i2))
			{

				// w2 = interaction width in seq2
				w2 = j2-i2+1;

				// reference access to cell value
				E_type &curMinE = hybridE_pq(i1,i2);
				E_type &curMinEseed = hybridE_pq_seed(i1,i2);

				// either interaction initiation
				if ( i1==j1 && i2==j2 )  {
					// single base pair
					if (noLpShift == 0) {
						curMinE = energy.getE_init();
					}
				} else
				// or at least two base pairs possible
				if ( w1 > 1 && w2 > 1) { // or more complex stuff

					// no-seed base case:
					// if lonely bps are allowed
					if (noLpShift == 0) {
						// test full-width internal loop energy (nothing between i and j)
						// will be E_INF if loop is too large
						curMinE = energy.getE_interLeft(i1,j1,i2,j2) + hybridE_pq(j1,j2);
					} else {
						// no lp allowed
						// check if right-side stacking of (i1,i2) is possible
						if (energy.isAccessible1(i1+noLpShift)
							&& energy.isAccessible2(i2+noLpShift)
							&& energy.areComplementary(i1+noLpShift,i2+noLpShift))
						{
							// get stacking term to avoid recomputation
							iStackE = energy.getE_interLeft(i1,i1+noLpShift,i2,i2+noLpShift);

							// init with stacking only
							curMinE = iStackE + ((w1==2&&w2==2) ? energy.getE_init() : hybridE_pq(i1+noLpShift, i2+noLpShift) );
						} else {
							//
							iStackE = E_INF;
						}
					}

					// seed base case = incorporate mfe seed starting at (i1,i2)
					//             + interaction on right side up to (p,q)
					// assumption: seed heeds noLP constraint
					if ( seedHandler.isSeedBound(i1,i2) ) {
						// decode right mfe boundary
						k1 = i1+seedHandler.getSeedLength1(i1,i2)-1;
						k2 = i2+seedHandler.getSeedLength2(i1,i2)-1;
						// compute overall energy of seed+upToPQ
						if ( k1 <= j1 && k2 <= j2 ) {
							// init = seed + something without seed
							if ( E_isNotINF(hybridE_pq(k1,k2)) ) {
								// Note: noLP-handling via hybridE_pq recursion
								curMinEseed = std::min( curMinEseed, seedHandler.getSeedE(i1,i2) + hybridE_pq(k1,k2) );
							} else
							// just the seed up to right boundary (explicit noLP handling)
							if (k1 == j1 && k2 == j2) {
								curMinEseed = std::min( curMinEseed, seedHandler.getSeedE(i1,i2) + energy.getE_init() );
							}
							// handle interior loops after seeds in noLP-mode
							if ( outConstraint.noLP) {
								// multi-bp seed can be followed directly by interior loop
								// check all combinations of decompositions into (k1,k2)..(l1,l2)-(j1,j2)
								if (k1 > i1 && j1-k1 > 2 && j2-k2 > 2) {
									for (size_t l1=std::min(j1-1,k1+energy.getMaxInternalLoopSize1()+1); l1>k1; l1--) {
									for (size_t l2=std::min(j2-1,k2+energy.getMaxInternalLoopSize2()+1); l2>k2; l2--) {
										// check if (k1,k2) are valid left boundary
										if ( E_isNotINF( hybridE_pq(l1,l2) ) ) {
											// get loop energy
											curEloop = energy.getE_interLeft(k1,l1,k2,l2);
											// update minimal value
											curMinEseed = std::min( curMinEseed, (seedHandler.getSeedE(i1,i2) + curEloop + hybridE_pq(l1,l2) ) );
										}
									}
									}
								}
							}
						}
					}

					// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
					if (w1 > 2 && w2 > 2) {
						for (k1=std::min(j1-1,i1+energy.getMaxInternalLoopSize1()+1+noLpShift); k1>i1+noLpShift; k1--) {
						for (k2=std::min(j2-1,i2+energy.getMaxInternalLoopSize2()+1+noLpShift); k2>i2+noLpShift; k2--) {
							// check if (k1,k2) are valid left boundary
							if ( E_isNotINF( hybridE_pq(k1,k2) ) ) {
								// get loop energy
								curEloop = energy.getE_interLeft(i1+noLpShift,k1,i2+noLpShift,k2);
								// update minimal value
								curMinE = std::min( curMinE, (iStackE + curEloop + hybridE_pq(k1,k2) ) );
								// check if (k1,k2) are valid left boundaries including a seed
								if ( E_isNotINF( hybridE_pq_seed(k1,k2) ) ) {
									curMinEseed = std::min( curMinEseed, (iStackE + curEloop + hybridE_pq_seed(k1,k2) ) );
								}
							}
						}
						}
					}
				}

				// update mfe if needed
				if (callUpdateOptima) {
					updateOptima( i1,j1,i2,j2, curMinEseed, true, true );
				}

			} // complementary base pair
		}
	}

}

//////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dSeed::
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
		throw std::runtime_error("PredictorMfe2dSeed::traceBack() : given interaction does not contain boundaries only");
	}
#endif

	// ensure sorting
	interaction.sort();

	// check for single base pair interaction
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

	// get indices in hybridE for boundary base pairs
	size_t	i1 = energy.getIndex1(interaction.basePairs.at(0)),
			j1 = energy.getIndex1(interaction.basePairs.at(1)),
			i2 = energy.getIndex2(interaction.basePairs.at(0)),
			j2 = energy.getIndex2(interaction.basePairs.at(1));

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
	fillHybridE( j1, j2, i1, i2, false );

	// the currently traced value for i1-j1, i2-j2
	E_type curE = hybridE_pq_seed(i1,i2);

	// determine whether or not lonely base pairs are allowed or if we have to
	// ensure a stacking to the right of the left boundary (i1,i2)
	const size_t noLpShift = outConstraint.noLP ? 1 : 0;
	E_type iStackE = E_type(0);

	// trace back
	bool seedNotTraced = true;
	while( i1 < j1 ) {

		// check if only 3 stack helix possible
		if (outConstraint.noLP && j1-i1+1 == 3) {
			interaction.basePairs.push_back( energy.getBasePair(i1+1,i2+1) );
			break;
		}


		// check if we still have to find the seed
		if (seedNotTraced) {

			// handle single-bp seed
			if (outConstraint.noLP) {
				iStackE = energy.getE_interLeft(i1,i1+noLpShift,i2,i2+noLpShift);
			}

			// check base case == seed only
			if ( seedHandler.isSeedBound(i1,i2) ) {

				// right boundary of seed
				k1 = i1 + seedHandler.getSeedLength1(i1,i2) -1;
				k2 = i2 + seedHandler.getSeedLength2(i1,i2) -1;
				// energy of seed
				const E_type seedE = seedHandler.getSeedE(i1,i2);

				if(k1 == j1 && k2 == j2 && E_equal(curE, seedE + energy.getE_init())) {
					// just seed
					// store seed information
					interaction.setSeedRange(
									energy.getBasePair(i1,i2),
									energy.getBasePair(k1,k2),
									energy.getE(i1,k1,i2,k2,seedE)+energy.getE_init());
					// trace back seed base pairs
					seedHandler.traceBackSeed( interaction, i1, i2 );
					seedNotTraced = false;
					// stop tracing
					break;
				} else
				// check if seed + something without seed
				if ( E_equal( curE,  seedE + hybridE_pq(k1,k2) ) ) {
					// store seed information
					interaction.setSeedRange(
									energy.getBasePair(i1,i2),
									energy.getBasePair(k1,k2),
									energy.getE(i1,k1,i2,k2,seedE)+energy.getE_init());
					// trace back seed base pairs
					seedHandler.traceBackSeed( interaction, i1, i2 );
					seedNotTraced = false;
					// continue after seed
					i1 = k1;
					i2 = k2;
					curE = hybridE_pq(i1,i2);
					continue;
				} else
				// check multi-bp seed followed by interior loop
				if (outConstraint.noLP && k1 > i1 && j1-k1 > 2 && j2-k2 > 2) {
					bool traceNotFound = true;
					for (size_t l1=std::min(j1-1,k1+energy.getMaxInternalLoopSize1()+1); traceNotFound && l1>k1; l1--) {
					for (size_t l2=std::min(j2-1,k2+energy.getMaxInternalLoopSize2()+1); traceNotFound && l2>k2; l2--) {
						// check if (k1,k2) are valid left boundaries including a seed
						if ( E_isNotINF( hybridE_pq(l1,l2) ) ) {
							// check if correct split
							if (E_equal ( curE,
									(seedE + energy.getE_interLeft(k1,l1,k2,l2)
											+ hybridE_pq(l1,l2) )
									) )
							{
								// store seed information
								interaction.setSeedRange(
												energy.getBasePair(i1,i2),
												energy.getBasePair(k1,k2),
												energy.getE(i1,k1,i2,k2,seedE)+energy.getE_init());
								// trace back seed base pairs
								seedHandler.traceBackSeed( interaction, i1, i2 );
								seedNotTraced = false;
								// store splitting base pair
								interaction.basePairs.push_back( energy.getBasePair(k1,k2) );
								// update trace back boundary
								i1=l1;
								i2=l2;
								curE = hybridE_pq(i1,i2);
								// stop search splits
								traceNotFound = false;
							}
						}
					} // l2
					} // l1
					// start next iteration if trace was found
					if (!traceNotFound) {
						continue;
					}
				}
			}

			// check all interval splits
			if ( (j1-i1) > 1 && (j2-i2) > 1) {
				// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
				// where k1..j1 contains a seed
				bool traceNotFound = true;
				for (k1=std::min(j1-seedHandler.getConstraint().getBasePairs()+1,i1+energy.getMaxInternalLoopSize1()+1+noLpShift); traceNotFound && k1>i1+noLpShift; k1--) {
				for (k2=std::min(j2-seedHandler.getConstraint().getBasePairs()+1,i2+energy.getMaxInternalLoopSize2()+1+noLpShift); traceNotFound && k2>i2+noLpShift; k2--) {
					// check if (k1,k2) are valid left boundaries including a seed
					if ( E_isNotINF( hybridE_pq_seed(k1,k2) ) ) {
						// check if correct split
						if (E_equal ( curE,
								(iStackE + energy.getE_interLeft(i1+noLpShift,k1,i2+noLpShift,k2)
										+ hybridE_pq_seed(k1,k2) )
								) )
						{
							if (outConstraint.noLP) {
								// store stacked base pair
								interaction.basePairs.push_back( energy.getBasePair(i1+noLpShift,i2+noLpShift) );
							}
							// update trace back boundary
							i1=k1;
							i2=k2;
							curE= hybridE_pq_seed(i1,i2);
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
			PredictorMfe2d::traceBack( rightSide );
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

	// add all seeds that are subsets of the interaction
	seedHandler.addSeeds( interaction );
}

////////////////////////////////////////////////////////////////////////////




} // namespace
