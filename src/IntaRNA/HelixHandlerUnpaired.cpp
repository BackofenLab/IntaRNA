#include "IntaRNA/HelixHandlerUnpaired.h"

namespace IntaRNA {

/////////////////////////////////////////////////////////////////////////////

size_t
HelixHandlerUnpaired::
fillHelix(const size_t i1min, const size_t i1max, const size_t i2min, const size_t i2max)
{
#if INTARNA_IN_DEBUG_MODE
	if ( i1min > i1max ) throw std::runtime_error("HelixHandlerUnpaired::fillHelix: i1min("+toString(i1min)+") > i1max("+toString(i1max)+")");
	if ( i2min > i2max ) throw std::runtime_error("HelixHandlerUnpaired::fillHelix: i2min("+toString(i2min)+") > i2max("+toString(i2max)+")");
	if ( i1max > energy.size1() ) throw std::runtime_error("HelixHandlerUnpaired::fillHelix: i1max("+toString(i1max)+") > energy.size1("+toString(energy.size1())+")");
	if ( i2max > energy.size2() ) throw std::runtime_error("HelixHandlerUnpaired::fillHelix: i2max("+toString(i2max)+") > energy.size2("+toString(energy.size2())+")");
	if ( helixConstraint.getMinBasePairs() > helixConstraint.getMaxBasePairs() )
		throw std::runtime_error("HelixHandlerUnpaired::fillHelix: bpMin("+toString(helixConstraint.getMinBasePairs()) +") > bpMax("+toString(helixConstraint.getMaxBasePairs())+")");
#endif

	// measure timing
	TIMED_FUNC_IF(timerObj,VLOG_IS_ON(9));

	helix.resize( i1max-i1min+1, i2max-i2min+1 );
	helixE_rec.resize( HelixIndex({{
				   (HelixRecMatrix::index)(helix.size1())
				   , (HelixRecMatrix::index)(helix.size2())
				   , (HelixRecMatrix::index)(getConstraint().getMaxBasePairs()+1)}}) );
	std::fill( helixE_rec.origin(), helixE_rec.origin() + helixE_rec.num_elements(), std::make_pair(E_INF,0));

	size_t minBP;
	E_type maxE;
	// if no seedHandler is given use given values,
	// else use default values to ensure they dont compromise the results of the helixSeed method
	if (seedHandler == NULL) {
		minBP = helixConstraint.getMinBasePairs();
		maxE = helixConstraint.getMaxE();
	} else {
		minBP = 2;
		maxE = 999;
	}

	// store index offset due to restricted matrix size generation
	offset1 = i1min;
	offset2 = i2min;

	// temporary variables
	size_t i1, i2, curBP, bestBP, u1, u2, j1, j2, hL1, hL2, k1, k2, u1best, u2best, bestL1, bestL2;
	E_type curE, bestE;

	size_t  helixCountNotInf = 0;

	// fill for all start indeices
	// in decreasing index order
	for (i1=i1max+1; i1-- > i1min;) {
	for (i2=i2max+1; i2-- > i2min;) {

		// init according to no helix interaction
		helix(i1-offset1,i2-offset2) = HelixMatrix::value_type( E_INF, 0, 0 );

		// skip non-complementary left helix boundaries
		if (!energy.areComplementary(i1,i2)) {
			continue; // go to next helixE index
		}

		// Calculate energy for all different numbers of base pairs (bpMin to bpMax)
		for (curBP=2; curBP < getConstraint().getMaxBasePairs()+1
					  && (i1+curBP-1-offset1) < helix.size1()
					  && (i2+curBP-1-offset2) < helix.size2(); curBP++) {

			// init current helix energy
			curE = E_INF;
			bestE = E_INF;
			bestL1 = 0;
			bestL2 = 0;

			// for feasible unpaired bases
			for (u1 = 0; u1 < getConstraint().getMaxIL()+1 && (i1+u1+1 -offset1)< helix.size1(); u1++) {
			for (u2 = 0; u2 < getConstraint().getMaxIL()+1 - u1 && (i2+u2+1 -offset2) < helix.size2(); u2++) {

				// get split base pair (right boundaries when curBP = 2)
				k1 = i1 + u1 +1;
				k2 = i2 + u2 +1;

				// check if split base pair is complementary
				if (energy.areComplementary(k1, k2)) {

					// base case: only left and right base pair present
					if (curBP == 2) {
						// energy for stacking/bulge/interior depending on u1/u2
						curE = energy.getE_interLeft(i1, k1, i2, k2);
						// save the best energy among all u1/u2 combinations
						if (curE < bestE) {
							bestE = curE;
							bestL1 = curBP +u1;
							bestL2 = curBP +u2;
						}
					} else {
						// check if recursed entry is < E_INF
						if (E_isINF(getHelixE(k1 - offset1, k2 - offset2, curBP - 1))) {
							continue; // invalid entry -> skip
						}

						// check right boundaries
						if (i1 + u1 + getHelixLength1(k1-offset1, k2-offset2, curBP-1) -offset1 >= helix.size1()
							|| i2 + u2 + getHelixLength2(k1-offset1, k2-offset2, curBP-1)-offset2 >= helix.size2()) {
							continue; // not within the boundaries -> skip
						}

						// update mfe for split at k1,k2
						curE = energy.getE_interLeft(i1, k1, i2, k2) + getHelixE(k1 - offset1, k2 - offset2, curBP - 1);

						// store best energy only
						if (curE < bestE) {
							bestE = curE;
							bestL1 = u1 + getHelixLength1(k1 - offset1, k2 - offset2, curBP - 1)+1;
							bestL2 = u2 + getHelixLength2(k1 - offset1, k2 - offset2, curBP - 1)+1;

						}
					} // more than two base pairs
				} // (j1, j2) complementary


			} // u2
			} // u1
			// store helix energy and length
			setHelixPair(i1 - offset1, i2 - offset2, curBP, bestE, encodeHelixLength(bestL1, bestL2));

		} // curBP

		// TODO: Possible runtime improvement by checking this while calculating
		// find best unpaired combination in helx for i1,i2,bp
		bestBP = 0;
		bestE = E_INF;
		E_type bestEfullDelta = E_type(0), curEfullDelta = E_type(0);


		// Calculate energy for all different numbers of base pairs
		// Ensuring minimum number of base pairs here
		for (curBP = minBP; curBP < helixConstraint.getMaxBasePairs() + 1
						&& (i1 + curBP - 1 - offset1) < helix.size1()
						&& (i2 + curBP - 1 - offset2) < helix.size2(); curBP++) {

			if (E_isINF(getHelixE(i1-offset1, i2-offset2, curBP))) {
				continue;
			}
			// get right helix boundaries
			j1 = i1 + getHelixLength1(i1-offset1, i2-offset2, curBP)-1;
			j2 = i2 + getHelixLength2(i1-offset1, i2-offset2, curBP)-1;

			if (energy.getED1(i1, j1) > helixConstraint.getMaxED()
				|| energy.getED2(i2, j2) > helixConstraint.getMaxED()) {
				continue;
			}

			curE = getHelixE(i1 - offset1, i2 - offset2, curBP);
			// get additional energy terms
			if (helixConstraint.evalFullE()) {
				curEfullDelta = energy.getE(i1, j1, i2, j2, E_type(0)) + energy.getE_init();
			}

			// check if within energy bound
			// check if better than what is known so far
			if ( (curE+curEfullDelta) < maxE
					&& (curE+curEfullDelta) < (bestE+bestEfullDelta)
					&& !E_equal(curE+curEfullDelta, bestE+bestEfullDelta))
			{
				bestE = curE;
				bestEfullDelta = curEfullDelta;
				bestBP = curBP;
			}
		} // curBP

		// reduce bestE to hybridization energy only (init+loops)
		if (E_isNotINF(bestE)) {
			// count true helix
			helixCountNotInf++;
		}


		// store best (mfe) helix for all u1/u2
		helix(i1 - offset1, i2 - offset2) = HelixMatrix::value_type(bestE,
																	E_isINF(bestE) ? 0 : encodeHelixLength(
																			getHelixLength1(i1-offset1, i2-offset2, bestBP),
																			getHelixLength2(i1-offset1, i2-offset2, bestBP)),
																	E_isINF(bestE) ? 0: bestBP);


	} // i2
	} // i1

#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) << "valid helices = " << helixCountNotInf; }

	return helixCountNotInf;
}

//////////////////////////////////////////////////////////////////////////

void
HelixHandlerUnpaired::
traceBackHelix( Interaction & interaction
			, const size_t i1_
			, const size_t i2_
			, const size_t bp)
{
	// get boundaries
	size_t 	  i1 = i1_
			, i2 = i2_
			, u1, u2
			, k1, k2
			;

	// get energy of provided seed
	E_type curE = getHelixE(i1_,i2_,bp);
	// trace helices
	// trace each helix base pair (excluding right most)
	for ( size_t curBP=1+bp; curBP-- > 2; ) {
 		// base case: only left and right base pair present
		if (curBP==2) {
			// add left base pair if not left helix boundary
			if (i1 != i1_) {
				interaction.basePairs.push_back( energy.getBasePair( i1+offset1, i2+offset2 ) );
			}

		} else {
			// split helix recursively into all possible leading interior loops
			// i1 .. i1+u1+1 .. j1
			// i2 .. i2+u2+1 .. j2
			bool traceNotFound = true;
			for (u1=0; u1 < getConstraint().getMaxIL()+1 && traceNotFound; u1++) {
			for (u2=0; u2 < getConstraint().getMaxIL()+1-u1 && traceNotFound; u2++) {
				k1 = i1+u1+1;
				k2 = i2+u2+1;

				// check right boundary
				if ( k1 >= helix.size1() || k2 >= helix.size2()) {
					continue;
				}

				// check if valid trace
				if ( E_isNotINF( getHelixE( k1, k2, curBP-1) ) ) {

					// check if correct trace
					if ( E_equal( curE, energy.getE_interLeft(i1+offset1, k1+offset1, i2+offset2, k2+offset2)
										+ getHelixE( k1, k2, curBP-1 )) )
					{
						// store left base pair if not left helix boundary
						if (i1 != i1_) {
							interaction.basePairs.push_back( energy.getBasePair(i1+offset1, i2+offset2) );
						}

						// store next energy value in trace
						curE = getHelixE( k1, k2, curBP-1 );
						// reset for next trace step
						i1 = k1;
						i2 = k2;

						// mark trace step done
						traceNotFound = false;
					}
				}

			} // u2
			} // u1
			assert( !traceNotFound ); // sanity check
		} // more than two base pairs

	} // bp
}

/////////////////////////////////////////////////////////////////////////////

size_t
HelixHandlerUnpaired::
fillHelixSeed(const size_t i1min, const size_t i1max, const size_t i2min, const size_t i2max)
{
#if INTARNA_IN_DEBUG_MODE
	if ( i1min > i1max ) throw std::runtime_error("HelixHandlerUnpaired::fillHelixSeed: i1min("+toString(i1min)+") > i1max("+toString(i1max)+")");
	if ( i2min > i2max ) throw std::runtime_error("HelixHandlerUnpaired::fillHelixSeed: i2min("+toString(i2min)+") > i2max("+toString(i2max)+")");
	if ( i1max > energy.size1() ) throw std::runtime_error("HelixHandlerUnpaired::fillHelixSeed: i1max("+toString(i1max)+") > energy.size1("+toString(energy.size1())+")");
	if ( i2max > energy.size2() ) throw std::runtime_error("HelixHandlerUnpaired::fillHelixSeed: i2max("+toString(i2max)+") > energy.size2("+toString(energy.size2())+")");
	if ( helixConstraint.getMinBasePairs() > helixConstraint.getMaxBasePairs() )
		throw std::runtime_error("HelixHandlerUnpaired::fillHelixSeed: bpMin("+toString(helixConstraint.getMinBasePairs()) +") > bpMax("+toString(helixConstraint.getMaxBasePairs())+")");
	if ( seedHandler == NULL ) throw std::runtime_error("HelixHandlerUnpaired::fillHelixSeed() no SeedHandler available");
#endif

	// measure timing
	TIMED_FUNC_IF(timerObj,VLOG_IS_ON(9));

	// check if we can abort since no seed fits into any allowed helix
	if (seedHandler->getConstraint().getBasePairs() > helixConstraint.getMaxBasePairs()) {
		return 0;
	}

	// resize matrix for filling
	helixSeed.resize( i1max-i1min+1, i2max-i2min+1, false );

	// store index offset due to restricted matrix size generation
	offset1 = i1min;
	offset2 = i2min;

	// temporary variables
	size_t i1, i2, seedStart1, seedStart2, seedEnd1, seedEnd2, j1, j2, bestL1, bestL2, possibleBasePairs;
	size_t  helixCountNotInf = 0, helixCount = 0;

	E_type curE_withED, curE, bestE_withED, bestE;

	// fill for all start indices
	// in decreasing index order
	for (i1=i1max+1; i1-- > i1min;) {
	for (i2=i2max+1; i2-- > i2min;) {

		// count possible helices
		helixCount++;

		// init according to no helix interaction
		helixSeed(i1 - offset1, i2 - offset2) = HelixSeedMatrix::value_type(E_INF, 0);

		// skip non-complementary left helix boundaries
		if (!energy.areComplementary(i1, i2)) {
			continue; // go to next helixSeedE index
		}

		// Check if a seed can fit and calculate how many base pairs are theoretically possible
		// Note: If seedHandler allows unpaired positions this check is not enough, check happens in loop
		if (std::min(helixSeed.size1()-i1+offset1, helixSeed.size2()-i2+offset2) < seedHandler->getConstraint().getBasePairs()) {
			continue;
		} else {
			// Seed fits, check how many bases are possible around
			possibleBasePairs = std::min(std::min(helixSeed.size1()-i1+offset1, helixSeed.size2()-i2+offset2), helixConstraint.getMaxBasePairs())-seedHandler->getConstraint().getBasePairs();
		}

		// Initialuze variables
		curE = E_INF;
		curE_withED = E_INF;
		bestE = E_INF;
		bestE_withED = E_INF;
		bestL1 = 0;
		bestL2 = 0;

		// screen over all possible leading and trailing base pair combinations
		for (size_t leadingBP=possibleBasePairs+1;  leadingBP-- > 0;) {

			// check if boundaries are met
			if ((i1+leadingBP-offset1) >= helixSeed.size1() && (i2+leadingBP-offset2) >= helixSeed.size2()) {
				continue;
			}

			// If leading base pairs exist and helixE = E_INF -> skip to the next leadingBP
			if (leadingBP != 0 && E_isINF(getHelixE(i1-offset1,i2-offset2,leadingBP+1))) {
				continue;
			}

			// the start positions for the seed
			if (leadingBP != 0) {
				seedStart1 = i1 + getHelixLength1(i1 - offset1, i2 - offset2, leadingBP + 1) - 1;
				seedStart2 = i2 + getHelixLength2(i1 - offset1, i2 - offset2, leadingBP + 1) - 1;
			} else {
				seedStart1 = i1;
				seedStart2 = i2;
			}

			// check whether the right boundaries are broken
			if (seedStart1-offset1 >= helixSeed.size1() || seedStart2-offset2 >= helixSeed.size2()) {
				continue;
			}

			// If no seed is possible here, skip to next leading base pair number
			if (!(seedHandler->isSeedBound(seedStart1, seedStart2))) {
				continue;
			}

			// the end positions of the seed
			seedEnd1 = seedStart1+seedHandler->getSeedLength1(seedStart1, seedStart2)-1;
			seedEnd2 = seedStart2+seedHandler->getSeedLength2(seedStart1, seedStart2)-1;

			// Run over all trailing base pairs
			for (size_t trailingBP=possibleBasePairs+1 - leadingBP; trailingBP-- > 0;) {

				// check if boundaries are met
				if ((seedEnd1+trailingBP-offset1) >= helixSeed.size1() && (seedEnd2+trailingBP-offset2) >= helixSeed.size2()) {
					continue;
				}
				// If trailing base pairs exist and helixE = E_INF -> skip to the next leadingBP
				if (trailingBP != 0 && E_isINF(getHelixE(seedEnd1 - offset1, seedEnd2 - offset2, trailingBP + 1))) {
					continue;
				}

				// right boundary
				if (trailingBP != 0) {
					j1 = seedEnd1+getHelixLength1(seedEnd1-offset1,seedEnd2-offset2,trailingBP+1)-1;
					j2 = seedEnd2+getHelixLength2(seedEnd1-offset1,seedEnd2-offset2,trailingBP+1)-1;
				} else {
					j1 = seedEnd1;
					j2 = seedEnd2;
				}

				// check whether the right boundaries are broken
				if (j1-offset1 >= helixSeed.size1()
					|| j2-offset2 >= helixSeed.size2()) {
					continue;
				}

				// ensure that ED-values are within the boundaries (default 999)
				if (energy.getED1(i1, j1) > helixConstraint.getMaxED()
					|| energy.getED2(i2, j2) > helixConstraint.getMaxED()) {
					continue;
				}


				// energy without contributions
				curE = getHelixE(i1 - offset1, i2 - offset2, leadingBP + 1)
					   + seedHandler->getSeedE(seedStart1, seedStart2)
					   + getHelixE(seedEnd1 - offset1, seedEnd2 - offset2, trailingBP + 1);
				// energy value
				curE_withED = energy.getE(i1,j1,i2,j2, curE) + energy.getE_init();

				// If no ED-values are wanted, remove them
				if (!helixConstraint.evalFullE())
					curE_withED -= (energy.getED1(i1,j1) + energy.getED2(i2, j2));

				if ((curE_withED < bestE_withED || E_equal(curE_withED, bestE_withED))
					&& leadingBP + trailingBP+ seedHandler->getConstraint().getBasePairs() >= getConstraint().getMinBasePairs()) {

					bestE_withED = curE_withED;
					bestE = curE;
					// current lengths
					bestL1 = seedHandler->getSeedLength1(seedStart1, seedStart2);
					bestL2 = seedHandler->getSeedLength2(seedStart1, seedStart2);

					// Add leadingBP length contribution
					if (leadingBP != 0) {
						bestL1 += getHelixLength1(i1-offset1, i2-offset2, leadingBP+1)-1;
						bestL2 += getHelixLength2(i1-offset1, i2-offset2, leadingBP+1)-1;
					}
					// Add trailingBP length contribution
					if (trailingBP != 0) {
						bestL1 += getHelixLength1(seedEnd1-offset1, seedEnd2-offset2, trailingBP+1)-1;
						bestL2 += getHelixLength2(seedEnd1-offset1, seedEnd2-offset2, trailingBP+1)-1;
					}

				}
			} // trailingBP
		} // leadingBP

		// Ensures that the helixCount is only increased for the mfe helix.
		if (E_isNotINF(bestE_withED)) {
			// overwrite all helices with too high energy -> infeasible start interactions
			if (bestE_withED > helixConstraint.getMaxE()) {
				bestE = E_INF;
			} else {
				helixCountNotInf++;
			}
		}

		// set best energy and according lengths, or (INF, 0) if not valid
		helixSeed(i1-offset1, i2-offset2) = HelixSeedMatrix::value_type(bestE, E_isINF(bestE) ? 0: encodeHelixSeedLength(bestL1,bestL2));


	} // i2
	} // i1

#if INTARNA_MULITHREADING
#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) << "valid helices with seed = " << helixCountNotInf ; }

	return helixCountNotInf;
}

//////////////////////////////////////////////////////////////////////////

void
HelixHandlerUnpaired::
traceBackHelixSeed( Interaction & interaction
		, const size_t i1_
		, const size_t i2_)
{
	// check if we can abort since no seed fits into any allowed helix
	if (seedHandler->getConstraint().getBasePairs() > helixConstraint.getMaxBasePairs()) {
		return;
	}

	size_t i1 = i1_
	, i2 = i2_
	, seedStart1, seedEnd1
	, seedStart2, seedEnd2
	, j1, j2
	, curL1, curL2;

	bool traceNotFound = true;

	E_type curE = getHelixSeedE(i1_,i2_);

	size_t bestL1 = getHelixSeedLength1(i1_, i2_);
	size_t bestL2 = getHelixSeedLength2(i1_, i2_);
	// No traceback possible for current boundary
	if (E_isINF(curE)) {
		return;
	}

	// Calculate how many base pairs are possible allongside the seed.
	// Note: If seedHandler allows unpaired positions this check is not enough, check happens in loop
	size_t possibleBasePairs = std::min(std::min(helixSeed.size1()-i1 +offset1, helixSeed.size2()-i2+offset2), helixConstraint.getMaxBasePairs())-seedHandler->getConstraint().getBasePairs();

	// screen over all possible leading and trailing base pair combinations
	for (size_t leadingBP=0; traceNotFound
							 && leadingBP <= possibleBasePairs
							 && i1 + leadingBP-offset1 < helixSeed.size1()
							 && i2 + leadingBP-offset2 < helixSeed.size2(); leadingBP++) {

		// If leading base pairs exist and helixE = E_INF -> skip to the next leadingBP
		if (leadingBP != 0 && E_isINF(getHelixE(i1-offset1,i2-offset2,leadingBP+1))) {
			continue;
		}

		// the start positions for the seed
		if (leadingBP != 0) {
			seedStart1 = i1 + getHelixLength1(i1 - offset1, i2 - offset2, leadingBP + 1) - 1;
			seedStart2 = i2 + getHelixLength2(i1 - offset1, i2 - offset2, leadingBP + 1) - 1;
		} else {
			seedStart1 = i1;
			seedStart2 = i2;
		}
		// check whether the right boundaries are broken
		if (seedStart1-offset1 >= helixSeed.size1() || seedStart2-offset2 >= helixSeed.size2()) {
			continue;
		}

		// Check whether seed is possible for this starting position
		if (!(seedHandler->isSeedBound(seedStart1, seedStart2))) {
			continue;
		}

		// end positions of the seed
		seedEnd1 = seedStart1 + seedHandler->getSeedLength1(seedStart1, seedStart2) - 1;
		seedEnd2 = seedStart2 + seedHandler->getSeedLength2(seedStart1, seedStart2) - 1;

		// Trailing base pairs
		for (size_t trailingBP = 0; traceNotFound
									&& trailingBP <= possibleBasePairs - leadingBP
									&& (seedEnd1 + trailingBP - offset1) < helixSeed.size1()
									&& (seedEnd2 + trailingBP - offset2) < helixSeed.size2(); trailingBP++) {


			// If trailing base pairs exist and helixE = E_INF -> skip to the next leadingBP
			if (trailingBP != 0 && E_isINF(getHelixE(seedEnd1 - offset1, seedEnd2 - offset2, trailingBP + 1))) {
				continue;
			}

			// right boundary
			if (trailingBP != 0) {
				j1 = seedEnd1+getHelixLength1(seedEnd1-offset1,seedEnd2-offset2,trailingBP+1)-1;
				j2 = seedEnd2+getHelixLength2(seedEnd1-offset1,seedEnd2-offset2,trailingBP+1)-1;
			} else {
				j1 = seedEnd1;
				j2 = seedEnd2;
			}

			// check whether the right boundaries are broken
			if (j1-offset1 >= helixSeed.size1()
				|| j2-offset2 >= helixSeed.size2()) {
				continue;
			}

			// check whether trace is found
			if (E_equal(curE, getHelixE(i1 - offset1, i2 - offset2, leadingBP + 1)
							  + seedHandler->getSeedE(seedStart1, seedStart2)
							  + getHelixE(seedEnd1 - offset1, seedEnd2 - offset2, trailingBP + 1))) {

				// calculate the current length
				curL1 = seedHandler->getSeedLength1(seedStart1, seedStart2);
				curL2 = seedHandler->getSeedLength2(seedStart1, seedStart2);

				// Add leadingBP length contribution
				if (leadingBP != 0) {
					curL1 += getHelixLength1(i1-offset1, i2-offset2, leadingBP+1)-1;
					curL2 += getHelixLength2(i1-offset1, i2-offset2, leadingBP+1)-1;
				}
				// Add trailingBP length contribution
				if (trailingBP != 0) {
					curL1 += getHelixLength1(seedEnd1-offset1, seedEnd2-offset2, trailingBP+1)-1;
					curL2 += getHelixLength2(seedEnd1-offset1, seedEnd2-offset2, trailingBP+1)-1;
				}

				// ensure that the lengths are correct
				if ( curL1 == bestL1  && curL2 == bestL2 ) {
					// Trace the leading part if existing
					if (leadingBP != 0) {
						traceBackHelix(interaction, i1 - offset1, i2 - offset2, leadingBP + 1);
						interaction.basePairs.push_back(energy.getBasePair(seedStart1, seedStart2));
					}
					// Trace the seed
					seedHandler->traceBackSeed(interaction, seedStart1, seedStart2);
					// Trace the trailing part if existing
					if (trailingBP != 0) {
						interaction.basePairs.push_back(energy.getBasePair(seedEnd1, seedEnd2));
						traceBackHelix(interaction, seedEnd1 - offset1, seedEnd2 - offset2, trailingBP + 1);
					}
					traceNotFound = false;
				}
			}
		} // trailing
	} // leading
	assert(!traceNotFound);

} // traceback

//////////////////////////////////////////////////////////////////////////


} // namespace


