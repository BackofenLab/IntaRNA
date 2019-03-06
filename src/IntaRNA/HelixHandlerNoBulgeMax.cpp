#include "IntaRNA/HelixHandlerNoBulgeMax.h"

namespace IntaRNA {

/////////////////////////////////////////////////////////////////////////////

size_t
HelixHandlerNoBulgeMax::
fillHelix(const size_t i1min, const size_t i1max, const size_t i2min, const size_t i2max)
{

#if INTARNA_IN_DEBUG_MODE
	if ( i1min > i1max ) throw std::runtime_error("HelixHandlerNoBulgeMax::fillHelix: i1min("+toString(i1min)+") > i1max("+toString(i1max)+")");
	if ( i2min > i2max ) throw std::runtime_error("HelixHandlerNoBulgeMax::fillHelix: i2min("+toString(i2min)+") > i2max("+toString(i2max)+")");
	if ( i1max > energy.size1() ) throw std::runtime_error("HelixHandlerNoBulgeMax::fillHelix: i1max("+toString(i1max)+") > energy.size1("+toString(energy.size1())+")");
	if ( i2max > energy.size2() ) throw std::runtime_error("HelixHandlerNoBulgeMax::fillHelix: i2max("+toString(i2max)+") > energy.size2("+toString(energy.size2())+")");
	if ( helixConstraint.getMinBasePairs() > helixConstraint.getMaxBasePairs() )
		throw std::runtime_error("HelixHandlerNoBulgeMax::fillHelix: bpMin("+toString(helixConstraint.getMinBasePairs()) +") > bpMax("+toString(helixConstraint.getMaxBasePairs())+")");
#endif

	// measure timing
	TIMED_FUNC_IF(timerObj,VLOG_IS_ON(9));

	// TODO replace with hash
	// reset storage of maximal helix information
	helix.clear();

	const size_t minBP = helixConstraint.getMinBasePairs();
	const E_type maxE = helixConstraint.getMaxE();

	//! container to store stacking energies of canonical helix base pairs
	typedef std::vector<E_type> StackingEnergyList;
	// left-stacking energies of all helix base pairs
	// (excludes right-most bp with E=E_init)
	// data structure is used as a ring list to avoid energy recomputation
	StackingEnergyList bpE(helixConstraint.getMaxBasePairs()-1,(E_type)0);
	// end index of current seed information within bpE ringlist
	size_t bpE_end = 0;

	const size_t i1maxBound=std::min( i1max+1, energy.size1() );
	const size_t i2maxBound=std::min( i2max+1, energy.size2() );
	size_t i1,i2,maxOffset,o,curHelixLength,curBestLength;
	bool validHelix = false;
	E_type curHelixE = E_INF, curBestE = E_INF;
	// for all start indices shifts in seq1 and seq2
	for (size_t i1start = i1min; i1start < i1maxBound; i1start++) {
	for (size_t i2start = i2min; (i2start==i2min) || (i1start==i1min && i2start < i2maxBound); i2start++) {

		// shift to right-most bp parallel to (i1start,i2start)
		maxOffset = std::min(i1maxBound-i1start, i2maxBound-i2start);
		curHelixE = E_INF;
		curHelixLength = 0;

		// iterate over all left end offsets of a helix
		for( o = maxOffset; o-->0;) {

			// putative base pair indices to check
			i1 = i1start+o;
			i2 = i2start+o;

			// check if valid base pair
			if( energy.isAccessible1( i1 )
				&& energy.isAccessible2( i2 )
				&& energy.areComplementary( i1, i2 ))
			{
				// start new canonical helix information
				if (E_isINF(curHelixE)) {
					curHelixE = E_type(0);
					curHelixLength = 1;
				}
				// update helix information
				else {
					// extended already to maxBP
					if (curHelixLength > bpE.size()) {
						// correct curHelixE
						curHelixE -= bpE[(bpE_end)%bpE.size()];
					} else {
						// increase helix length
						curHelixLength++;
					}
					// compute and store new stacking energy
					bpE[(bpE_end)%bpE.size()] = energy.getE_interLeft(i1,i1+1,i2,i2+1);
					// update overall helix energy
					curHelixE += bpE[(bpE_end)%bpE.size()];

					// shift storage information
					bpE_end++;
				}

				// find best helix length for current left boundary
				curBestE = E_INF;
				E_type curBestEfullDelta = E_type(0);
				curBestLength = 0;
				E_type curE = curHelixE, curEFullDelta = E_type(0); // init with loop energies only
				for (size_t l = curHelixLength; l>=minBP; l--) {
					validHelix = true;
					// ensure that ED-values are within the boundaries
					validHelix = validHelix
							&& energy.getED1(i1, i1+l-1) <= helixConstraint.getMaxED()
							&& energy.getED2(i2, i2+l-1) <= helixConstraint.getMaxED();
					// check how to evaluate energy
					if (helixConstraint.evalFullE()) {
						// get additional energy terms
						curEFullDelta = energy.getE(i1,i1+l-1,i2,i2+l-1,E_type(0)) + energy.getE_init();
					}
					// check maxE
					validHelix = validHelix && (curE+curEFullDelta) < maxE;
					// update best helix information
					if ( validHelix && (curE+curEFullDelta) < (curBestE+curBestEfullDelta) ) {
						curBestE = curE;
						curBestEfullDelta = curEFullDelta;
						curBestLength = l;
					}
					// update curE
					curE -= bpE[(bpE_end+1-l)%bpE.size()];
				}

				// store helix information
				if (E_isNotINF(curBestE)) {
					// store helix information
					helix[BP(i1,i2)] = HelixData(curBestE,curBestLength);
				}
			} else {
				// reset data structures
				for(auto e=bpE.begin(); e!=bpE.end(); e++) { *e = (E_type)0; }
				bpE_end = 0; // set to start position
				curHelixE = E_INF; // reset helix energy
				curHelixLength = 0; // reset helix length
			}

		}

	} // i2start
	} // i2start

#if INTARNA_MULITHREADING
#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) << "valid helices = " << helix.size() ; }

	return helix.size();
}

////////////////////////////////////////////////////////////////////////////

void
HelixHandlerNoBulgeMax::
traceBackHelix(Interaction &interaction
		, const size_t i1
		, const size_t i2
)
{
#if INTARNA_IN_DEBUG_MODE
	if ( E_isINF( getHelixE(i1,i2) ) ) throw std::runtime_error("HelixHandlerNoBulgeMax::traceBackHelix(i1="+toString(i1)+",i2="+toString(i2)+") no helix known (E_INF)");
#endif
	// get number of base pairs for best helix
	const size_t bp = getHelixLength1(i1,i2);

	// trace each helix base pair (excluding left and right most)
	for (size_t curBP = 1; curBP < bp-1; curBP++) {
		interaction.basePairs.push_back(energy.getBasePair(i1+curBP, i2+curBP));
	}
}

/////////////////////////////////////////////////////////////////////////////

size_t
HelixHandlerNoBulgeMax::
fillHelixSeed(const size_t i1min, const size_t i1max, const size_t i2min, const size_t i2max)
{
#if INTARNA_IN_DEBUG_MODE
	if ( i1min > i1max ) throw std::runtime_error("HelixHandlerNoBulgeMax::fillHelixSeed: i1min("+toString(i1min)+") > i1max("+toString(i1max)+")");
	if ( i2min > i2max ) throw std::runtime_error("HelixHandlerNoBulgeMax::fillHelixSeed: i2min("+toString(i2min)+") > i2max("+toString(i2max)+")");
	if ( i1max > energy.size1() ) throw std::runtime_error("HelixHandlerNoBulgeMax::fillHelixSeed: i1max("+toString(i1max)+") > energy.size1("+toString(energy.size1())+")");
	if ( i2max > energy.size2() ) throw std::runtime_error("HelixHandlerNoBulgeMax::fillHelixSeed: i2max("+toString(i2max)+") > energy.size2("+toString(energy.size2())+")");
	if ( helixConstraint.getMinBasePairs() > helixConstraint.getMaxBasePairs() )
		throw std::runtime_error("HelixHandlerNoBulgeMax::fillHelixSeed: bpMin("+toString(helixConstraint.getMinBasePairs()) +") > bpMax("+toString(helixConstraint.getMaxBasePairs())+")");
	if ( seedHandler == NULL ) throw std::runtime_error("HelixHandlerNoBulgeMax::fillHelixSeed() no SeedHandler available");
#endif

	// measure timing
	TIMED_FUNC_IF(timerObj,VLOG_IS_ON(9));

	// clear data structure
	helixSeed.clear();

	// temporary variables
	size_t i1, i2,j1,j2, seedStart1, seedStart2, seedEnd1, seedEnd2, bestL1, bestL2, possibleBasePairs;

	E_type bestE, curE, bestEfullDelta, curEfullDelta;

	// fill for all start indices
	// in decreasing index order
	for (i1=i1max+1; i1-- > i1min;) {
	for (i2=i2max+1; i2-- > i2min;) {

		// skip non-complementary left helix boundaries
		if (!(energy.isAccessible1(i1)
			&& energy.isAccessible2(i2)
			&& energy.areComplementary(i1, i2)))
		{
			continue; // go to next helixSeedE index
		}

		// Check if a seed can fit given the left boundaries
		// Note: If seedHandler allows unpaired positions this check is not enough, check happens in loop
		if (std::min(i1max+1-i1, i2max+1-i2) < seedHandler->getConstraint().getBasePairs()) {
			continue;
		} else {
			// Seed fits, check how many bases are possible around
			possibleBasePairs = std::min(std::min(i1max+1-i1, i2max+1-i2), helixConstraint.getMaxBasePairs())-seedHandler->getConstraint().getBasePairs();
		}

		// Initialuze variables
		curEfullDelta = E_type(0);
		curE = E_INF;
		bestEfullDelta = E_type(0);
		bestE = E_INF;
		bestL1 = 0;
		bestL2 = 0;

		// stores energy of leading base pairs to avoid recomputation
		E_type leadingE = E_type(0);
		size_t leadingL = 0;


		// screen over all possible leading and trailing base pair combinations
		for (size_t leadingBP=0; leadingBP <= possibleBasePairs
								 && (i1+leadingBP) <= i1max
								 && (i2+leadingBP) <= i2max; leadingBP++)
		{
			// the start positions for the seed
			seedStart1 = i1+leadingBP;
			seedStart2 = i2+leadingBP;

			// check if leading based pairs are possible, otherwise stop computation
			if (!(energy.isAccessible1(seedStart1)
				&& energy.isAccessible2(seedStart2)
				&& energy.areComplementary(seedStart1,seedStart2)))
			{
				break;
			}

			// If no seed is possible here, skip to next leading base pair number
			if (E_isINF(seedHandler->getSeedE(seedStart1, seedStart2))) {
				continue;
			}

			// the end positions of the seed
			seedEnd1 = seedStart1+seedHandler->getSeedLength1(seedStart1, seedStart2)-1;
			seedEnd2 = seedStart2+seedHandler->getSeedLength2(seedStart1, seedStart2)-1;

			// update leadingE to cover all (missing) leading base pairs
			for (; leadingL < leadingBP; leadingL++) {
				leadingE += energy.getE_interLeft(i1+leadingL, i1+leadingL+1, i2+leadingL, i2+leadingL+1);
			}

			E_type trailingE = E_type(0);
			size_t trailingL = 0;

			// Run over all trailing base pairs
			for (size_t trailingBP=0; trailingBP <= possibleBasePairs - leadingBP
									  && (seedEnd1+trailingBP) <= i2max
									  && (seedEnd2+trailingBP) <= i2max; trailingBP++) {

				j1 = seedEnd1+trailingBP;
				j2 = seedEnd2+trailingBP;

				// check if trailing based pairs are possible, otherwise extension computation
				if (trailingBP>0
					&& !(energy.isAccessible1(j1)
						&& energy.isAccessible2(j2)
						&& energy.areComplementary(j1,j2)))
				{
					break;
				}

				// ensure that ED-values are within the boundaries (default 999)
				if (energy.getED1(i1, j1) > helixConstraint.getMaxED()
					|| energy.getED2(i2, j2) > helixConstraint.getMaxED()) {
					continue;
				}

				// update trailingE to cover all (missing) trailing base pairs
				for (; trailingL < trailingBP; trailingL++ ) {
					trailingE += energy.getE_interLeft(seedEnd1+trailingL-1, seedEnd1+trailingL, seedEnd2+trailingL-1, seedEnd2+trailingL);
				}

				// energy value without dangling ends and E_init
				curE = leadingE +
					   seedHandler->getSeedE(seedStart1, seedStart2) +
					   trailingE;

				// check how to evaluate energy terms
				if (helixConstraint.evalFullE()) {
					// get additional energy terms
					curEfullDelta = energy.getE(i1, j1, i2, j2, 0.0) + energy.getE_init();
				}

				// check maxE constraint
				// Check whether new combination is better
				// check minBP constraint
				if ( (curE+curEfullDelta) < helixConstraint.getMaxE()
					&& (curE+curEfullDelta) < (bestE+bestEfullDelta)
					&& (leadingBP + trailingBP + seedHandler->getConstraint().getBasePairs()) >= getConstraint().getMinBasePairs())
				{
					bestE = curE;
					bestEfullDelta = curEfullDelta;
					bestL1 = leadingBP + seedHandler->getSeedLength1(seedStart1,seedStart2) + trailingBP;
					bestL2 = leadingBP + seedHandler->getSeedLength2(seedStart1,seedStart2) + trailingBP;
				}

			} // trailingBP
		} // leadingBP

		// check if valid helix found
		if (E_isNotINF( bestE )) {
			// store best (mfe) helix for left boundary i1, i2
			helixSeed[BP(i1,i2)] = HelixData(bestE, encodeHelixSeedLength(bestL1,bestL2));
		}


	} // i2
	} // i1

#if INTARNA_MULITHREADING
#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) << "valid helices with seed = " << helixSeed.size() ; }


	return helixSeed.size();
}

//////////////////////////////////////////////////////////////////////////

void
HelixHandlerNoBulgeMax::
traceBackHelixSeed( Interaction & interaction
		, const size_t i1
		, const size_t i2)
{
#if INTARNA_IN_DEBUG_MODE
	if ( seedHandler == NULL ) throw std::runtime_error("HelixHandlerNoBulgeMax::traceBackHelixSeed() no SeedHandler available");
//	if ( helixSeed.find(BP(i1,i2)) == helixSeed.end() ) throw std::runtime_error("HelixHandlerNoBulgeMax::traceBackHelixSeed() no helix available for left bp "+toString(i1)+","+toString(i2));
#endif
	// check if something to trace
	if ( helixSeed.find(BP(i1,i2)) == helixSeed.end() ) {
		return;
	}

	size_t seedStart1, seedEnd1
	     , seedStart2, seedEnd2;

	bool traceNotFound = true;

	// get energy to be traced
	E_type curE = helixSeed.at(BP(i1,i2)).first;

	const size_t bestL1 = getHelixSeedLength1(i1,i2);
	const size_t bestL2 = getHelixSeedLength2(i1,i2);

	const size_t j1 = i1 + bestL1 -1;
	const size_t j2 = i2 + bestL2 -1;

	// No traceback possible for current boundary
	if (E_isINF(curE)) {
		return;
	}

	// Calculate how many base pairs are possible allongside the seed.
	// Note: If seedHandler allows unpaired positions this check is not enough, check happens in loop
	size_t possibleBasePairs = std::min(std::min(j1+1-i1, j2+1-i2), helixConstraint.getMaxBasePairs())-seedHandler->getConstraint().getBasePairs();

	// stores energy of leading base pairs to avoid recomputation
	E_type leadingE = E_type(0);
	size_t leadingL = 0;

	// screen over all possible leading and trailing base pair combinations
	for (size_t leadingBP=0; traceNotFound
							 && leadingBP <= possibleBasePairs
							 && (i1 + leadingBP) <= j1
							 && (i2 + leadingBP) <= j2; leadingBP++) {

		// check if leading based pairs are possible, otherwise stop computation
		if (!energy.areComplementary(i1+leadingBP,i2+leadingBP))
		{
			break;
		}

		// start positions of the seed
		seedStart1 = i1 + leadingBP;
		seedStart2 = i2 + leadingBP;

		// Check whether seed is possible for this starting position
		if (E_isINF(seedHandler->getSeedE(seedStart1, seedStart2))) {
			continue;
		}

		// end positions of the seed
		seedEnd1 = seedStart1+seedHandler->getSeedLength1(seedStart1,seedStart2)-1;
		seedEnd2 = seedStart2+seedHandler->getSeedLength2(seedStart1,seedStart2)-1;

		// check if seed is within helix boundary
		// check if trailing gap can be filled with canonical base pairs
		if ( seedEnd1>j1 || seedEnd2>j2 || (j1-seedEnd1) != (j2-seedEnd2) ) {
			continue;
		}

		// update leadingE to cover all (missing) leading base pairs
		for (; leadingL < leadingBP; leadingL++) {
			leadingE += energy.getE_interLeft(i1+leadingL, i1+leadingL+1, i2+leadingL, i2+leadingL+1);
		}

		E_type trailingE = E_type(0);
		const size_t trailingBP = (j1-seedEnd1);
		// update trailingE to cover all (missing) trailing base pairs
		for (size_t trailingL = 0; trailingL < trailingBP; trailingL++) {
			// check if trailing based pairs are possible, otherwise stop computation
			if (!energy.areComplementary(seedEnd1+trailingL,seedEnd2-trailingL))
			{
				break;
			}
			// update energy
			trailingE += energy.getE_interLeft(seedEnd1+trailingL-1, seedEnd1+trailingL, seedEnd2+trailingL-1, seedEnd2+trailingL);
		}



		if (E_equal(curE, leadingE
						  + seedHandler->getSeedE(seedStart1, seedStart2)
						  + trailingE))
		{
			// length check
			if (bestL1 == (leadingBP + seedHandler->getSeedLength1(seedStart1,seedStart2) + trailingBP)
				&& bestL2 == (leadingBP + seedHandler->getSeedLength2(seedStart1,seedStart2) + trailingBP)) {

				// Trace the first part if existing
				for (size_t l = 1; l<=leadingBP; l++) {
					interaction.basePairs.push_back(energy.getBasePair(i1+l, i2+l));
				}
				// Trace the seed
				seedHandler->traceBackSeed(interaction,seedStart1,seedStart2);
				// store seed information
				interaction.setSeedRange(
								energy.getBasePair(seedStart1,seedStart2),
								energy.getBasePair(seedEnd1,seedEnd2),
								energy.getE(seedStart1,seedEnd1,seedStart2,seedEnd2,seedHandler->getSeedE(seedStart1,seedStart2))+energy.getE_init());
				// Trace the last part if existing
				for (size_t l = 0; l<trailingBP; l++) {
					interaction.basePairs.push_back(energy.getBasePair(seedEnd1+l, seedEnd2+l));
				}
				traceNotFound = false;
			}
		} // energy check

	} // leading bp
	assert(!traceNotFound);

} // traceback

} // namespace
