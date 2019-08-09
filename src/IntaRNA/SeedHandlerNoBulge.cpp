
#include "IntaRNA/SeedHandlerNoBulge.h"

namespace IntaRNA {

//////////////////////////////////////////////////////////////////////////

size_t
SeedHandlerNoBulge::
fillSeed( const size_t i1min, const size_t i1max, const size_t i2min, const size_t i2max)
{

#if INTARNA_IN_DEBUG_MODE
	if ( i1min > i1max ) throw std::runtime_error("SeedHandlerNoBulge::fillSeed: i1min("+toString(i1min)+") > i1max("+toString(i1max)+")");
	if ( i2min > i2max ) throw std::runtime_error("SeedHandlerNoBulge::fillSeed: i2min("+toString(i2min)+") > i2max("+toString(i2max)+")");
	if ( i1max > energy.size1() ) throw std::runtime_error("SeedHandlerNoBulge::fillSeed: i1max("+toString(i1max)+") > energy.size1("+toString(energy.size1())+")");
	if ( i2max > energy.size2() ) throw std::runtime_error("SeedHandlerNoBulge::fillSeed: i2max("+toString(i2max)+") > energy.size2("+toString(energy.size2())+")");
#endif
	// measure timing
	TIMED_FUNC_IF(timerObj,VLOG_IS_ON(9));

	// reset data
	seedForLeftEnd.clear();

	// temporary data structures
	const size_t seedBP = seedConstraint.getBasePairs();
	const size_t i1maxBound=std::min( i1max+1, energy.size1() );
	const size_t i2maxBound=std::min( i2max+1, energy.size2() );

	// left-stacking energies of all seed base pairs
	// (excludes right-most bp with E=E_init)
	// data structure is used as a ring list to avoid energy recomputation
	StackingEnergyList bpE(seedBP-1,(E_type)0);
	// start index of current seed information within bpE ringlist
	size_t bpE_start = 0;


	// for all start indices shifts in seq1 and seq2
	for (size_t i1start = i1min; i1start < i1maxBound; i1start++) {
	for (size_t i2start = i2min; (i2start==i2min) || (i1start==i1min && i2start < i2maxBound); i2start++) {


		// iterate over all right end offsets of a seed
		for( size_t o = seedBP-1; (o+i1start)<i1maxBound && (o+i2start)<i2maxBound;) {

			// get right end of seed
			size_t j1 = i1start + o;
			size_t j2 = i2start + o;

			// check if right end is a valid seed base pair
			if ( !isFeasibleSeedBasePair(j1,j2) )
			{
				// jump to next possible right seed end
				o += seedBP;
				// skip remaining processing of this right seed end
				continue;
			} // else right end is complementary

			// reset energy data structure
			for(auto e=bpE.begin(); e!=bpE.end(); e++) { *e = (E_type)0; }
			bpE_start = bpE.size(); // set to invalid position

			// extend seed to the left as far as possible
			size_t seedLength=1; // so far only right-most bp
			for (seedLength = 1; seedLength<seedBP; seedLength++) {
				// get current left-most seed bp
				size_t i1=j1-seedLength, i2=j2-seedLength;
				// check if not complementary etc.
				if ( !isFeasibleSeedBasePair(i1,i2) )
				{
					break; // stop left extension
				} // else can be extended to the left to (i1,i2)
				// update bpE (shift seed start to the left)
				bpE_start = bpE.size()-seedLength;
				bpE[bpE_start] = energy.getE_interLeft(i1,i1+1,i2,i2+1);
			}

			// extend to the right if needed
			if (bpE_start != 0) {
				for ( ; seedLength<seedBP && (j1+1)<i1maxBound && (j2+1)<i2maxBound; seedLength++) {
					// check if right extension possible
					if ( isFeasibleSeedBasePair(j1+1,j2+1) ) {
						// get left stacking energy
						bpE[(bpE_start+seedLength-1)%bpE.size()] = energy.getE_interLeft(j1,j1+1,j2,j2+1);
						// update right boundary information
						j1++;
						j2++;
						o++;
					} else {
						break; // stop right extension
					}
				}
				// check if full seed extension was NOT possible
				if (seedLength < seedBP) {
					// jump to next possible right seed end
					o += seedBP + 1;
					// skip further processing of this right seed end
					continue;
				}
			}

			// store current seed
			storeSeed( j1, j2, bpE );

			// slide seed along
			j1++; j2++; o++;
			while( j1<i1maxBound && j2<i2maxBound && isFeasibleSeedBasePair(j1,j2) ) {
				// update energy information of last left-stacking bp (overwrite old start pos)
				bpE[(bpE_start)%bpE.size()] = energy.getE_interLeft(j1-1,j1,j2-1,j2);
				// shift start index in ring list
				bpE_start++;
				// store seed information
				storeSeed( j1, j2, bpE );
				// update next right end possibility
				j1++; j2++; o++;
			}

			// jump to next possible right seed end
			o += seedBP;

		} // for( o )

	} // for( i2start )
	} // for( i1start )


#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) <<"valid seeds = "<<seedForLeftEnd.size(); }

	// get final number of identified seeds
	return seedForLeftEnd.size();
}

//////////////////////////////////////////////////////////////////////////

bool
SeedHandlerNoBulge::
updateToNextSeed( size_t & i1_out, size_t & i2_out
		, const size_t i1min, const size_t i1max
		, const size_t i2min, const size_t i2max
		) const
{
	// ensure we do have any seed
	if (seedForLeftEnd.empty()) {
		return false;
	}

	size_t i1 = i1_out, i2 = i2_out;
	// find true max value
	const size_t i1maxVal = std::min(energy.size1()-1,i1max)
				, i2maxVal = std::min(energy.size2()-1,i2max);

	// find current seed
	auto curSeedData = seedForLeftEnd.find( Interaction::BasePair(i1,i2) );
	// check if we have to provide first seed (out of bound or no seed start)
	if (curSeedData == seedForLeftEnd.end()) {
		curSeedData = seedForLeftEnd.begin();
	} else {
		// go to next seed
		curSeedData++;
	}

	// find next seed within range
	while (curSeedData != seedForLeftEnd.end()
			&& (curSeedData->first.first < i1min
				|| curSeedData->first.first > i1maxVal
				|| curSeedData->first.second < i2min
				|| curSeedData->first.second > i2maxVal
			))
	{
		curSeedData++;
	}
	// ensure we have a valid new seed
	if (curSeedData != seedForLeftEnd.end()) {
		// copy data
		i1_out = curSeedData->first.first;
		i2_out = curSeedData->first.second;
		return true;
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////

} // namespace
