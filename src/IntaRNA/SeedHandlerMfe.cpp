
#include "IntaRNA/SeedHandlerMfe.h"

namespace IntaRNA {

//////////////////////////////////////////////////////////////////////////

size_t
SeedHandlerMfe::
fillSeed( const size_t i1min, const size_t i1max, const size_t i2min, const size_t i2max)
{

#if INTARNA_IN_DEBUG_MODE
	if ( i1min > i1max ) throw std::runtime_error("SeedHandlerMfe::fillSeed: i1min("+toString(i1min)+") > i1max("+toString(i1max)+")");
	if ( i2min > i2max ) throw std::runtime_error("SeedHandlerMfe::fillSeed: i2min("+toString(i2min)+") > i2max("+toString(i2max)+")");
	if ( i1max > energy.size1() ) throw std::runtime_error("SeedHandlerMfe::fillSeed: i1max("+toString(i1max)+") > energy.size1("+toString(energy.size1())+")");
	if ( i2max > energy.size2() ) throw std::runtime_error("SeedHandlerMfe::fillSeed: i2max("+toString(i2max)+") > energy.size2("+toString(energy.size2())+")");
#endif
	// measure timing
	TIMED_FUNC_IF(timerObj,VLOG_IS_ON(9));

	// resize matrizes
	seed.resize( i1max-i1min+1, i2max-i2min+1 );
	seedE_rec.resize( SeedIndex({{ // setup ring-list data for seed computation
					  (SeedRecMatrix::index)(seed.size1())
					, (SeedRecMatrix::index)(seed.size2())
					, (SeedRecMatrix::index)(seedConstraint.getBasePairs()+1-2) // +1 for size and -2 to encode at least 2 bps or more
					, (SeedRecMatrix::index)(seedConstraint.getMaxUnpaired1()+1) // +1 for size
					, (SeedRecMatrix::index)(seedConstraint.getMaxUnpaired2()+1) // +1 for size
				}}));

	// store index offset due to restricted matrix size generation
	offset1 = i1min;
	offset2 = i2min;

	// temporary variables
	size_t i1, i2, bpIn, u1, u2, j1, j2, u1p, u2p, k1,k2, u1best, u2best;
	E_type curE, bestE;

	// determine whether or not lonely base pairs are allowed or if we have to
	// ensure a stacking to the right of the left boundary (i1,i2)
	const size_t noLpShift = seedConstraint.isLpAllowed() ? 0 : 1;
	E_type iStackE = E_type(0);

	size_t seedCountNotInf = 0, seedCount = 0;

	// fill for all start indices
	// in decreasing index order
	for (i1=i1max+1; i1-- > i1min;) {
	for (i2=i2max+1; i2-- > i2min;) {

		// count seed possibility
		seedCount++;

		// init according to no seed interaction
		seed(i1-offset1,i2-offset2) = SeedMatrix::value_type( E_INF, 0 );

		// skip non-complementary or infeasible seed base pair
		if (!isFeasibleSeedBasePair(i1,i2)) {
			continue; // go to next seedE index
		}

		// for feasible number of base pairs (bp+1) in increasing order
		// bp=0 encodes 2 base pairs
		for (bpIn=0; bpIn<seedE_rec.shape()[2] && (i1+bpIn+1-offset1)<seed.size1() && (i2+bpIn+1-offset2)<seed.size2(); bpIn++) {

			bool validLeftEnd = true;

			// if no LP : check for direct right-stacking of i
			if (noLpShift > 0) {
				// check if feasible extension
				if (isFeasibleSeedBasePair(i1+1,i2+1)) {
					// get stacking energy
					iStackE = energy.getE_interLeft(i1,i1+1,i2,i2+1);
				} else {
					// no valid noLP extension possible
					validLeftEnd = false;
				}
			}

			// for feasible unpaired in seq1 in increasing order
			for (u1=0; u1<seedE_rec.shape()[3] && (i1+bpIn+1+u1-offset1) < seed.size1(); u1++) {

				// get right seed boundaries
				// check if this index range is to be considered for seed search

			// for feasible unpaired in seq2 in increasing order
			for (u2=0; u2<seedE_rec.shape()[4] && (u1+u2)<=seedConstraint.getMaxUnpairedOverall() && (i2+bpIn+1+u2-offset2) < seed.size2(); u2++) {

				// get right seed boundaries
				j1 = i1+bpIn+1+u1;
				j2 = i2+bpIn+1+u2;

				// init current seed energy
				curE = E_INF;

				// check if this index range is to be considered for seed search
				// check if boundaries are complementary
				if (validLeftEnd
					&& isFeasibleSeedBasePair(j1,j2,true)
					&& (noLpShift==0 || isFeasibleSeedBasePair(j1-1,j2-1)))
				{

					// base case: only left and right base pair present
					if (bpIn==0) {
						// if lonely bps are allowed or no bulge
						if (noLpShift == 0 || (u1==0 && u2==0)) {
							// energy for stacking/bulge/interior depending on u1/u2
							curE = energy.getE_interLeft(i1,j1,i2,j2);
						}

					} else {

						// explicitly check direct stacking extension in noLP mode
						if (noLpShift > 0 && E_isNotINF( getSeedE( i1+1, i2+1, bpIn-1, u1, u2 ) )) {
							curE = std::min( curE, iStackE + getSeedE( i1+1, i2+1, bpIn-1, u1, u2 ) );
						}

						// if enough interior base pairs left
						if (bpIn >= 1+noLpShift) {
							// split seed recursively into all possible leading interior loops
							// i1 .. i1+u1p+1 .. j1
							// i2 .. i2+u2p+1 .. j2
							for (u1p=1+std::min(u1,energy.getMaxInternalLoopSize1()); u1p-- > 0;) {
							for (u2p=1+std::min(u2,energy.getMaxInternalLoopSize2()); u2p-- > 0;) {

								// skip stacked extension for noLP since already covered above
								if (u1p+u2p < noLpShift) { continue; }

								k1 = i1+u1p+1+noLpShift;
								k2 = i2+u2p+1+noLpShift;
								// check if split pair is complementary
								// and recursed entry is < E_INF
								if (! (isFeasibleSeedBasePair(k1,k2) && E_isNotINF( getSeedE( k1, k2, bpIn-1-noLpShift, u1-u1p, u2-u2p ) ) ) ) {
									continue; // not complementary -> skip
								}

								// update mfe for split at k1,k2
								curE = std::min( curE,
										iStackE
										+ energy.getE_interLeft(i1+noLpShift,k1,i2+noLpShift,k2)
										+ getSeedE( k1, k2, bpIn-1-noLpShift, u1-u1p, u2-u2p )
										);
							} // u2p
							} // u1p
						}
					} // more than two base pairs

				} // (j1,j2) complementary

				// store seed energy
				setSeedE( i1, i2, bpIn, u1, u2, curE );
			} // u2
			} // u1

			// check if full base pair number reached
			// check if feasible left seed boundaries
			if ( bpIn+1==seedE_rec.shape()[2] && isFeasibleSeedBasePair(i1,i2,true) ) {

				// find best unpaired combination in seed for i1,i2,bp
				u1best = 0;
				u2best = 0;
				bestE = E_INF;

				// for feasible unpaired in seq1 in increasing order
				for (u1=0; u1<seedE_rec.shape()[3] && (i1+bpIn+1+u1-offset1) < seed.size1(); u1++) {
				// for feasible unpaired in seq2 in increasing order
				for (u2=0; u2<seedE_rec.shape()[4] && (u1+u2)<=seedConstraint.getMaxUnpairedOverall() && (i2+bpIn+1+u2-offset2) < seed.size2(); u2++) {

					// get right seed boundaries
					j1 = i1+bpIn+1+u1;
					j2 = i2+bpIn+1+u2;

					// skip if ED boundary exceeded
					if (energy.getED1(i1,j1) >= seedConstraint.getMaxED()
						|| energy.getED2(i2,j2) >= seedConstraint.getMaxED() )
					{
						continue;
					}

					// get overall interaction energy
					E_type curEhyb = getSeedE( i1, i2, bpIn, u1, u2 ) + energy.getE_init();
					curE = energy.getE( i1, j1, i2, j2, curEhyb );

					// check hybrid energy bound including Einit
					// check if better than what is known so far
					if ( curEhyb <= seedConstraint.getMaxEhybrid()
						&& curE <= seedConstraint.getMaxE()
						&& curE < bestE )
					{
						bestE = curE;
						u1best = u1;
						u2best = u2;
					}
				} // u2
				} // u1

				// reduce bestE to hybridization energy only (init+loops)
				if (E_isNotINF( bestE )) {
					// get seed's hybridization loop energies only
					bestE = getSeedE( i1, i2, bpIn, u1best, u2best );
					// count true seed
					seedCountNotInf++;
				}

				// store best (mfe) seed for all u1/u2
				seed(i1-offset1,i2-offset2) = SeedMatrix::value_type( bestE
						, E_isINF(bestE)?0:encodeSeedLength(bpIn+2+u1best,bpIn+2+u2best) );

			} // store best seed

		} // bp
	} // i2
	} // i1

#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) <<"valid seeds = "<<seedCountNotInf <<" ("<<(seedCountNotInf/seedCount)<<"% of start index combinations)"; }

	// return final number of valid seeds
	return seedCountNotInf;
}

//////////////////////////////////////////////////////////////////////////

void
SeedHandlerMfe::
traceBackSeed( Interaction & interaction
		, const size_t i1_
		, const size_t i2_
		, const size_t bpInbetween
		, const size_t u1_
		, const size_t u2_
		) const
{

	// get boundaries
	size_t 	  i1 = i1_
			, i2 = i2_
			, j1 = i1+bpInbetween+1+u1_
			, j2 = i2+bpInbetween+1+u2_
			, u1max = u1_
			, u2max = u2_
			, uMax = u1_+u2_
			, u1, u2
			, k1, k2
	   ;

	// get energy of provided seed
	E_type curE = getSeedE(i1_,i2_,bpInbetween,u1_,u2_);

	// determine whether or not lonely base pairs are allowed or if we have to
	// ensure a stacking to the right of the left boundary (i1,i2)
	const size_t noLpShift = seedConstraint.isLpAllowed() ? 0 : 1;
	E_type iStackE = E_type(0);

	// trace seed
	// trace each seed base pair (excluding right most)
	for( size_t bpIn=1+bpInbetween; bpIn-- > 0; ) {

		// base case: only left and right base pair present
		if (bpIn==0) {
			// add left base pair if not left seed boundary
			if (i1 != i1_) {
				interaction.basePairs.push_back( energy.getBasePair(i1,i2) );
			}

		} else {

			// if no LP : check for direct right-stacking of i
			if (noLpShift > 0) {
				// check if feasible extension
				assert(isFeasibleSeedBasePair(i1+1,i2+1));
				// get stacking energy
				iStackE = energy.getE_interLeft(i1,i1+1,i2,i2+1);
				// noLP : check stacking of i
				if ( E_equal( curE, iStackE + getSeedE( i1+1, i2+1, bpIn-1, u1max, u2max )) ) {
					// store left base pair if not left seed boundary
					if (i1 != i1_) {
						interaction.basePairs.push_back( energy.getBasePair(i1,i2) );
					}
					i1++;
					i2++;
					curE = getSeedE( i1, i2, bpIn-1, u1max, u2max );
					continue;
				}
				// sanity check for noLP mode
				assert( bpIn >= 1+noLpShift );
			}

			// split seed recursively into all possible leading interior loops
			// i1 .. i1+u1p+1 .. j1
			// i2 .. i2+u2p+1 .. j2
			bool traceNotFound = true;
			for (u1=0; traceNotFound && u1<=u1max ; u1++) {
			for (u2=0; traceNotFound && u2<=u2max && (u1+u2)<=uMax ; u2++) {
				// check if overall number of unpaired is not exceeded
				// or skip stacked extension since covered above
				if (u1+u2 < noLpShift) {
					continue;
				}

				k1 = i1+u1+1+noLpShift;
				k2 = i2+u2+1+noLpShift;

				// check if valid trace
				if ( isFeasibleSeedBasePair(k1, k2) && E_isNotINF( getSeedE( k1, k2, bpIn-1-noLpShift, u1max-u1, u2max-u2 ) ) ) {

					// check if correct trace
					if ( E_equal( curE, iStackE
										+ energy.getE_interLeft(i1+noLpShift,k1,i2+noLpShift,k2)
										+ getSeedE( k1, k2, bpIn-1-noLpShift, u1max-u1, u2max-u2 )) )
					{
						// store next energy value to trace
						curE = getSeedE( k1, k2, bpIn-1-noLpShift, u1max-u1, u2max-u2 );
						// store left base pair if not left seed boundary
						if (i1 != i1_) {
							interaction.basePairs.push_back( energy.getBasePair(i1,i2) );
						}
						if (noLpShift > 0) {
							interaction.basePairs.push_back( energy.getBasePair(i1+noLpShift,i2+noLpShift) );
							// reflect additional base pair
							bpIn--;
						}
						// reset for next trace step
						i1 = k1;
						i2 = k2;
						// update boundaries for unpaired positions to reduce trace effort
						u1max -= u1;
						u2max -= u2;
						uMax -= (u1 + u2);
						// mark trace step done
						traceNotFound = false;
					}
				}

			} // u2
			} // u1
			assert( !traceNotFound ); // sanity check
		} // more than two base pairs

	} // bpIn

}

////////////////////////////////////////////////////////////////////////////

} // namespace
