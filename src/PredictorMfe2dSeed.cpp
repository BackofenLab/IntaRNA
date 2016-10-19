
#include "PredictorMfe2dSeed.h"

#include "easylogging++.h"

//////////////////////////////////////////////////////////////////////////

PredictorMfe2dSeed::
PredictorMfe2dSeed(
		const InteractionEnergy & energy
		, OutputHandler & output
		, const SeedConstraint& seedConstraint )
 :
	PredictorMfe2d(energy,output)
	, seedConstraint(seedConstraint)
	, hybridE_pq_seed()
	, seedE_rec( SeedIndex({{ // setup ring-list data for seed computation
			  (SeedRecMatrix::index)(seedConstraint.getBasePairs()+seedConstraint.getMaxUnpaired1()+1)
			, (SeedRecMatrix::index)(seedConstraint.getBasePairs()+seedConstraint.getMaxUnpaired2()+1)
			, (SeedRecMatrix::index)(seedConstraint.getBasePairs()+1-2) // +1 for size and -2 to encode at least 2 bps or more
			, (SeedRecMatrix::index)(seedConstraint.getMaxUnpaired1()+1) // +1 for size
			, (SeedRecMatrix::index)(seedConstraint.getMaxUnpaired2()+1) // +1 for size
		}}))
	, seed()
{
	assert( seedConstraint.getBasePairs() > 1 );
}

//////////////////////////////////////////////////////////////////////////

PredictorMfe2dSeed::
~PredictorMfe2dSeed()
{
}

//////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dSeed::
predict( const IndexRange & r1, const IndexRange & r2 )
{

	VLOG(2) <<"predicting mfe interactions in O(n^2) space using seed constraint...";

#if IN_DEBUG_MODE
	// measure timing
	TIMED_FUNC(timerObj);
	// check indices
	if (!(r1.isAscending() && r2.isAscending()) )
		throw std::runtime_error("PredictorMfe2d::predict("+toString(r1)+","+toString(r2)+") is not sane");
#endif

	// resize matrix
	hybridE_pq.resize( std::min( energy.getAccessibility1().getSequence().size()
						, (r1.to==RnaSequence::lastPos?energy.getAccessibility1().getSequence().size()-1:r1.to)-r1.from+1 )
				, std::min( energy.getAccessibility2().getSequence().size()
						, (r2.to==RnaSequence::lastPos?energy.getAccessibility2().getSequence().size()-1:r2.to)-r2.from+1 ) );
	hybridE_pq_seed.resize( hybridE_pq.size1(), hybridE_pq.size2() );

	i1offset = r1.from;
	i2offset = r2.from;

	// compute seeds
	seed.resize( hybridE_pq.size1(), hybridE_pq.size2() );
	fillSeed( 0, seed.size1()-1, 0, seed.size2()-1 );

	// initialize mfe interaction for updates
	initMfe();

	// for all right ends j1
	for (size_t j1 = hybridE_pq.size1(); j1-- > 0; ) {
		// check if j1 is accessible
		if (energy.getAccessibility1().getAccConstraint().isBlocked(j1+i1offset))
			continue;
		// iterate over all right ends j2
		for (size_t j2 = hybridE_pq.size2(); j2-- > 0; ) {
			// check if j2 is accessible
			if (energy.getAccessibility2().getAccConstraint().isBlocked(j2+i2offset))
				continue;
			// check if base pair (j1,j2) possible
			if (!energy.areComplementary( j1+i1offset, j2+i2offset ))
				continue;

			// compute hybridE_pq
			fillHybridE( j1, j2 );

			// compute hybridE_pq_seed and update mfe via PredictorMfe2d::updateMfe()
			fillHybridE_seed( j1, j2 );
		}
	}

	// check if interaction is better than no interaction (E==0)
	if (mfeInteraction.energy < 0.0) {
		// fill mfe interaction with according base pairs
		traceBack( mfeInteraction );
	} else {
		// replace mfeInteraction with no interaction
		mfeInteraction.clear();
		mfeInteraction.energy = 0.0;
	}

	// report mfe interaction
	output.add( mfeInteraction );
}

//////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dSeed::
updateMfe( const size_t i1, const size_t j1
		, const size_t i2, const size_t j2
		, const E_type energy )
{
	// do nothing and ignore calls from fillHybridE()
}

//////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dSeed::
fillSeed( const size_t i1min, const size_t i1max, const size_t i2min, const size_t i2max)
{

#if IN_DEBUG_MODE
	if ( i1min > i1max ) throw std::runtime_error("PredictorMfe2dSeed::fillSeed: i1min("+toString(i1min)+") > i1max("+toString(i1max)+")");
	if ( i2min > i2max ) throw std::runtime_error("PredictorMfe2dSeed::fillSeed: i2min("+toString(i2min)+") > i2max("+toString(i2max)+")");
	if ( i1min > seed.size1() ) throw std::runtime_error("PredictorMfe2dSeed::fillSeed: i1min("+toString(i1min)+") > seed.size1("+toString(seed.size1())+")");
	if ( i2min > seed.size2() ) throw std::runtime_error("PredictorMfe2dSeed::fillSeed: i2min("+toString(i2min)+") > seed.size2("+toString(seed.size2())+")");
#endif

	// temporary variables
	size_t i1, i2, bp, u1, u2, j1, j2, u1p, u2p, k1,k2, u1best, u2best;
	E_type curE, bestE;

	// fill for all start indices
	// in decreasing index order
	for (i1=i1max+1; i1-- > i1min;) {
	for (i2=i2max+1; i2-- > i2min;) {

		// init according to no seed interaction
		seed(i1,i2) = SeedMatrix::value_type( E_INF, 0 );

		// skip non-complementary left seed boundaries
		if (!energy.areComplementary(i1,i2)) {
			continue; // go to next seedE index
		}

		// for feasible number of base pairs (bp+1) in increasing order
		// bp=0 encodes 2 base pairs
		for (bp=0; bp<seedE_rec.shape()[2] && (i1+bp+1)<seed.size1() && (i2+bp+1)<seed.size2(); bp++) {

			// for feasible unpaired in seq1 in increasing order
			for (u1=0; u1<seedE_rec.shape()[3] && (i1+bp+1+u1) < seed.size1(); u1++) {
			// for feasible unpaired in seq2 in increasing order
			for (u2=0; u2<seedE_rec.shape()[4] && (u1+u2)<=seedConstraint.getMaxUnpairedOverall() && (i2+bp+1+u2) < seed.size2(); u2++) {

				// get right seed boundaries
				j1 = i1+bp+1+u1;
				j2 = i2+bp+1+u2;

				// init current seed energy
				curE = E_INF;

				// check if right boundary is complementary
				if (energy.areComplementary(j1,j2)) {

					// base case: only left and right base pair present
					if (bp==0) {
						// energy for stacking/bulge/interior depending on u1/u2
						curE = energy.getE_interLeft(i1,j1,i2,j2);

					} else {
						// split seed recursively into all possible leading interior loops
						// i1 .. i1+u1p+1 .. j1
						// i2 .. i2+u2p+1 .. j2
						for (u1p=1+std::min(u1,energy.getMaxInternalLoopSize1()); u1p-- > 0;) {
						for (u2p=1+std::min(u2,energy.getMaxInternalLoopSize2()); u2p-- > 0;) {

							k1 = i1+u1p+1;
							k2 = i2+u2p+1;
							// check if split pair is complementary
							// and recursed entry is < E_INF
							if (! (energy.areComplementary(k1,k2) && E_isNotINF( getSeedE( k1, k2, bp-1, u1-u1p, u2-u2p ) ) ) ) {
								continue; // not complementary -> skip
							}

							// update mfe for split at k1,k2
							curE = std::min( curE,
									energy.getE_interLeft(i1+i1offset,k1+i1offset,i2+i2offset,k2+i2offset)
									+ getSeedE( k1, k2, bp-1, u1-u1p, u2-u2p )
									);
						} // u2p
						} // u1p
					} // more than two base pairs

				} // (j1,j2) complementary

				// store seed energy
				setSeedE( i1, i2, bp, u1, u2, curE );

			} // u2
			} // u1

			// check if full base pair number reached
			if (bp+1==seedE_rec.shape()[2]) {

				// find best unpaired combination in seed seed for i1,i2,bp
				u1best = 0;
				u2best = 0;
				bestE = E_INF;

				// for feasible unpaired in seq1 in increasing order
				for (u1=0; u1<seedE_rec.shape()[3] && (i1+bp+1+u1) < seed.size1(); u1++) {
				// for feasible unpaired in seq2 in increasing order
				for (u2=0; u2<seedE_rec.shape()[4] && (u1+u2)<=seedConstraint.getMaxUnpairedOverall() && (i2+bp+1+u2) < seed.size2(); u2++) {

					// get right seed boundaries
					j1 = i1+bp+1+u1;
					j2 = i2+bp+1+u2;

					// get overall interaction energy
					curE = energy.getE( i1, j1, i2, j2, getSeedE( i1, i2, bp, u1, u2 ) ) + energy.getE_init();

					// check if better than what is known so far
					if ( curE < bestE ) {
						bestE = curE;
						u1best = u1;
						u2best = u2;
					}
				} // u2
				} // u1

				// reduce bestE to hybridization energy only (init+loops)
				if (E_isNotINF( bestE )) {
					// overwrite all seeds with too high energy -> infeasible start interactions
					if (bestE > seedConstraint.getMaxE()) {
						bestE = E_INF;
					} else {
						// get seed's hybridization loop energies only
						bestE -= energy.getE( i1, i1+bp+1+u1best, i2, i2+bp+1+u2best, 0.0 );
						bestE -= energy.getE_init();
					}
				}

				// store best (mfe) seed for all u1/u2
				seed(i1,i2) = SeedMatrix::value_type( bestE
						, E_isINF(bestE)?0:encodeSeedLength(bp+2+u1best,bp+2+u2best) );

			} // store best seed

		} // bp
	} // i2
	} // i1

}

//////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dSeed::
fillHybridE_seed( const size_t j1, const size_t j2, const size_t i1init, const size_t i2init )
{
	assert(j1<=hybridErange.r1.to);
	assert(j2<=hybridErange.r2.to);
	assert(i1init <= j1);
	assert(i2init <= j2);
	assert(j1<hybridE_pq.size1());
	assert(j2<hybridE_pq.size2());
	assert(seed.size1() == hybridE_pq.size1());
	assert(seed.size2() == hybridE_pq.size2());


	// check if it is possible to have a seed ending on the right at (j1,j2)
	if (j1+1 < seedConstraint.getBasePairs() || j2+1 < seedConstraint.getBasePairs()) {
		// no seed possible, abort computation
		return;
	}

	// global vars to avoid reallocation
	size_t i1,i2,k1,k2;

	// get i1/i2 index boundaries for computation
	const IndexRange i1range( std::max(hybridErange.r1.from,i1init), j1+1-seedConstraint.getBasePairs() );
	const IndexRange i2range( std::max(hybridErange.r2.from,i2init), j2+1-seedConstraint.getBasePairs() );

	//////////  FIRST ROUND : COMPUTE HYBRIDIZATION ENERGIES ONLY  ////////////

	// current minimal value
	E_type curMinE = E_INF;
	// iterate over all window starts i1 (seq1) and i2 (seq2)
	// TODO PARALLELIZE THIS DOUBLE LOOP ?!
	for (i1=1+i1range.to; i1-- > i1range.from; ) {
		// screen for left boundaries i2 in seq2
		for (i2=1+i2range.to; i2-- > i2range.from; ) {

			// check if this cell is to be computed (!=E_INF)
			if( E_isNotINF( hybridE_pq(i1,i2) ) ) {

				// compute entry
				curMinE = E_INF;

				// base case = incorporate mfe seed starting at (i1,i2)
				//             + interaction on right side up to (p,q)
				const SeedMatrix::value_type & seedMfe = seed(i1,i2);
				if ( E_isNotINF( seedMfe.first ) && seedMfe.second > 0) {
					// decode right mfe boundary
					k1 = i1+decodeSeedLength1(seedMfe.second)-1;
					k2 = i2+decodeSeedLength2(seedMfe.second)-1;
					// compute overall energy of seed+upToPQ
					if ( k1 <= j1 && k2 <= j2 ) {
						curMinE = seedMfe.first + hybridE_pq(k1,k2);
					}
				}

				// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
				// where k1..j1 contains a seed
				for (k1=std::min(i1range.to,i1+energy.getMaxInternalLoopSize1()+1); k1>i1; k1--) {
				for (k2=std::min(i2range.to,i2+energy.getMaxInternalLoopSize2()+1); k2>i2; k2--) {
					// check if (k1,k2) are valid left boundaries including a seed
					if ( E_isNotINF( hybridE_pq_seed(k1,k2) ) ) {
						curMinE = std::min( curMinE,
								(energy.getE_interLeft(i1+i1offset,k1+i1offset,i2+i2offset,k2+i2offset)
										+ hybridE_pq_seed(k1,k2) )
								);
					}
				}
				}
				// store value
				hybridE_pq_seed(i1,i2) = curMinE;
				// update mfe if needed (call super class)
				PredictorMfe2d::updateMfe( i1,j1,i2,j2, hybridE_pq_seed(i1,i2) );
				continue;
			}
		}
	}

}

//////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dSeed::
traceBack( Interaction & interaction )
{
	// check if something to trace
	if (interaction.basePairs.size() < 2) {
		return;
	}

#if IN_DEBUG_MODE
	// sanity checks
	if ( ! interaction.isValid() ) {
		throw std::runtime_error("PredictorMfe2dSeed::traceBack() : given interaction not valid");
	}
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

	// ensure sorting
	interaction.sort();
	// get indices in hybridE for boundary base pairs
	size_t	i1 = interaction.basePairs.at(0).first - i1offset,
			j1 = interaction.basePairs.at(1).first - i1offset,
			i2 = energy.getAccessibility2().getReversedIndex(interaction.basePairs.at(0).second) - i2offset,
			j2 = energy.getAccessibility2().getReversedIndex(interaction.basePairs.at(1).second) - i2offset
			;

#if IN_DEBUG_MODE
	// check if it is possible to have a seed ending on the right at (j1,j2)
	if (j1+1 < seedConstraint.getBasePairs() || j2+1 < seedConstraint.getBasePairs()) {
		// no seed possible, abort computation
		throw std::runtime_error("PredictorMfe2dSeed::traceBack() : given boundaries "+toString(interaction)+" can not hold a seed of "+toString(seedConstraint.getBasePairs())+" base pairs");
	}
#endif

	// temp variables
	size_t k1,k2;


	// refill submatrices of mfe interaction
	fillHybridE( j1, j2, i1, i2 );
	fillHybridE_seed( j1, j2, i1, i2 );

	// the currently traced value for i1-j1, i2-j2
	E_type curE = hybridE_pq_seed(i1,i2);

	// trace back
	bool seedNotTraced = true;
	while( i1 != j1 ) {

		// check if we still have to find the seed
		if (seedNotTraced) {

			// check base case == seed only
			if ( E_isNotINF( seed(i1,i2).first ) ) {

				// right boundary of seed
				k1 = i1 + decodeSeedLength1( seed(i1,i2).second ) -1;
				k2 = i2 + decodeSeedLength2( seed(i1,i2).second ) -1;

				// check if correct trace
				if ( E_equal( curE, seed(i1,i2).first + hybridE_pq(k1,k2) ) ) {
					// trace back seed base pairs
					traceBackSeed( interaction, i1, i2, seedConstraint.getBasePairs()-2, k1-i1+1-seedConstraint.getBasePairs(), k2-i2+1-seedConstraint.getBasePairs() );
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
				for (k1=std::min(j1-seedConstraint.getBasePairs()+1,i1+energy.getMaxInternalLoopSize1()+1); traceNotFound && k1>i1; k1--) {
				for (k2=std::min(j2-seedConstraint.getBasePairs()+1,i2+energy.getMaxInternalLoopSize2()+1); traceNotFound && k2>i2; k2--) {
					// check if (k1,k2) are valid left boundaries including a seed
					if ( E_isNotINF( hybridE_pq_seed(k1,k2) ) ) {
						// check if correct split
						if (E_equal ( curE,
								(energy.getE_interLeft(i1+i1offset,k1+i1offset,i2+i2offset,k2+i2offset)
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
							interaction.addInteraction( k1+i1offset, energy.getAccessibility2().getReversedIndex(k2+i2offset) );
						}
					}
				}
				}
				assert(!traceNotFound);
			}
		}
		// seed was already traced, do "normal" interaction trace
		else {
			// create temporary data structure to be filed
			Interaction rightSide( *interaction.s1, *interaction.s2 );
			rightSide.addInteraction( i1+i1offset, energy.getAccessibility2().getReversedIndex(i2+i2offset) );
			rightSide.addInteraction( j1+i1offset, energy.getAccessibility2().getReversedIndex(j2+i2offset) );
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
		bps.rbegin()->first = j1+i1offset;
		bps.rbegin()->second = energy.getAccessibility2().getReversedIndex(j2+i2offset);
	}
}

//////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dSeed::
traceBackSeed( Interaction & interaction
		, const size_t i1_
		, const size_t i2_
		, const size_t bpInbetween
		, const size_t u1_
		, const size_t u2_
		)
{

	// get boundaries
	size_t 	  i1 = i1_
			, i2 = i2_
			, j1 = i1+bpInbetween+1+u1_
			, j2 = i2+bpInbetween+1+u2_
			, u1max = u1_
			, u2max = u2_
			, uMax = u1_+u2_
			, u1, u2, k1, k2
	   ;

	// recompute seedE_rec for this seed region (exclude right boundary, not needed)
	fillSeed( i1, j1-1, i2, j2-1 );

	// get energy of provided seed
	E_type curE = getSeedE(i1_,i2_,bpInbetween,u1_,u2_);

	// trace seed
	// trace each seed base pair (excluding right most)
	for( size_t bpIn=1+bpInbetween; bpIn-- > 0; ) {

		// base case: only left and right base pair present
		if (bpIn==0) {
			// add left base pair if not left seed boundary
			if (i1 != i1_) {
				interaction.addInteraction( i1+i1offset, energy.getAccessibility2().getReversedIndex(i2+i2offset) );
			}

		} else {
			// split seed recursively into all possible leading interior loops
			// i1 .. i1+u1p+1 .. j1
			// i2 .. i2+u2p+1 .. j2
			bool traceNotFound = true;
			for (u1=1+u1max; traceNotFound && u1-- > 0;) {
			for (u2=1+u2max; traceNotFound && u2-- > 0;) {
				// check if overall number of unpaired is not exceeded
				if (u1+u2 > uMax) {
					continue;
				}

				k1 = i1+u1+1;
				k2 = i2+u2+1;
				// check if base pair is possible
				if ( E_isINF( hybridE_pq(k1,k2) ) ){
					continue;
				}

				// check if valid trace
				if ( E_isNotINF( getSeedE( k1, k2, bpIn-1, u1max-u1, u2max-u2 ) ) ) {

					// check if correct trace
					if ( E_equal( curE, energy.getE_interLeft(i1,k1,i2,k2)
										+ getSeedE( k1, k2, bpIn-1, u1max-u1, u2max-u2 )) )
					{
						// store left base pair if not left seed boundary
						if (i1 != i1_) {
							interaction.addInteraction( i1+i1offset, energy.getAccessibility2().getReversedIndex(i2+i2offset) );
						}
						// store next energy value to trace
						curE = getSeedE( k1, k2, bpIn-1, u1max-u1, u2max-u2 );
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
		} // more than two base pairs

	} // bpIn

}

//////////////////////////////////////////////////////////////////////////

E_type
PredictorMfe2dSeed::
getSeedE( const size_t i1, const size_t i2, const size_t bpInbetween, const size_t u1, const size_t u2 )
{
	// access i1/i2 via modulo operation to get a ring-list like access behavior with mem-reusage
	return seedE_rec[i1 % seedE_rec.shape()[0]][i2 % seedE_rec.shape()[1]][bpInbetween][u1][u2];
//	return seedE_rec( SeedIndex({{ i1 % seedE_rec.shape()[0], i2 % seedE_rec.shape()[1], bpInbetween, u1, u2 }}) );
}

//////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dSeed::
setSeedE( const size_t i1, const size_t i2, const size_t bpInbetween, const size_t u1, const size_t u2, const E_type E )
{
	// access i1/i2 via modulo operation to get a ring-list like access behavior with mem-reusage
	seedE_rec[i1 % seedE_rec.shape()[0]][i2 % seedE_rec.shape()[1]][bpInbetween][u1][u2] = E;
//	seedE_rec( SeedIndex({{ i1 % seedE_rec.shape()[0], i2 % seedE_rec.shape()[1], bpInbetween, u1, u2 }}) ) = E;
}

//////////////////////////////////////////////////////////////////////////

size_t
PredictorMfe2dSeed::
encodeSeedLength( const size_t l1, const size_t l2 ) const
{
	return l1 + l2*(seedConstraint.getMaxLength1()+1);
}

//////////////////////////////////////////////////////////////////////////

size_t
PredictorMfe2dSeed::
decodeSeedLength1( const size_t code ) const
{
	return code % (seedConstraint.getMaxLength1()+1);
}

//////////////////////////////////////////////////////////////////////////

size_t
PredictorMfe2dSeed::
decodeSeedLength2( const size_t code ) const
{
	return code / (seedConstraint.getMaxLength1()+1);
}

//////////////////////////////////////////////////////////////////////////



