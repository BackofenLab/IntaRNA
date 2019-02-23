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


	helix.resize( i1max-i1min+1, i2max-i2min+1 );
	helixE_rec.resize( HelixIndex({{
				   (HelixRecMatrix::index)(helix.size1())
				   , (HelixRecMatrix::index)(helix.size2())
				   , (HelixRecMatrix::index)(getConstraint().getMaxBasePairs()+1)}}) );
	std::fill( helixE_rec.origin(), helixE_rec.origin() + helixE_rec.num_elements(), std::make_pair(E_INF,0));

	size_t minBP;
	E_type maxE;
	// TODO: Fix until better solution is found
	// if no seedHandler is given use given values, else use default values to ensure they dont compromise the results of the helixSeed method
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

	size_t  helixCountNotInf = 0, helixCount = 0;

	// fill for all start indeices
	// in decreasing index order
	for (i1=i1max+1; i1-- > i1min;) {
	for (i2=i2max+1; i2-- > i2min;) {
		// count possible helices
		helixCount++;
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

			// get overall interaction energy
			curE = energy.getE(i1, j1, i2, j2, getHelixE(i1 - offset1, i2 - offset2, curBP)) + energy.getE_init();

			// skip if ED boundary exceeded and ED value computation is disabled
			if (helixConstraint.useNoED())
			{
				curE -= (energy.getED1(i1,j1) + energy.getED2(i2,j2));
			}

			// check if better than what is known so far
			if (curE < bestE && !E_equal(curE, bestE)) {
				bestE = curE;
				bestBP = curBP;
			}
		} // curBP

		// reduce bestE to hybridization energy only (init+loops)
		if (E_isNotINF(bestE)) {
			// overwrite all helices with too high energy -> infeasible start interactions
			if (bestE > maxE) {
				bestE = E_INF;
			} else {
				// get helix hybridization loop energies only
				bestE = getHelixE(i1 - offset1, i2 - offset2, bestBP);
				// count true helix
				helixCountNotInf++;
			}
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
	{ VLOG(2) << "valid helices = " << helixCountNotInf << " (" << (helixCountNotInf/helixCount) << "% of start index combinations)"; }

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


} // namespace


