#include "IntaRNA/HelixHandlerStackingOnly.h"

namespace IntaRNA {

/////////////////////////////////////////////////////////////////////////////

size_t
HelixHandlerStackingOnly::
fillHelix(const size_t i1min, const size_t i1max, const size_t i2min, const size_t i2max)
{

#if INTARNA_IN_DEBUG_MODE
	if ( i1min > i1max ) throw std::runtime_error("HelixHandlerStackingOnly::fillHelix: i1min("+toString(i1min)+") > i1max("+toString(i1max)+")");
	if ( i2min > i2max ) throw std::runtime_error("HelixHandlerStackingOnly::fillHelix: i2min("+toString(i2min)+") > i2max("+toString(i2max)+")");
	if ( i1max > energy.size1() ) throw std::runtime_error("HelixHandlerStackingOnly::fillHelix: i1max("+toString(i1max)+") > energy.size1("+toString(energy.size1())+")");
	if ( i2max > energy.size2() ) throw std::runtime_error("HelixHandlerStackingOnly::fillHelix: i2max("+toString(i2max)+") > energy.size2("+toString(energy.size2())+")");
	if ( helixConstraint.getMinBasePairs() > helixConstraint.getMaxBasePairs() )
		throw std::runtime_error("HelixHandlerStackingOnly::fillHelix: bpMin("+toString(helixConstraint.getMinBasePairs()) +") > bpMax("+toString(helixConstraint.getMaxBasePairs())+")");
#endif

	helix.resize( i1max-i1min+1, i2max-i2min+1 );
	helixE_rec.resize( HelixIndex({{(HelixRecMatrix::index)(helix.size1())
										   , (HelixRecMatrix::index)(helix.size2())
										   , (HelixRecMatrix::index)(getConstraint().getMaxBasePairs()+1)}}));

	// Initialize helixE_rec with E_INF in order to ensure that skipped values are not considered.
	std::fill( helixE_rec.origin(), helixE_rec.origin() + helixE_rec.num_elements(), E_INF);

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
	size_t i1, i2, curBP, j1, j2, k1, k2, bestBP;
	E_type curE, bestE;

	size_t helixCountNotInf = 0, helixCount = 0;

	// fill for all start indices
	// in decreasing index order
	for (i1=i1max+1; i1-- > i1min;) {
	for (i2=i2max+1; i2-- > i2min;) {

		// count possible helices
		helixCount++;

		// init according to no helix interaction
		helix(i1-offset1,i2-offset2) = HelixMatrix::value_type( E_INF, 0 );

		// skip non-complementary left helix boundaries
		if (!energy.areComplementary(i1,i2)) {
			continue; // go to next helixE index
		}

		// screen over all possible base pair combinations (starting at 2)
		for (curBP=2; curBP < helixConstraint.getMaxBasePairs()+1 && (i1+curBP-1-offset1<helix.size1())
					  											  && (i2+curBP-1-offset2<helix.size2()); curBP++) {

			// right helix boundaries
			j1 = i1+curBP-1;
			j2 = i2+curBP-1;

			// init current helix energy
			curE = E_INF;

			// Check whether right boundaries are complementary
			if (energy.areComplementary(j1,j2)) {
				// base case: only left and right base pair present
				if (curBP == 2) {
					curE = energy.getE_interLeft(i1,j1,i2,j2);
				} else {
					// splitting base pair
					k1 = i1 + 1;
					k2 = i2 + 1;

					// check if split pair is complementary
					// and recursed entry is < E_INF
					if (!(energy.areComplementary(k1,k2)) && E_isNotINF(getHelixE(k1-offset1, k2-offset2, curBP-1))) {
						continue; // not complementary -> skip
					}

					// update mfe for split at k1,k2
					curE = energy.getE_interLeft(i1,k1,i2,k2) + getHelixE(k1-offset1, k2-offset2, curBP-1);
				}
			} else {
				break;
			}

			// store helix energy
			setHelixE(i1-offset1, i2-offset2, curBP, curE);

		} // curBP

		// find best combination in helix for i1, i2, bp
		bestBP = 0;
		bestE = E_INF;

		// Calculate energy for all different number of base pairs
		// Ensure minimum number of base pairs
		for (curBP = minBP; curBP < helixConstraint.getMaxBasePairs() + 1
														&& (i1+curBP-1-offset1) < helix.size1()
														&& (i2+curBP-1-offset2) < helix.size2(); curBP++) {
			// right helix boundaries
			j1 = i1 + curBP-1;
			j2 = i2 + curBP-1;

			if (E_isINF(getHelixE(i1-offset1, i2-offset2, curBP))) {
				continue;
			}

			// ensure that ED-values are within the boundaries (default 999)
			if (energy.getED1(i1, j1) > helixConstraint.getMaxED()
				|| energy.getED2(i2, j2) > helixConstraint.getMaxED())
			{
				continue;
			}

			// calculate energy with ED values
			curE = energy.getE(i1,j1,i2,j2,getHelixE(i1-offset1, i2-offset2, curBP)) + energy.getE_init();

			// if noED option is set, ED-values are skipped, remove them.
			if (helixConstraint.useNoED())
				curE -= (energy.getED1(i1, j1) + energy.getED2(i2, j2));

			// check if better than what is know so far
			if (curE < bestE) {
				bestE = curE;
				bestBP = curBP;
			}
		} // curBP

		// reduce bestE to hybridization energy only (init+loops)
		if (E_isNotINF( bestE )) {
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

		// store best (mfe) helix for left boundary i1, i2
		helix(i1-offset1, i2-offset2) = HelixMatrix::value_type(bestE, E_isINF(bestE) ? 0: encodeHelixLength(bestBP, bestBP));
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
HelixHandlerStackingOnly::
traceBackHelix( Interaction & interaction
		, const size_t i1_
		, const size_t i2_
		, const size_t bp)
{

		// get boundaries
	size_t 	  i1 = i1_
	, i2 = i2_;

	// trace helices
	// trace each helix base pair (excluding right most)
	for (size_t curBP = 0; curBP < bp-1; curBP++) {
		if (i1 != i1_) {
			interaction.basePairs.push_back(energy.getBasePair(i1+offset1, i2+offset2));
		}
		i1++;
		i2++;
	}

}

} // namespace