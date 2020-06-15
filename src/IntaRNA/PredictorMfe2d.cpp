
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
		, const IndexRange & r2 )
{
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) <<"predicting mfe interactions in O(n^2) space..."; }
	// measure timing
	TIMED_FUNC_IF(timerObj,VLOG_IS_ON(9));

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

			// fill matrix and store best interaction
			fillHybridE( j1, j2, 0, 0, true );

		}
	}

	// report mfe interaction
	reportOptima();
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2d::
fillHybridE( const size_t j1, const size_t j2
			, const size_t i1init, const size_t i2init
			, const bool callUpdateOptima )
{
	// temporary access
	const OutputConstraint & outConstraint = output.getOutputConstraint();
#if INTARNA_IN_DEBUG_MODE
	if (i1init > j1)
		throw std::runtime_error("PredictorMfe2d::fillHybridE() : i1init > j1 : "+toString(i1init)+" > "+toString(j1));
	if (i2init > j2)
		throw std::runtime_error("PredictorMfe2d::fillHybridE() : i2init > j2 : "+toString(i2init)+" > "+toString(j2));
#endif

	// get minimal start indices heeding max interaction length
	const size_t i1start = std::max(i1init,j1-std::min(j1,energy.getAccessibility1().getMaxLength()+1));
	const size_t i2start = std::max(i2init,j2-std::min(j2,energy.getAccessibility2().getMaxLength()+1));

	// global vars to avoid reallocation
	size_t i1,i2,w1,w2,k1,k2;

	// determine whether or not lonely base pairs are allowed or if we have to
	// ensure a stacking to the right of the left boundary (i1,i2)
	const size_t noLpShift = outConstraint.noLP ? 1 : 0;
	E_type iStackE = E_type(0);

	//////////  COMPUTE HYBRIDIZATION ENERGIES  ////////////

	// iterate over all window starts i1 (seq1) and i2 (seq2)
	for (i1=j1+1; i1-- > i1start; ) {
		// w1 = interaction width in seq1
		w1 = j1-i1+1;
		// screen for left boundaries in seq2
		for (i2=j2+1; i2-- > i2start; ) {

			// init: mark as invalid boundary
			hybridE_pq(i1,i2) = E_INF;

			// check if this cell is to be computed (!=E_INF)
			if( energy.isAccessible1(i1)
				&& energy.isAccessible2(i2)
				&& energy.areComplementary(i1,i2)
			)
			{
				// w2 = interaction width in seq2
				w2 = j2-i2+1;

				// reference access to cell value
				E_type &curMinE = hybridE_pq(i1,i2);

				// either interaction initiation
				if ( i1==j1 && i2==j2)  {
					if (noLpShift == 0) {
						// single base pair
						curMinE = energy.getE_init();
					}
				}
				else
				// or at least two base pairs possible
				if ( w1 > 1 && w2 > 1) {
					// init curMinE
					// if lonely bps are allowed
					if (noLpShift == 0) {
						// test full-width internal loop energy (nothing between i and j)
						// will be E_INF if loop is too large
						curMinE = energy.getE_interLeft(i1,j1,i2,j2)
								+ hybridE_pq(j1,j2);
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
					// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
					// ensure stacking is possible if no LP allowed
					if (w1 > 2 && w2 > 2 && E_isNotINF(iStackE)) {
						for (k1=std::min(j1-1,i1+energy.getMaxInternalLoopSize1()+1+noLpShift); k1>i1+noLpShift; k1--) {
						for (k2=std::min(j2-1,i2+energy.getMaxInternalLoopSize2()+1+noLpShift); k2>i2+noLpShift; k2--) {
							// check if (k1,k2) are valid left boundary
							if ( E_isNotINF( hybridE_pq(k1,k2) ) ) {
								// update minimal value
								curMinE = std::min( curMinE,
										( iStackE + energy.getE_interLeft(i1+noLpShift,k1,i2+noLpShift,k2)
												+ hybridE_pq(k1,k2) )
										);
							}
						}
						}
					}
				}

				// update mfe if needed
				if (callUpdateOptima) {
					updateOptima( i1,j1,i2,j2, curMinE, true, true );
				}

			} // complementary base pair
		}
	}

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2d::
traceBack( Interaction & interaction  )
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
		throw std::runtime_error("PredictorMfe2d::traceBack() : given interaction does not contain boundaries only");
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
		throw std::runtime_error("PredictorMfe2d::traceBack() : given interaction is not valid");
	}
#endif

	// get indices in hybridE for boundary base pairs
	size_t	i1 = energy.getIndex1(interaction.basePairs.at(0)),
			j1 = energy.getIndex1(interaction.basePairs.at(1)),
			i2 = energy.getIndex2(interaction.basePairs.at(0)),
			j2 = energy.getIndex2(interaction.basePairs.at(1));

	// refill submatrix of mfe interaction
	fillHybridE( j1, j2, i1, i2, false );

	// the currently traced value for i1-j1, i2-j2
	E_type curE = hybridE_pq(i1,i2);

	// determine whether or not lonely base pairs are allowed or if we have to
	// ensure a stacking to the right of the left boundary (i1,i2)
	const size_t noLpShift = outConstraint.noLP ? 1 : 0;
	E_type iStackE = E_type(0);

	// trace back
	while( i1+noLpShift < j1 ) {

		// check if just internal loop
		if ( noLpShift == 0 )
		{
			if ( E_equal( curE, (energy.getE_interLeft(i1,j1,i2,j2) + hybridE_pq(j1,j2)) ) )
				break;
		} else {
			// check if only 3 stack helix possible
			if (j1-i1+1 == 3) {
				interaction.basePairs.push_back( energy.getBasePair(i1+1,i2+1) );
				break;
			}
			// get stacking term
			iStackE = energy.getE_interLeft(i1,i1+noLpShift, i2,i2+noLpShift);

			// check just stacking
			if ( E_equal( curE, iStackE + hybridE_pq(i1+noLpShift,i2+noLpShift) )) {
				// store bp
				interaction.basePairs.push_back( energy.getBasePair(i1+noLpShift,i2+noLpShift) );
				// update trace
				i1 += noLpShift;
				i2 += noLpShift;
				curE = hybridE_pq(i1,i2);
				// next trace step
				continue;
			}

		}
		// check all interval splits
		if ( (j1-i1) > 1 && (j2-i2) > 1) {
			// temp variables
			size_t k1,k2;
			bool traceNotFound = true;
			// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
			for (k1=std::min(j1-1,i1+energy.getMaxInternalLoopSize1()+1+noLpShift); traceNotFound && k1>i1+noLpShift; k1--) {
			for (k2=std::min(j2-1,i2+energy.getMaxInternalLoopSize2()+1+noLpShift); traceNotFound && k2>i2+noLpShift; k2--) {
				// check if (k1,k2) are valid left boundary
				if ( E_isNotINF( hybridE_pq(k1,k2) ) ) {
					if ( E_equal( curE,
							(iStackE + energy.getE_interLeft(i1+noLpShift,k1,i2+noLpShift,k2) + hybridE_pq(k1,k2)) ) )
					{
						// stop searching
						traceNotFound = false;
						// if no lonely base pairs : store stacking base pair
						if (noLpShift > 0) {
							interaction.basePairs.push_back( energy.getBasePair(i1+noLpShift,i2+noLpShift) );
						}
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
			// avoid infinity loops
			if (traceNotFound) {
				throw std::runtime_error("PredictorMfe2d::traceBack() could not traceback base pair ("+toString(energy.getBasePair(i1,i2).first)+","+toString(energy.getBasePair(i1,i2).second)+") of interaction "+toString(interaction));
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


} // namespace
