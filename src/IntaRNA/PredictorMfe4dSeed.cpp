
#include "IntaRNA/PredictorMfe4dSeed.h"

#include <stdexcept>
#include <algorithm>

namespace IntaRNA {

////////////////////////////////////////////////////////////////////////////

PredictorMfe4dSeed::
PredictorMfe4dSeed( const InteractionEnergy & energy
					, OutputHandler & output
					, PredictionTracker * predTracker
					, SeedHandler * seedHandlerInstance )
 : PredictorMfe4d(energy,output,predTracker)
	, seedHandler(seedHandlerInstance)
	, hybridE_seed(0,0)
{
}


////////////////////////////////////////////////////////////////////////////

PredictorMfe4dSeed::
~PredictorMfe4dSeed()
{
	// clean up
	this->clear();
}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfe4dSeed::
predict( const IndexRange & r1
		, const IndexRange & r2
		, const OutputConstraint & outConstraint
		)
{
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) <<"predicting mfe interactions with seed in O(n^4) space and time..."; }
	// measure timing
	TIMED_FUNC_IF(timerObj,VLOG_IS_ON(9));

#if INTARNA_IN_DEBUG_MODE
	// check indices
	if (!(r1.isAscending() && r2.isAscending()) )
		throw std::runtime_error("PredictorMfe4dSeed::predict("+toString(r1)+","+toString(r2)+") is not sane");
#endif

	// clear data
	clear();

	// setup index offset
	energy.setOffset1(r1.from);
	energy.setOffset2(r2.from);
	seedHandler.setOffset1(r1.from);
	seedHandler.setOffset2(r2.from);

	const size_t hybridEsize1 = std::min( energy.size1()
			, (r1.to==RnaSequence::lastPos?energy.size1()-1:r1.to)-r1.from+1 );
	const size_t hybridEsize2 = std::min( energy.size2()
			, (r2.to==RnaSequence::lastPos?energy.size2()-1:r2.to)-r2.from+1 );


	// compute seed interactions for whole range
	// and check if any seed possible
	if (seedHandler.fillSeed( 0, hybridEsize1-1, 0, hybridEsize2-1 ) == 0) {
		// trigger empty interaction reporting
		initOptima(outConstraint);
		reportOptima(outConstraint);
		// stop computation
		return;
	}


	// resize matrix
	hybridE.resize( hybridEsize1, hybridEsize2 );
	hybridE_seed.resize( hybridEsize1, hybridEsize2 );

	size_t debug_count_cells_null=0
			, debug_count_cells_nonNull = 0
			, debug_count_cells_inf = 0
			, debug_cellNumber=0
			, w1, w2;


	size_t maxWidthFori1i2 = 0;

	bool i1blocked, i1or2blocked, skipw1w2;
	// initialize 3rd and 4th dimension of the matrix
	for (size_t i1=0; i1<hybridEsize1; i1++) {
		// check if i1 is blocked for interaction
		i1blocked = !energy.isAccessible1(i1);
		for (size_t i2=0; i2<hybridEsize2; i2++) {
			// check whether i1 or i2 is blocked for interaction
			i1or2blocked = i1blocked || !energy.isAccessible2(i2);

			if (hybridEsize1-i1 < hybridEsize2-i2) {
				maxWidthFori1i2 = getMaxInteractionWidth( hybridEsize1-i1, energy.getMaxInternalLoopSize1() );
			} else {
				maxWidthFori1i2 = getMaxInteractionWidth( hybridEsize2-i2, energy.getMaxInternalLoopSize2() );
			}

			debug_cellNumber =
					/*w1 = */ std::min(energy.getAccessibility1().getMaxLength(), std::min( hybridEsize1-i1, maxWidthFori1i2) )
				*	/*w2 = */ std::min(energy.getAccessibility2().getMaxLength(), std::min( hybridEsize2-i2, maxWidthFori1i2) );

			// check if i1 and i2 are not blocked and can form a base pair
			if ( ! i1or2blocked
				&& energy.areComplementary( i1, i2 ))
			{
				// create new 2d matrix for different interaction site widths
				hybridE(i1,i2) = new E2dMatrix(
					/*w1 = */ std::min(energy.getAccessibility1().getMaxLength(), std::min( hybridEsize1-i1, maxWidthFori1i2) ),
					/*w2 = */ std::min(energy.getAccessibility2().getMaxLength(), std::min( hybridEsize2-i2, maxWidthFori1i2) ));
				hybridE_seed(i1,i2) = new E2dMatrix( hybridE(i1,i2)->size1(), hybridE(i1,i2)->size2() );

				debug_count_cells_nonNull += debug_cellNumber;

				// screen for cells that can be skipped from computation (decreasing window sizes)
				for (size_t w1x = (*hybridE(i1,i2)).size1(); w1x>0; w1x--) {
					w1 = w1x-1;

				for (size_t w2x = (*hybridE(i1,i2)).size2(); w2x>0; w2x--) {
					w2 = w2x-1;

					// check if window size too large
					skipw1w2 = false;

					// check if ED penalty exceeds maximal energy gain
					if (w1 > 0 && w2 > 0){
						// check if all larger windows needing this site are already set to INF
						bool largerWindowsINF = w1x==(*hybridE(i1,i2)).size1() && w2x==(*hybridE(i1,i2)).size2();
						// check all larger windows w1 + w2p (that might need this window for computation)
						for (size_t w2p=(*hybridE(i1,i2)).size2()-1; largerWindowsINF && w2p>w2; w2p++) {
							// check if larger window is E_INF
							largerWindowsINF = (std::numeric_limits<E_type>::max() < (*hybridE(i1,i2))(w1+1,w2p));
						}
						// check all larger windows w2 + w1p (that might need this window for computation)
						for (size_t w1p=(*hybridE(i1,i2)).size1()-1; largerWindowsINF && w1p>w1; w1p++) {
							// check if larger window is E_INF
							largerWindowsINF = (std::numeric_limits<E_type>::max() < (*hybridE(i1,i2))(w1p,w2+1));
						}

						// if it holds for all w'>=w: ED1(i1+w1')+ED2(i2+w2')-outConstraint.maxE > -1*(min(w1',w2')*EmaxStacking + Einit + 2*Edangle + 2*Eend)
						// ie. the ED values exceed the max possible energy gain of an interaction
						skipw1w2 = skipw1w2
								|| ( largerWindowsINF &&
										( -1.0*(std::min(w1,w2)*minStackingEnergy + minInitEnergy + 2.0*minDangleEnergy + 2.0*minEndEnergy) <
											(energy.getED1(i1,i1+w1) + energy.getED2(i2,i2+w2) - outConstraint.maxE))
									)
									;
					}

					if (skipw1w2) {
						// init with infinity to mark that this cell is not to be computed later on
						(*hybridE(i1,i2))(w1,w2) = E_INF;
						(*hybridE_seed(i1,i2))(w1,w2) = E_INF;
						debug_count_cells_inf++;
					}

				}
				}

			} else {
				// reduce memory consumption and avoid computation for this start index combination
				hybridE(i1,i2) = NULL;
				hybridE_seed(i1,i2) = NULL;
				debug_count_cells_null += debug_cellNumber;
			}
		}
	}

#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_logOutput)
#endif
	{ LOG(DEBUG) <<"init 2x 4d matrix : "<<(debug_count_cells_nonNull-debug_count_cells_inf)<<" (-"<<debug_count_cells_inf <<") to be filled ("
				<<((double)(debug_count_cells_nonNull-debug_count_cells_inf)/(double)(debug_count_cells_nonNull+debug_count_cells_null))
				<<"%) and "<<debug_count_cells_null <<" not allocated ("
				<<((double)(debug_count_cells_null)/(double)(debug_count_cells_nonNull+debug_count_cells_null))
				<<"%)"; }

	// init mfe without seed condition
	OutputConstraint tmpOutConstraint(1, outConstraint.reportOverlap, outConstraint.maxE, outConstraint.deltaE);
	initOptima( tmpOutConstraint );

	// fill matrix
	// compute hybridization energies WITHOUT seed condition
	// sets also -energy -hybridE
	// -> no tracker update since updateOptima overwritten
	fillHybridE( );

	// check if any interaction possible
	// if not no seed-containing interaction is possible neither
	if (!(this->mfeInteractions.begin()->energy < tmpOutConstraint.maxE)) {
		// stop computation since no favorable interaction found
		reportOptima(tmpOutConstraint);
		return;
	}

	// initialize mfe interaction with seed for updates
	initOptima( outConstraint );
	// fill final matrix
	fillHybridE_seed( );

	// report mfe interaction
	reportOptima( outConstraint );

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe4dSeed::
clear()
{
	// delete 3rd and 4th dimension of the matrix
	for (E4dMatrix::iterator1 iRows = hybridE_seed.begin1(); iRows != hybridE_seed.end1(); iRows++) {
		for (E4dMatrix::iterator2 ijEntry = iRows.begin(); ijEntry != iRows.end(); ijEntry++) {
			// delete 2d matrix for current ij
			 INTARNA_CLEANUP(*ijEntry);
		}
	}
	// clear matrix, free data
	hybridE_seed.clear();

	// clean up super class data structures
	PredictorMfe4d::clear();
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe4dSeed::
fillHybridE_seed( )
{

	// global vars to avoid reallocation
	size_t i1,i2,j1,j2,w1,w2,k1,k2;

	//////////  FIRST ROUND : COMPUTE HYBRIDIZATION ENERGIES ONLY  ////////////

	// current minimal value
	E_type curMinE = E_INF;
	// iterate increasingly over all window sizes w1 (seq1) and w2 (seq2)
	// minimal size == number of seed base pairs
	for (w1=0; w1<energy.getAccessibility1().getMaxLength(); w1++) {
	for (w2=0; w2<energy.getAccessibility2().getMaxLength(); w2++) {
		// iterate over all window starts i1 (seq1) and i2 (seq2)
		// TODO PARALLELIZE THIS DOUBLE LOOP ?!
		for (i1=0; i1+w1<hybridE_seed.size1(); i1++) {
		for (i2=0; i2+w2<hybridE_seed.size2(); i2++) {
			// check if left boundary is complementary
			if (hybridE_seed(i1,i2) == NULL) {
				// interaction not possible: nothing to do, since no storage reserved
				continue;
			}

			// check if widths are possible (ie available within data structure)
			if ( w1 >= hybridE_seed(i1,i2)->size1() || w2 >= hybridE_seed(i1,i2)->size2() ) {
				// interaction not possible: nothing to do, since no storage reserved
				continue;
			}
			assert(hybridE(i1,i2)->size1() > w1);
			assert(hybridE(i1,i2)->size2() > w2);
			// check if no interaction without seed possible -> none with neither
			if ( E_isINF((*hybridE(i1,i2))(w1,w2))
					// check if seed base pairs fit not into interaction window
					|| w1+1 < seedHandler.getConstraint().getBasePairs()
					|| w2+1 < seedHandler.getConstraint().getBasePairs() )
			{
				// ignore this entry
				(*hybridE_seed(i1,i2))(w1,w2) = E_INF;
				continue;
			}

			// get window ends j1 (seq1) and j2 (seq2)
			j1=i1+w1;
			j2=i2+w2;

			// compute entry  ///////////////////////////
			curMinE = E_INF;

			// base case = incorporate mfe seed starting at (i1,i2)
			//             + interaction on right side up to (j1,j2)
			// -> check if widths are wide enough to contain a seed
			if ( E_isNotINF( seedHandler.getSeedE(i1,i2) ) ) {
				// decode right mfe boundary
				k1 = i1+seedHandler.getSeedLength1(i1,i2)-1;
				k2 = i2+seedHandler.getSeedLength2(i1,i2)-1;
				// compute overall energy of seed+upToPQ
				if ( k1 <= j1 && k2 <= j2
						&& hybridE(k1,k2) != NULL
						&& j1-k1 < hybridE(k1,k2)->size1()
						&& j2-k2 < hybridE(k1,k2)->size2()
						&& E_isNotINF((*hybridE(k1,k2))(j1-k1,j2-k2)))
				{
					curMinE = seedHandler.getSeedE(i1,i2) + (*hybridE(k1,k2))(j1-k1,j2-k2);
				}
			}

			// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
			// where k1..j1 contains a seed
			for (k1=std::min(j1-seedHandler.getConstraint().getBasePairs()+1,i1+energy.getMaxInternalLoopSize1()+1); k1>i1; k1--) {
			for (k2=std::min(j2-seedHandler.getConstraint().getBasePairs()+1,i2+energy.getMaxInternalLoopSize2()+1); k2>i2; k2--) {
				// check if (k1,k2) are valid left boundaries including a seed
				if ( hybridE_seed(k1,k2) != NULL
						&& j1-k1 < hybridE_seed(k1,k2)->size1()
						&& j2-k2 < hybridE_seed(k1,k2)->size2()
						&& E_isNotINF( (*hybridE_seed(k1,k2))(j1-k1,j2-k2) ) )
				{
					curMinE = std::min( curMinE,
							(energy.getE_interLeft(i1,k1,i2,k2)
									+ (*hybridE_seed(k1,k2))(j1-k1,j2-k2) )
						);
				}
			}
			}

			// store computed value
			(*hybridE_seed(i1,i2))(w1,w2) = curMinE;

			// update mfe if needed (call super class)
			if (E_isNotINF(curMinE)) {
				PredictorMfe4d::updateOptima( i1,j1,i2,j2, curMinE, true );
			}

		}
		}
	}
	}

}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfe4dSeed::
traceBack( Interaction & interaction, const OutputConstraint & outConstraint  )
{
	// check if something to trace
	if (interaction.basePairs.size() < 2) {
		return;
	}

#if INTARNA_IN_DEBUG_MODE
	// sanity checks
	if ( ! interaction.isValid() ) {
		throw std::runtime_error("PredictorMfe4dSeed::traceBack() : given interaction not valid");
	}
	if ( interaction.basePairs.size() != 2 ) {
		throw std::runtime_error("PredictorMfe4dSeed::traceBack() : given interaction does not contain boundaries only");
	}
#endif


	// check for single interaction (start==end)
	if (interaction.basePairs.begin()->first == interaction.basePairs.rbegin()->first) {
		// delete second boundary (identical to first)
		interaction.basePairs.resize(1);
		// update done
		return;
	}

	// ensure sorting
	interaction.sort();
	// get indices in hybridE_seed for boundary base pairs
	size_t	i1 = energy.getIndex1(interaction.basePairs.at(0)),
			j1 = energy.getIndex1(interaction.basePairs.at(1)),
			i2 = energy.getIndex2(interaction.basePairs.at(0)),
			j2 = energy.getIndex2(interaction.basePairs.at(1)),
			k1, k2;

	// the currently traced value for i1-j1, i2-j2
	E_type curE = (*hybridE_seed(i1,i2))(j1-i1,j2-i2);

	// trace back
	// do until only right boundary is left over
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
				if ( hybridE_seed(k1,k2) != NULL
						&& E_equal( curE, seedHandler.getSeedE(i1,i2) + (*hybridE(k1,k2))(j1-k1,j2-k2) ) )
				{
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
					curE = (*hybridE(k1,k2))(j1-k1,j2-k2);
					seedNotTraced = false;
					continue;
				}
			}
			// check all interval splits if no trace already found
			if ( seedNotTraced )
			{
				// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
				// where k1..j1 contains a seed
				bool traceNotFound = true;
				for (k1=std::min(j1-seedHandler.getConstraint().getBasePairs()+1,i1+energy.getMaxInternalLoopSize1()+1); traceNotFound && k1>i1; k1--) {
				for (k2=std::min(j2-seedHandler.getConstraint().getBasePairs()+1,i2+energy.getMaxInternalLoopSize2()+1); traceNotFound && k2>i2; k2--) {
					// check if (k1,k2) are valid left boundaries including a seed
					if ( 	hybridE_seed(k1,k2) != NULL
							&& j1-k1 < hybridE_seed(k1,k2)->size1()
							&& j2-k2 < hybridE_seed(k1,k2)->size2()
							&& E_isNotINF( (*hybridE_seed(k1,k2))(j1-k1,j2-k2) ) )
					{
						// check if correct split
						if (E_equal ( curE,
								(energy.getE_interLeft(i1,k1,i2,k2)
										+ (*hybridE_seed(k1,k2))(j1-k1,j2-k2) )
								) )
						{
							// update trace back boundary
							i1=k1;
							i2=k2;
							curE= (*hybridE_seed(k1,k2))(j1-k1,j2-k2);
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
			PredictorMfe4d::traceBack( rightSide, outConstraint );
			// copy base pairs (excluding last)
			for (size_t i=0; i+1<rightSide.basePairs.size(); i++) {
				interaction.basePairs.push_back( rightSide.basePairs.at(i) );
			}
			i1 = j1;
			i2 = j2;
			// stop traceback
			break;
		}
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
PredictorMfe4dSeed::
getNextBest( Interaction & curBest )
{

	// store current best energy.second
	const E_type curBestE = curBest.energy;

	// overwrite energy for update
	curBest.energy = E_INF;
	curBest.basePairs.resize(2);

	// TODO replace index iteration with something based on ranges from reportedInteractions

	// identify cell with next best non-overlapping interaction site
	// iterate (decreasingly) over all left interaction starts
	E2dMatrix * curTable = NULL;
	IndexRange r1,r2;
	size_t d1 = 0, d2 = 0;  // temp vars to deals with possible interaction lengths
	E_type curE = E_INF;
	for (r1.from=hybridE_seed.size1(); r1.from-- > 0;) {

		// ensure interaction site start is not covered
		if (reportedInteractions.first.covers(r1.from)) {
			continue;
		}

		for (r2.from=hybridE_seed.size2(); r2.from-- > 0;) {

			// ensure interaction site start is not covered
			if (reportedInteractions.second.covers(r2.from)) {
				continue;
			}
			// check if left boundary is complementary
			if (hybridE_seed(r1.from,r2.from) == NULL) {
				// interaction not possible: nothing to do, since no storage reserved
				continue;
			}

			// access energy table for left-most interaction base pair
			curTable = hybridE_seed(r1.from,r2.from);

			// iterate over all available interaction site lengths in seq1
			for (d1 = 0; d1<curTable->size1(); d1++) {

				// set according right interaction boundary in seq1
				r1.to = r1.from + d1;
				// check of overlapping
				if (reportedInteractions.first.overlaps(r1)) {
					// stop since all larger sites will overlap as well
					break;;
				}

				// iterate over all available interaction site lengths in seq2
				for (d2=0; d2<curTable->size2(); d2++) {

					// set according right interaction boundary in seq2
					r2.to = r2.from + d2;
					// check of overlapping
					if (reportedInteractions.second.overlaps(r2)) {
						// stop since all larger sites will overlap as well
						break;;
					}

					// get overall energy of entry
					curE = energy.getE( r1.from, r1.to, r2.from, r2.to, (*curTable)(d1,d2));

					// skip sites with energy too low
					// or higher than current best found so far
					if (  curE< curBestE || curE >= curBest.energy ) {
						continue;
					}

					//// FOUND THE NEXT BETTER SOLUTION
					// overwrite current best found so far
					curBest.energy = curE;
					curBest.basePairs[0] = energy.getBasePair( r1.from, r2.from );
					curBest.basePairs[1] = energy.getBasePair( r1.to, r2.to );

				} // j2
			} // j1


		} // i2
	} // i1

}


////////////////////////////////////////////////////////////////////////////


} // namespace
