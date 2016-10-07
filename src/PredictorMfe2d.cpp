
#include "PredictorMfe2d.h"

#include <stdexcept>


////////////////////////////////////////////////////////////////////////////

PredictorMfe2d::
PredictorMfe2d( const InteractionEnergy & energy, OutputHandler & output )
 : PredictorMfe(energy,output)
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
#ifdef NDEBUG
	// no check
#else
	// check indices
	if (!(r1.isAscending() && r2.isAscending()) )
		throw std::runtime_error("PredictorMfe2d::predict("+toString(r1)+","+toString(r2)+") is not sane");
#endif


	// resize matrix
	hybridE_pq.resize( std::min( energy.getAccessibility1().getSequence().size()
						, (r1.to==RnaSequence::lastPos?energy.getAccessibility1().getSequence().size()-1:r1.to)-r1.from+1 )
				, std::min( energy.getAccessibility2().getSequence().size()
						, (r2.to==RnaSequence::lastPos?energy.getAccessibility2().getSequence().size()-1:r2.to)-r2.from+1 ) );

	i1offset = r1.from;
	i2offset = r2.from;

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

			// fill matrix and store best interaction
			fillHybridE( j1, j2 );

		}
	}

	// check if interaction is better than no interaction (E==0)
	if (mfeInteraction.energy < 0.0) {
		// fill mfe interaction with according base pairs
		traceBack( mfeInteraction );
	} else {
		// TODO : check if better to skip output handler report instead of overwrite
		// replace mfeInteraction with no interaction
		mfeInteraction.clear();
		mfeInteraction.energy = 0.0;
	}

	// report mfe interaction
	output.add( mfeInteraction );
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2d::
initHybridE( const size_t j1, const size_t j2 )
{

	// global vars to avoid reallocation
	size_t i1,i2,w1,w2,k1,k2;


	const E_type E_MAX = std::numeric_limits<E_type>::max();

	// init used part with E_INF - 1 and E_INF if ED test fails
	// iterating decreasing window size seq1 (only sane windows)
	const size_t minI1 = j1-std::min(j1,energy.getAccessibility1().getMaxLength()+2);
	const size_t minI2 = j2-std::min(j2,energy.getAccessibility2().getMaxLength()+2);
	for (i1=minI1; i1<=j1; i1++ ) {
		w1 = j1-i1+1;
		// iterating decreasing window size seq2 (only sane windows)
		for (i2=minI2; i2<=j2; i2++) {
			w2 = j2-i2+1;
			// check if ED penalty exceeds maximal energy gain
			{
				// check if all larger windows needing this site are already set to INF
				bool largerWindowsINF = i1==minI1 && i2==minI2;
				// check all larger windows w1 + i2p..j2 (that might need this window for computation)
				for (size_t i2p=minI2; largerWindowsINF && i2p>i2; i2p++) {
					// check if larger window is E_INF
					largerWindowsINF = (hybridE_pq(i1,i2p) < E_INF);
				}
				// check all larger windows w2 + w1p (that might need this window for computation)
				for (size_t i1p=minI1; largerWindowsINF && i1p>i1; i1p++) {
					// check if larger window is E_INF
					largerWindowsINF = (hybridE_pq(i1,i2) < E_INF);
				}

				// if it holds for all w'>=w: ED1(i1+w1')+ED2(i2+w2') > -1*(min(w1',w2')*EmaxStacking + Einit + 2*Edangle + 2*Eend)
				// ie. the ED values exceed the max possible energy gain of an interaction
				if( largerWindowsINF &&
								( -1.0*(std::min(w1,w2)*minStackingEnergy + minInitEnergy + 2.0*minDangleEnergy + 2.0*minEndEnergy) >
									(energy.getAccessibility1().getED(i1+i1offset,j1+i1offset)
											+ energy.getAccessibility2().getED(i2+i2offset,j2+i2offset)))
							)
				{
					// mark as NOT to be computed
					hybridE_pq(i1,i2) = E_INF;
					continue;
				}
			}

			// mark as to be computed (has to be < E_INF)
			hybridE_pq(i1,i2) = E_MAX;
		}
	}
}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2d::
fillHybridE( const size_t j1, const size_t j2 )
{

	// init for right interaction end (j1,j2)
	initHybridE( j1, j2 );

	// global vars to avoid reallocation
	size_t i1,i2,w1,w2,k1,k2;

	//////////  FIRST ROUND : COMPUTE HYBRIDIZATION ENERGIES ONLY  ////////////

	// current minimal value
	E_type curMinE = E_INF;
	// iterate over all window starts i1 (seq1) and i2 (seq2)
	// TODO PARALLELIZE THIS DOUBLE LOOP ?!
	for (i1=j1+1; i1-->0; ) {
		w1 = j1-i1+1;
		// check if maximal interaction length exceeded (seq1)
		if ( w1 > energy.getAccessibility1().getMaxLength() ) {
			// skip all remaining i1, since only increasing width
			break;
		}
		// get maximal w2 for this w1
		const size_t maxW2 = getMaxInteractionWidth( w1, energy.getMaxInternalLoopSize1());
		// screen for left boundaries in seq2
		for (i2=j2+1; i2-->0; ) {
			w2 = j2-i2+1;
			// check if maximal interaction length exceeded (seq2)
			if ( w2 > energy.getAccessibility2().getMaxLength() ) {
				// skip all remaining i2, since only increasing width
				break;
			}
			// check if widths' combination possible
			if ( w2 > maxW2 || w1 > getMaxInteractionWidth( w2, energy.getMaxInternalLoopSize2()) )
			{
				// combination not possible
				hybridE_pq(i1,i2) = E_INF;
				continue;
			}
			// check if left boundary (i1,i2) is complementary
			if (!energy.areComplementary(i1+i1offset,i2+i2offset)) {
				// interaction not possible
				hybridE_pq(i1,i2) = E_INF;
				continue;
			}

			// check if this cell is to be computed (!=E_INF)
			if( hybridE_pq(i1,i2) < E_INF) {

				// compute entry

				// either interaction initiation (i1==j2)
				// or full internal loop energy (nothing between i and j)
				curMinE = energy.getInterLoopE(i1+i1offset,j1+i1offset,i2+i2offset,j2+i2offset)
						+ (i1<j1 ? hybridE_pq(j1,j2) : 0.0 ) ;

				// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
				if (w1 > 2 && w2 > 2) {
					for (k1=std::min(j1-1,i1+energy.getMaxInternalLoopSize1()+1); k1>i1; k1--) {
					for (k2=std::min(j2-1,i2+energy.getMaxInternalLoopSize2()+1); k2>i2; k2--) {
						// check if (k1,k2) are valid left boundary
						if (hybridE_pq(k1,k2) < E_INF) {
							curMinE = std::min( curMinE,
									(energy.getInterLoopE(i1+i1offset,k1+i1offset,i2+i2offset,k2+i2offset)
											+ hybridE_pq(k1,k2) )
									);
						}
					}
					}
				}
				// store value
				hybridE_pq(i1,i2) = curMinE;
				// update mfe if needed
				updateMfe( i1,j1,i2,j2, hybridE_pq(i1,i2) );
				continue;
			}
		}
	}

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2d::
traceBack( Interaction & interaction )
{
	// check if something to trace
	if (interaction.basePairs.size() < 2) {
		return;
	}

#ifdef NDEBUG           /* required by ANSI standard */
	// no check
#else
	// sanity checks
	if ( ! interaction.isValid() ) {
		throw std::runtime_error("PredictorMfe2d::traceBack() : given interaction not valid");
	}
	if ( interaction.basePairs.size() != 2 ) {
		throw std::runtime_error("PredictorMfe2d::traceBack() : given interaction does not contain boundaries only");
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
			i2 = energy.getAccessibility2().getReversedIndex(interaction.basePairs.at(0).second - i2offset),
			j2 = energy.getAccessibility2().getReversedIndex(interaction.basePairs.at(1).second - i2offset);


	// TODO RESTRICT RECOMPUTATION AND TRACEBACK TO I1-J1 I2-J2

	// fill matrix and store best interaction
	fillHybridE( j1, j2 );


	// the currently traced value for i1-j1, i2-j2
	E_type curE = hybridE_pq(i1,i2);

	// trace back
	do {
		// check if just internal loop
		if (curE == (energy.getInterLoopE(i1+i1offset,j1+i1offset,i2+i2offset,j2+i2offset)
				+ (i1<j1 ? hybridE_pq(j1,j2) : 0.0 )) )
		{
			break;
		}
		// check all interval splits
		if ( (j1-i1) > 1 && (j2-i2) > 1) {
			// temp variables
			size_t k1,k2;
			bool traceNotFound = true;
			// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
			for (k1=std::min(j1-1,i1+energy.getMaxInternalLoopSize1()+1); traceNotFound && k1>i1; k1--) {
			for (k2=std::min(j2-1,i2+energy.getMaxInternalLoopSize2()+1); traceNotFound && k2>i2; k2--) {
				// check if (k1,k2) are valid left boundary
				if (hybridE_pq(k1,k2) < E_INF) {
					if (curE == (energy.getInterLoopE(i1+i1offset,k1+i1offset,i2+i2offset,k2+i2offset)
							+ hybridE_pq(k1,k2)) )
					{
						// stop searching
						traceNotFound = false;
						// store splitting base pair
						interaction.addInteraction( k1+i1offset, energy.getAccessibility2().getReversedIndex(k2+i2offset) );
						// trace right part of split
						i1=k1;
						i2=k2;
						curE = hybridE_pq(i1,i2);
					}
				}
			}
			}
		}
	// do until only right boundary is left over
	} while( i1 != j1 );

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

////////////////////////////////////////////////////////////////////////////

