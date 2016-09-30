
#include "PredictorRNAup.h"

#include <stdexcept>


////////////////////////////////////////////////////////////////////////////

PredictorMfeRNAup::
PredictorMfeRNAup( const InteractionEnergy & energy, OutputHandler & output )
 : Predictor(energy,output)
	, hybridE( 0,0 )
	, i1offset(0)
	, i2offset(0)
	, mfeInteraction(energy.getAccessibility1().getSequence()
							,energy.getAccessibility2().getAccessibilityOrigin().getSequence())
{

}


////////////////////////////////////////////////////////////////////////////

PredictorMfeRNAup::
~PredictorMfeRNAup()
{
	// clean up
	this->clear();
}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfeRNAup::
predict( const size_t i1, const size_t j1
		, const size_t i2, const size_t j2)
{
	// TODO check indices

	// clear data
	this->clear();

	// resize matrix
	hybridE.resize( std::min( energy.getAccessibility1().getSequence().size()
						, (j1==RnaSequence::lastPos?energy.getAccessibility1().getSequence().size()-1:j1)-i1+1 )
				, std::min( energy.getAccessibility2().getSequence().size()
						, (j2==RnaSequence::lastPos?energy.getAccessibility2().getSequence().size()-1:j2)-i2+1 ) );
//	hybridE.resize( energy.getAccessibility1().getSequence().size()
//				, energy.getAccessibility2().getSequence().size() );

	i1offset = i1;
	i2offset = i2;

	// initialize 3rd and 4th dimension of the matrix
	for (size_t i1=0; i1<hybridE.size1(); i1++) {
	for (size_t i2=0; i2<hybridE.size2(); i2++) {
		// check if i and k can form a base pair
		if (RnaSequence::areComplementary(
				  energy.getAccessibility1().getSequence()
				, energy.getAccessibility2().getSequence()
				, i1+i1offset, i2+i2offset ))
//				, i1, i2 ))
		{
			// create new 2d matrix for different interaction site widths
			hybridE(i1,i2) = new E2dMatrix(
				/*w1 = */ std::min(energy.getAccessibility1().getMaxLength(), hybridE.size1()-i1),
				/*w2 = */ std::min(energy.getAccessibility2().getMaxLength(), hybridE.size2()-i2));
//				/*w1 = */ std::min(energy.getAccessibility1().getMaxLength(), energy.getAccessibility1().getSequence().size()-i1),
//				/*w2 = */ std::min(energy.getAccessibility2().getMaxLength(), energy.getAccessibility2().getSequence().size()-i2));
		} else {
			// reduce memory consumption
			hybridE(i1,i2) = NULL;
		}
	}
	}

	// fill matrix and find mfe
	fillHybridE( energy );

	// get mfe interaction details
	traceBack( mfeInteraction );

	// report mfe interaction
	output.add( mfeInteraction );

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeRNAup::
clear()
{
	// delete 3rd and 4th dimension of the matrix
	for (E4dMatrix::iterator1 ijEntry = hybridE.begin1(); ijEntry != hybridE.end1(); ijEntry++) {
		if (*ijEntry != NULL) {
			// delete 2d matrix for current ij
			delete (*ijEntry);
			*ijEntry = NULL;
		}
	}
	// clear matrix, free data
	hybridE.clear();
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeRNAup::
fillHybridE( const InteractionEnergy & energy )
{

	// global vars to avoid reallocation
	size_t i1,i2,j1,j2,w1,w2,k1,k2;

	//////////  FIRST ROUND : COMPUTE HYBRIDIZATION ENERGIES ONLY  ////////////

	// current minimal value
	E_type curMinE = E_MAX;
	// iterate increasingly over all window sizes w1 (seq1) and w2 (seq2)
	for (w1=0; w1<energy.getAccessibility1().getMaxLength(); w1++) {
	for (w2=0; w2<energy.getAccessibility2().getMaxLength(); w2++) {
		// iterate over all window starts i (seq1) and k (seq2)
		// TODO PARALLELIZE THIS DOUBLE LOOP ?!
		for (i1=0; i1+w1<hybridE.size1(); i1++) {
		for (i2=0; i2+w2<hybridE.size2(); i2++) {
			// check if left boundary is complementary
			if (hybridE(i1,i2) == NULL) {
				// nothing to do
				continue;
			}
			// get window ends j (seq1) and l (seq2)
			j1=i1+w1;
			j2=i2+w2;
			// check if right boundary is complementary
			if (hybridE(j1,j2) == NULL) {
				// not complementary -> ignore this entry
				(*hybridE(i1,i2))(w1,w2) = E_MAX;
			} else {
				// compute entry
				// get full internal loop energy (nothing between i and j)
				curMinE = energy.getInterLoopE(i1+i1offset,j1+i1offset,i2+i2offset,j2+i2offset) + energy.getInterLoopE(j1+i1offset,j1+i1offset,j2+i2offset,j2+i2offset);
//				curMinE = (energy.getInterLoopE(i1,j1,i2,j2) + energy.getInterLoopE(j1,j1,j2,j2));


				if (w1 > 1 && w2 > 1) {
					// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
					for (k1=std::min(j1-1,i1+energy.getMaxInternalLoopSize1()+1); k1>i1; k1--) {
					for (k2=std::min(j2-1,i2+energy.getMaxInternalLoopSize2()+1); k2>i2; k2--) {
						// check if (k1,k2) are complementary
						if (hybridE(k1,k2) != NULL) {
							curMinE = std::min( curMinE,
									(energy.getInterLoopE(i1+i1offset,k1+i1offset,i2+i2offset,k2+i2offset) + (*hybridE(k1,k2))(j1-k1,j2-k2))
//									(energy.getInterLoopE(i1,k1,i2,k2) + (*hybridE(k1,k2))(j1-k1,j2-k2))
									);
						}
					}
					}
				}
				// store value
				(*hybridE(i1,i2))(w1,w2) = curMinE;
			}
		}
		}
	}
	}

	//////////  SECOND ROUND : COMPUTE FINAL ENERGIES AND MFE  ////////////

	// set initial global minimal value
	initMfe();

	// iterate increasingly over all window sizes w1 (seq1) and w2 (seq2)
	for (w1=0; w1<energy.getAccessibility1().getMaxLength(); w1++) {
	for (w2=0; w2<energy.getAccessibility2().getMaxLength(); w2++) {
		// iterate over all window starts i (seq1) and k (seq2)
		for (i1=0; i1+w1<hybridE.size1(); i1++) {
		for (i2=0; i2+w2<hybridE.size2(); i2++) {
			// check if left boundary is complementary
			if (hybridE(i1,i2) == NULL) {
				// nothing to do
				continue;
			}
			// check if reasonable entry
			if ((*hybridE(i1,i2))(w1,w2) < E_MAX) {
				// get window ends j (seq1) and l (seq2)
				j1=i1+w1;
				j2=i2+w2;

				// update mfe if needed
				updateMfe( i1,j1,i2,j2
						, (*hybridE(i1,i2))(w1,w2)
						, energy.getDanglingLeft(i1+i1offset,i2+i2offset)
							+ energy.getDanglingRight(j1+i1offset,j2+i2offset)
						, energy.getAccessibility1().getED(i1+i1offset,j1+i1offset)
							+ energy.getAccessibility2().getED(i2+i2offset,j2+i2offset)
//						, energy.getDanglingLeft(i1,i2)
//							+ energy.getDanglingRight(j1,j2)
//						, energy.getAccessibility1().getED(i1,j1)
//							+ energy.getAccessibility2().getED(i2,j2)
						);

			}

		} // i2
		} // i1
	} // w2
	} // w1

}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfeRNAup::
initMfe()
{
	// delete stored global E min interaction
	mfeInteraction.clear();
	// initialize global E minimum
	mfeInteraction.setEnergy(E_MAX);
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeRNAup::
updateMfe( const size_t i1, const size_t j1
		, const size_t i2, const size_t j2
		, const E_type eH, const E_type eE, const E_type eD )
{
//				std::cerr <<"#DEBUG : energy( "<<i1<<"-"<<j1<<", "<<i2<<"-"<<j2<<" ) = "
//						<<eH<<" + "<<eE <<" + "<<eD
//						<<" = " <<(eH + eE + eD)
//						<<std::endl;

	if (eH+eE+eD < mfeInteraction.getEnergy()) {
		// store new global min
		mfeInteraction.setEnergy(eH+eE+eD);
		// store interaction boundaries
		mfeInteraction.clear();
		mfeInteraction.addInteraction(i1+i1offset, energy.getAccessibility2().getReversedIndex(i2+i2offset));
		mfeInteraction.addInteraction(j1+i1offset, energy.getAccessibility2().getReversedIndex(j2+i2offset));
//		mfeInteraction.addInteraction(i1, energy.getAccessibility2().getReversedIndex(i2));
//		mfeInteraction.addInteraction(j1, energy.getAccessibility2().getReversedIndex(j2));
	}
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeRNAup::
traceBack( Interaction & interaction ) const
{
	// check if something to trace
	if (interaction.getBasePairs().size() < 2) {
		return;
	}

#ifdef NDEBUG           /* required by ANSI standard */
	// no check
#else
	// sanity checks
	if ( ! interaction.isValid() ) {
		throw std::runtime_error("PredictorRNAup::traceBack() : given interaction not valid");
	}
	if ( interaction.getBasePairs().size() != 2 ) {
		throw std::runtime_error("PredictorRNAup::traceBack() : given interaction does not contain boundaries only");
	}
#endif

	// ensure sorting
	interaction.sort();
	// get indices in hybridE for boundary base pairs
	size_t	i1 = interaction.getBasePairs().at(0).first - i1offset,
			j1 = interaction.getBasePairs().at(1).first - i1offset,
			i2 = energy.getAccessibility2().getReversedIndex(interaction.getBasePairs().at(0).second - i2offset),
			j2 = energy.getAccessibility2().getReversedIndex(interaction.getBasePairs().at(1).second - i2offset);
//	size_t	i1 = interaction.getBasePairs().at(0).first,
//			j1 = interaction.getBasePairs().at(1).first,
//			i2 = energy.getAccessibility2().getReversedIndex(interaction.getBasePairs().at(0).second),
//			j2 = energy.getAccessibility2().getReversedIndex(interaction.getBasePairs().at(1).second);

	// the currently traced value for i1-j1, i2-j2
	E_type curE = (*hybridE(i1,i2))(j1-i1,j2-i2);

	// trace back
	do {
		// check if just internal loop
		if (curE == (energy.getInterLoopE(i1+i1offset,j1+i1offset,i2+i2offset,j2+i2offset) + energy.getInterLoopE(j1+i1offset,j1+i1offset,j2+i2offset,j2+i2offset))) {
//		if (curE == (energy.getInterLoopE(i1,j1,i2,j2) + energy.getInterLoopE(j1,j1,j2,j2))) {
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
				// check if (k1,k2) are complementary
				if (hybridE(k1,k2) != NULL) {
					if (curE == (energy.getInterLoopE(i1+i1offset,k1+i1offset,i2+i2offset,k2+i2offset) + (*hybridE(k1,k2))(j1-k1,j2-k2)) ) {
//					if (curE == (energy.getInterLoopE(i1,k1,i2,k2) + (*hybridE(k1,k2))(j1-k1,j2-k2)) ) {
						// stop searching
						traceNotFound = false;
						// store splitting base pair
						interaction.addInteraction( k1+i1offset, energy.getAccessibility2().getReversedIndex(k2+i2offset) );
//						interaction.addInteraction( k1, energy.getAccessibility2().getReversedIndex(k2) );
						// trace right part of split
						i1=k1;
						i2=k2;
						curE = (*hybridE(i1,i2))(j1-i1,j2-i2);
					}
				}
			}
			}
		}
	// do until only right boundary is left over
	} while( i1 != j1 );

	// sort final interaction (to make valid) (faster than calling sort())
	if (interaction.getBasePairs().size() > 2) {
		Interaction::PairingVec & bps = interaction.getBasePairs();
		// shift all added base pairs to the front
		for (size_t i=2; i<bps.size(); i++) {
			bps.at(i-1).first = bps.at(i).first;
			bps.at(i-1).second = bps.at(i).second;
		}
		// set last to j1-j2
		bps.rbegin()->first = j1+i1offset;
		bps.rbegin()->second = energy.getAccessibility2().getReversedIndex(j2+i2offset);
//		bps.rbegin()->first = j1;
//		bps.rbegin()->second = energy.getAccessibility2().getReversedIndex(j2);
	}

}

////////////////////////////////////////////////////////////////////////////

