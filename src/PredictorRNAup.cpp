
#include "PredictorRNAup.h"

////////////////////////////////////////////////////////////////////////////

PredictorRNAup::
PredictorRNAup( const InteractionEnergy & energy, const OutputHandler & output )
 : Predictor(energy,output)
	, hybridE( energy.getAccessibility1().getSequence().size()
			, energy.getAccessibility2().getSequence().size() )
	, hybridEmin(E_MAX)
{

	// initialize 3rd and 4th dimension of the matrix
	for (size_t i1=0; i1<energy.getAccessibility1().getSequence().size(); i1++) {
	for (size_t i2=0; i2<energy.getAccessibility2().getSequence().size(); i2++) {
		// check if i and k can form a base pair
		if (RnaSequence::areComplementary(
				  energy.getAccessibility1().getSequence()
				, energy.getAccessibility2().getSequence()
				, i1, i2 ))
		{
			// create new 2d matrix for different interaction site widths
			hybridE(i1,i2) = new E2dMatrix(
				/*w1 = */ std::min(energy.getAccessibility1().getMaxLength(), energy.getAccessibility1().getSequence().size()-i1),
				/*w2 = */ std::min(energy.getAccessibility2().getMaxLength(), energy.getAccessibility2().getSequence().size()-i2));
		} else {
			// reduce memory consumption
			hybridE(i1,i2) = NULL;
		}
	}
	}

	// fill matrix
	fillHybridE( energy );

}


////////////////////////////////////////////////////////////////////////////

PredictorRNAup::
~PredictorRNAup()
{
	// delete 3rd and 4th dimension of the matrix
	for (E4dMatrix::iterator1 ijEntry = hybridE.begin1(); ijEntry != hybridE.end1(); ijEntry++) {
		if (*ijEntry != NULL) {
			// delete 2d matrix for current ij
			delete (*ijEntry);
			*ijEntry = NULL;
		}
	}
}


////////////////////////////////////////////////////////////////////////////

void
PredictorRNAup::
fillHybridE( const InteractionEnergy & energy )
{


	//////////  FIRST ROUND : COMPUTE HYBRIDIZATION ENERGIES ONLY  ////////////

	// current minimal value
	E_type curMinE = E_MAX;
	// iterate increasingly over all window sizes w1 (seq1) and w2 (seq2)
	for (size_t w1=0; w1<energy.getAccessibility1().getMaxLength(); w1++) {
	for (size_t w2=0; w2<energy.getAccessibility2().getMaxLength(); w2++) {
		// iterate over all window starts i (seq1) and k (seq2)
		for (size_t i1=0; i1+w1<hybridE.size1(); i1++) {
		for (size_t i2=0; i2+w2<hybridE.size2(); i2++) {
			// check if left boundary is complementary
			if (hybridE(i1,i2) == NULL) {
				// nothing to do
				continue;
			}
			// get window ends j (seq1) and l (seq2)
			size_t j1=i1+w1;
			size_t j2=i2+w2;
			// check if right boundary is complementary
			if (hybridE(j1,j2) == NULL) {
				// not complementary -> ignore this entry
				(*hybridE(i1,i2))(w1,w2) = E_MAX;
			} else {
				// compute entry
				// get full internal loop energy (nothing between i and j)
				curMinE = E_MAX;
				if ((i1==j1 && i2==j2) || (i1<j1 && i2<j2) ) {
					curMinE = energy.getInterLoopE(i1,j1,i2,j2);
				}
				if (w1 > 1 && w2 > 1) {
					// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
					for (size_t k1=std::min(j1-1,i1+energy.getMaxInternalLoopSize1()+1); k1>i1; k1--) {
					for (size_t k2=std::min(j2-1,i2+energy.getMaxInternalLoopSize2()+1); k2>i2; k2--) {
						// check if (k1,k2) are complementary
						if (hybridE(k1,k2) != NULL) {
							curMinE = std::min( curMinE,
									(energy.getInterLoopE(i1,k1,i2,k2) + (*hybridE(k1,k2))(j1-k1,j2-k2))
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

	//////////  SECOND ROUND : COMPUTE FINAL ENERGIES  ////////////


	// global minimal value
	hybridEmin = E_MAX;

	// iterate increasingly over all window sizes w1 (seq1) and w2 (seq2)
	for (size_t w1=0; w1<energy.getAccessibility1().getMaxLength(); w1++) {
	for (size_t w2=0; w2<energy.getAccessibility2().getMaxLength(); w2++) {
		// iterate over all window starts i (seq1) and k (seq2)
		for (size_t i1=0; i1+w1<hybridE.size1(); i1++) {
		for (size_t i2=0; i2+w2<hybridE.size2(); i2++) {
			// check if left boundary is complementary
			if (hybridE(i1,i2) == NULL) {
				// nothing to do
				continue;
			}
			// get window ends j (seq1) and l (seq2)
			size_t j1=i1+w1;
			size_t j2=i2+w2;
			// check if reasonable entry
			if ((*hybridE(i1,i2))(w1,w2) < E_MAX) {
				// correct energy with accessibility
				(*hybridE(i1,i2))(w1,w2)
						+= energy.getAccessibility1().getED(i1,j1);
						+  energy.getAccessibility2().getED(i2,j2);

				// update global min
				hybridEmin = std::min(hybridEmin, (*hybridE(i1,i2))(w1,w2));
			}

		} // i2
		} // i1
	} // w2
	} // w1


	std::cerr <<"\n#DEBUG : hybridEmin = " << hybridEmin <<"\n" <<std::endl;
}


////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////

