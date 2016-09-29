
#include "PredictorRNAup.h"

////////////////////////////////////////////////////////////////////////////

PredictorRNAup::
PredictorRNAup( const InteractionEnergy & energy, const OutputHandler & output )
 : Predictor(energy,output)
	, hybridE( energy.getAccessibility1().getSequence().size()
			, energy.getAccessibility1().getSequence().size()
			, 0
			, energy.getAccessibility1().getMaxLength() )
{

	// initialize 3rd and 4th dimension of the matrix
	for (E4dMatrix::iterator1 ijEntry = hybridE.begin1(); ijEntry != hybridE.end1(); ijEntry++) {
		// create new 2d matrix for kl of second sequence for current ij
		*ijEntry = new E2dMatrix( energy.getAccessibility2().getSequence().size()
								, energy.getAccessibility2().getSequence().size()
								, 0
								, energy.getAccessibility2().getMaxLength() );
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
		// delete 2d matrix for current ij
		delete (*ijEntry);
		*ijEntry = NULL;
	}
}


////////////////////////////////////////////////////////////////////////////

void
PredictorRNAup::
fillHybridE( const InteractionEnergy & energy )
{
	// iterate increasingly over all window sizes w1 (seq1) and w2 (seq2)
	for (size_t w1=0; w1<energy.getAccessibility1().getMaxLength(); w1++) {
	for (size_t w2=0; w2<energy.getAccessibility2().getMaxLength(); w2++) {
		// iterate over all window starts i (seq1) and k (seq2)
		for (size_t i=0; i+w1<hybridE.size1(); i++) {
		for (size_t k=0; k+w2<hybridE.size2(); k++) {
			// get window ends j (seq1) and l (seq2)
			size_t j=i+w1;
			size_t l=k+w2;
		}
		}
	}
	}
}


////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////

