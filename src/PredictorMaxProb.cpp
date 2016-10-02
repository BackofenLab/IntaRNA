
#include "PredictorMaxProb.h"

#include <stdexcept>


////////////////////////////////////////////////////////////////////////////

PredictorMaxProb::
PredictorMaxProb( const InteractionEnergy & energy, OutputHandler & output )
 : Predictor(energy,output)
	, hybridZ( 0,0 )
	, Z(0.0)
	, maxProbInteraction(energy.getAccessibility1().getSequence()
			,energy.getAccessibility2().getAccessibilityOrigin().getSequence())
	, i1offset(0)
	, i2offset(0)
{
}


////////////////////////////////////////////////////////////////////////////

PredictorMaxProb::
~PredictorMaxProb()
{
	// clean up
	this->clear();
}


////////////////////////////////////////////////////////////////////////////

void
PredictorMaxProb::
predict( const size_t i1, const size_t j1
				, const size_t i2, const size_t j2)
{
	// TODO check indices

	// clear data
	clear();

	// resize matrix
	hybridZ.resize( std::min( energy.getAccessibility1().getSequence().size()
						, (j1==RnaSequence::lastPos?energy.getAccessibility1().getSequence().size()-1:j1)-i1+1 )
				, std::min( energy.getAccessibility2().getSequence().size()
						, (j2==RnaSequence::lastPos?energy.getAccessibility2().getSequence().size()-1:j2)-i2+1 ) );

	i1offset = i1;
	i2offset = i2;

	// initialize 3rd and 4th dimension of the matrix
	for (size_t i1=0; i1<hybridZ.size1(); i1++) {
	for (size_t i2=0; i2<hybridZ.size2(); i2++) {
		// check if i and k can form a base pair
		if (RnaSequence::areComplementary(
				  energy.getAccessibility1().getSequence()
				, energy.getAccessibility2().getSequence()
				, i1+i1offset, i2+i2offset ))
		{
			// create new 2d matrix for different interaction site widths
			hybridZ(i1,i2) = new E2dMatrix(
				/*w1 = */ std::min(energy.getAccessibility1().getMaxLength(), hybridZ.size1()-i1),
				/*w2 = */ std::min(energy.getAccessibility2().getMaxLength(), hybridZ.size2()-i2));
		} else {
			// reduce memory consumption
			hybridZ(i1,i2) = NULL;
		}
	}
	}

	// fill matrix
	fillHybridZ( energy );

	// maximal probability is
	// double maxProb = (double)maxProbInteraction.getEnergy() / Z;

	// convert the partition function into an ensemble energy
	maxProbInteraction.setEnergy( - energy.getRT() * std::log( maxProbInteraction.getEnergy() ) );

	// report interaction site with maximal probability
	output.add( maxProbInteraction );
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMaxProb::
clear()
{
	// delete 3rd and 4th dimension of the matrix
	for (E4dMatrix::iterator1 ijEntry = hybridZ.begin1(); ijEntry != hybridZ.end1(); ijEntry++) {
		if (*ijEntry != NULL) {
			// delete 2d matrix for current ij
			delete (*ijEntry);
			*ijEntry = NULL;
		}
	}
	// clear matrix, free data
	hybridZ.clear();
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMaxProb::
fillHybridZ( const InteractionEnergy & energy )
{

	// global vars to avoid reallocation
	size_t i1,i2,j1,j2,w1,w2,k1,k2;

	//////////  FIRST ROUND : COMPUTE HYBRIDIZATION ENERGIES ONLY  ////////////

	// current minimal value
	E_type curZ = E_MAX;
	// iterate increasingly over all window sizes w1 (seq1) and w2 (seq2)
	for (w1=0; w1<energy.getAccessibility1().getMaxLength(); w1++) {
	for (w2=0; w2<energy.getAccessibility2().getMaxLength(); w2++) {
		// iterate over all window starts i (seq1) and k (seq2)
		// TODO PARALLELIZE THIS DOUBLE LOOP ?!
		for (i1=0; i1+w1<hybridZ.size1(); i1++) {
		for (i2=0; i2+w2<hybridZ.size2(); i2++) {
			// check if left boundary is complementary
			if (hybridZ(i1,i2) == NULL) {
				// nothing to do
				continue;
			}
			// get window ends j (seq1) and l (seq2)
			j1=i1+w1;
			j2=i2+w2;
			// check if right boundary is complementary
			if (hybridZ(j1,j2) == NULL) {
				// not complementary -> ignore this entry
				(*hybridZ(i1,i2))(w1,w2) = E_MAX;
			} else {
				// compute entry
				curZ = 0;
				// get full internal loop energy (nothing between i and j)
				// if allowed distance between i and j
				if ( (w1+1) <= energy.getMaxInternalLoopSize1() && (w2+1) <= energy.getMaxInternalLoopSize2()) {
					curZ += getBoltzmannWeight(energy.getInterLoopE(i1+i1offset,j1+i1offset,i2+i2offset,j2+i2offset))
						* getBoltzmannWeight(energy.getInterLoopE(j1+i1offset,j1+i1offset,j2+i2offset,j2+i2offset));
				}

				if (w1 > 1 && w2 > 1) {
					// sum all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
					for (k1=std::min(j1-1,i1+energy.getMaxInternalLoopSize1()+1); k1>i1; k1--) {
					for (k2=std::min(j2-1,i2+energy.getMaxInternalLoopSize2()+1); k2>i2; k2--) {
						// check if (k1,k2) are complementary
						if (hybridZ(k1,k2) != NULL) {
							curZ += getBoltzmannWeight(energy.getInterLoopE(i1+i1offset,k1+i1offset,i2+i2offset,k2+i2offset))
									* ((*hybridZ(k1,k2))(j1-k1,j2-k2));
						}
					}
					}
				}
				// store value
				(*hybridZ(i1,i2))(w1,w2) = curZ;
			}
		}
		}
	}
	}

	//////////  SECOND ROUND : COMPUTE FINAL ENERGIES AND MFE  ////////////

	// initialize max prob interaction for updates
	initMaxProbInteraction();
	// reset overall partition function
	Z = (E_type)0.0;

	// iterate increasingly over all window sizes w1 (seq1) and w2 (seq2)
	for (w1=0; w1<energy.getAccessibility1().getMaxLength(); w1++) {
	for (w2=0; w2<energy.getAccessibility2().getMaxLength(); w2++) {
		// iterate over all window starts i (seq1) and k (seq2)
		for (i1=0; i1+w1<hybridZ.size1(); i1++) {
		for (i2=0; i2+w2<hybridZ.size2(); i2++) {
			// check if left boundary is complementary
			if (hybridZ(i1,i2) == NULL) {
				// nothing to do
				continue;
			}
			// check if reasonable entry
			if ((*hybridZ(i1,i2))(w1,w2) < E_MAX) {
				// get window ends j (seq1) and l (seq2)
				j1=i1+w1;
				j2=i2+w2;

				// update if needed
				updateMaxProbInteraction( i1,j1,i2,j2
						, (*hybridZ(i1,i2))(w1,w2)
						 * getBoltzmannWeight(energy.getDanglingLeft(i1+i1offset,i2+i2offset))
						 * getBoltzmannWeight(energy.getDanglingRight(j1+i1offset,j2+i2offset))
						 * getBoltzmannWeight(energy.getAccessibility1().getED(i1+i1offset,j1+i1offset))
						 * getBoltzmannWeight(energy.getAccessibility2().getED(i2+i2offset,j2+i2offset))
						);

			}

		} // i2
		} // i1
	} // w2
	} // w1

}


////////////////////////////////////////////////////////////////////////////

void
PredictorMaxProb::
initMaxProbInteraction()
{
	// initialize max prob interaction (partition function value)
	maxProbInteraction.setEnergy(0.0);
	// ensure it holds only the boundary
	if (maxProbInteraction.getBasePairs().size()!=2) {
		maxProbInteraction.getBasePairs().resize(2);
	}
	// reset boundary base pairs
	maxProbInteraction.getBasePairs()[0].first = RnaSequence::lastPos;
	maxProbInteraction.getBasePairs()[0].second = RnaSequence::lastPos;
	maxProbInteraction.getBasePairs()[1].first = RnaSequence::lastPos;
	maxProbInteraction.getBasePairs()[1].second = RnaSequence::lastPos;
}
////////////////////////////////////////////////////////////////////////////

void
PredictorMaxProb::
updateMaxProbInteraction( const size_t i1, const size_t j1
		, const size_t i2, const size_t j2
		, const E_type curZ )
{
//				std::cerr <<"#DEBUG : Z( "<<i1<<"-"<<j1<<", "<<i2<<"-"<<j2<<" ) = "
//						<<curZ
//						<<" = " <<(eH + eE + eD)
//						<<std::endl;

	// update overall partition function
	Z += (double)curZ;

	if (curZ > maxProbInteraction.getEnergy()) {
		// store new global min
		maxProbInteraction.setEnergy(curZ);
		// store interaction boundaries
		// left
		maxProbInteraction.getBasePairs()[0].first = i1+i1offset;
		maxProbInteraction.getBasePairs()[0].second = energy.getAccessibility2().getReversedIndex(i2+i2offset);
		// right
		maxProbInteraction.getBasePairs()[1].first = j1+i1offset;
		maxProbInteraction.getBasePairs()[1].second = energy.getAccessibility2().getReversedIndex(j2+i2offset);
	}
}

////////////////////////////////////////////////////////////////////////////

E_type
PredictorMaxProb::
getBoltzmannWeight( const E_type e )
{
	// TODO can be optimized when using exp-energies from VRNA
	return std::exp( - e / energy.getRT() );
}


////////////////////////////////////////////////////////////////////////////

