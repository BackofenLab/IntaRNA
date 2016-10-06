
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
predict( const IndexRange & r1
		, const IndexRange & r2 )
{
#ifdef NDEBUG
	// no check
#else
	// check indices (both regions ascending due to reversing of seq2)
	if (!(r1.isAscending() && r2.isAscending()) )
		throw std::runtime_error("PredictorMaxProb::predict("+toString(r1)+","+toString(r2)+") is not sane");
#endif

	// clear data
	clear();

	// resize matrix
	hybridZ.resize( std::min( energy.getAccessibility1().getSequence().size()
						, (r1.to==RnaSequence::lastPos?energy.getAccessibility1().getSequence().size()-1:r1.to)-r1.from+1 )
				, std::min( energy.getAccessibility2().getSequence().size()
						, (r2.to==RnaSequence::lastPos?energy.getAccessibility2().getSequence().size()-1:r2.to)-r2.from+1 ) );

	i1offset = r1.from;
	i2offset = r2.from;

	size_t debug_count_cells_null=0, debug_count_cells_nonNull = 0, debug_cellNumber=0;

	size_t maxWidthFori1i2 = 0;

	bool i1blocked, i1or2blocked;
	// initialize 3rd and 4th dimension of the matrix
	for (size_t i1=0; i1<hybridZ.size1(); i1++) {
		// check if i1 is blocked for interaction
		i1blocked =
		// - shows ambiguous nucleotide encoding
			energy.getAccessibility1().getSequence().isAmbiguous( i1 )
		// - is blocked by an accessibility constraint
			|| energy.getAccessibility1().getAccConstraint().isBlocked(i1);
	for (size_t i2=0; i2<hybridZ.size2(); i2++) {
		// check whether i1 or i2 is blocked for interaction
		i1or2blocked = i1blocked
		// - shows ambiguous nucleotide encoding
			|| energy.getAccessibility2().getSequence().isAmbiguous( i2 )
		// - is blocked by an accessibility constraint
			|| energy.getAccessibility2().getAccConstraint().isBlocked(i2);

		if (hybridZ.size1()-i1 < hybridZ.size2()-i2) {
			maxWidthFori1i2 = getMaxInteractionWidth( hybridZ.size1()-i1, energy.getMaxInternalLoopSize1() );
		} else {
			maxWidthFori1i2 = getMaxInteractionWidth( hybridZ.size2()-i2, energy.getMaxInternalLoopSize2() );
		}

		debug_cellNumber =
				/*w1 = */ std::min(energy.getAccessibility1().getMaxLength(), std::min( hybridZ.size1()-i1, maxWidthFori1i2) )
			*	/*w2 = */ std::min(energy.getAccessibility2().getMaxLength(), std::min( hybridZ.size2()-i2, maxWidthFori1i2) );

		// check if i1 and i2 are not blocked and can form a base pair
		if ( ! i1or2blocked
			&& RnaSequence::areComplementary(
				  energy.getAccessibility1().getSequence()
				, energy.getAccessibility2().getSequence()
				, i1+i1offset, i2+i2offset ))
		{
			// create new 2d matrix for different interaction site widths
			hybridZ(i1,i2) = new E2dMatrix(
				/*w1 = */ std::min(energy.getAccessibility1().getMaxLength(), std::min( hybridZ.size1()-i1, maxWidthFori1i2) ),
				/*w2 = */ std::min(energy.getAccessibility2().getMaxLength(), std::min( hybridZ.size2()-i2, maxWidthFori1i2) ));

			debug_count_cells_nonNull += debug_cellNumber;
		} else {
			// reduce memory consumption and avoid computation for this start index combination
			hybridZ(i1,i2) = NULL;
			debug_count_cells_null += debug_cellNumber;
		}
	}
	}

	std::cerr <<"#DEBUG: init 4d matrix : "<<debug_count_cells_nonNull <<" to be filled ("
				<<((double)debug_count_cells_nonNull/(double)(debug_count_cells_nonNull+debug_count_cells_null))
				<<"%) and "<<debug_count_cells_null <<" not allocated" <<std::endl;

	// fill matrix
	fillHybridZ( energy );

	// maximal probability is
	// double maxProb = (double)maxProbInteraction.getEnergy() / Z;

	// convert the partition function into an ensemble energy
	maxProbInteraction.energy = ( - energy.getRT() * std::log( maxProbInteraction.energy ) );

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

	// current Z value
	E_type curZ = 0;
	// iterate increasingly over all window sizes w1 (seq1) and w2 (seq2)
	for (w1=0; w1<energy.getAccessibility1().getMaxLength(); w1++) {
	for (w2=0; w2<energy.getAccessibility2().getMaxLength(); w2++) {
		// iterate over all window starts i1 (seq1) and i2 (seq2)
		// TODO PARALLELIZE THIS DOUBLE LOOP ?!
		for (i1=0; i1+w1<hybridZ.size1(); i1++) {
		for (i2=0; i2+w2<hybridZ.size2(); i2++) {
			// check if left boundary is complementary
			// and widths are possible
			if (hybridZ(i1,i2) == NULL || hybridZ(i1,i2)->size1()<=w1 || hybridZ(i1,i2)->size2()<=w2) {
				// interaction not possible: nothing to do
				continue;
			}
			// check if interaction exceeds possible with due to max-loop-length
			if ( getMaxInteractionWidth( 1+w1, energy.getMaxInternalLoopSize1() ) < w2
				|| getMaxInteractionWidth( 1+w2, energy.getMaxInternalLoopSize2() ) < w1)
			{
				// ignore this entry
				(*hybridZ(i1,i2))(w1,w2) = 0;
				continue;
			}
			// check if right boundary is complementary
			if (hybridZ(j1,j2) == NULL) {
				// not complementary -> ignore this entry
				(*hybridZ(i1,i2))(w1,w2) = 0;
				continue;
			}

			// get window ends j (seq1) and l (seq2)
			j1=i1+w1;
			j2=i2+w2;
			// compute entry
			curZ = 0;

			// get full internal loop energy (nothing between i and j)
			// if allowed distance between i and j
			if ( (w1+1) <= energy.getMaxInternalLoopSize1() && (w2+1) <= energy.getMaxInternalLoopSize2()) {
				curZ += energy.getBoltzmannWeight(energy.getInterLoopE(i1+i1offset,j1+i1offset,i2+i2offset,j2+i2offset))
					* energy.getBoltzmannWeight(energy.getInterLoopE(j1+i1offset,j1+i1offset,j2+i2offset,j2+i2offset));
			}

			if (w1 > 1 && w2 > 1) {
				// sum all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
				for (k1=std::min(j1-1,i1+energy.getMaxInternalLoopSize1()+1); k1>i1; k1--) {
				for (k2=std::min(j2-1,i2+energy.getMaxInternalLoopSize2()+1); k2>i2; k2--) {
					// check if (k1,k2) are complementary
					if (hybridZ(k1,k2) != NULL) {
						curZ += energy.getBoltzmannWeight(energy.getInterLoopE(i1+i1offset,k1+i1offset,i2+i2offset,k2+i2offset))
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

	//////////  SECOND ROUND : COMPUTE FINAL ENERGIES AND MFE  ////////////

	// initialize max prob interaction for updates
	initMaxProbInteraction();
	// reset overall partition function
	Z = (E_type)0.0;

	// iterate increasingly over all window sizes w1 (seq1) and w2 (seq2)
	for (w1=0; w1<energy.getAccessibility1().getMaxLength(); w1++) {
	for (w2=0; w2<energy.getAccessibility2().getMaxLength(); w2++) {
		// iterate over all window starts i1 (seq1) and i2 (seq2)
		for (i1=0; i1+w1<hybridZ.size1(); i1++) {
		for (i2=0; i2+w2<hybridZ.size2(); i2++) {
			// check if left boundary is complementary
			// and widths are possible
			if (hybridZ(i1,i2) == NULL || hybridZ(i1,i2)->size1()<=w1 || hybridZ(i1,i2)->size2()<=w2) {
				// interaction not possible: nothing to do
				continue;
			}
			// check if reasonable entry
			if ((*hybridZ(i1,i2))(w1,w2) < E_INF) {
				// get window ends j (seq1) and l (seq2)
				j1=i1+w1;
				j2=i2+w2;

				// update if needed
				// TODO : intaRNA-1 uses dangling contribs only, if dangling base is in ED region, ie additional cases used
				// TODO AU-end penalty missing
				updateMaxProbInteraction( i1,j1,i2,j2
						, (*hybridZ(i1,i2))(w1,w2)
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
	maxProbInteraction.energy = 0.0;
	// reset boundary base pairs
	maxProbInteraction.r1.from = RnaSequence::lastPos;
	maxProbInteraction.r1.to = RnaSequence::lastPos;
	maxProbInteraction.r2.from = RnaSequence::lastPos;
	maxProbInteraction.r2.to = RnaSequence::lastPos;
}
////////////////////////////////////////////////////////////////////////////

void
PredictorMaxProb::
updateMaxProbInteraction( const size_t i1, const size_t j1
		, const size_t i2, const size_t j2
		, const E_type hybridZ )
{
//				std::cerr <<"#DEBUG : Z( "<<i1<<"-"<<j1<<", "<<i2<<"-"<<j2<<" ) = "
//						<<curZ
//						<<" = " <<(eH + eE + eD)
//						<<std::endl;

	// add Boltzmann weights of all penalties
	E_type curZ = hybridZ * energy.getBoltzmannWeight( energy.getE(i1,j1,i2,j2,0.0) );

	// update overall partition function
	Z += (double)curZ;

	if (curZ > maxProbInteraction.energy) {
		// store new global min
		maxProbInteraction.energy = curZ;
		// store interaction boundaries
		// left
		maxProbInteraction.r1.from = i1+i1offset;
		maxProbInteraction.r2.from = energy.getAccessibility2().getReversedIndex(i2+i2offset);
		// right
		maxProbInteraction.r1.to = j1+i1offset;
		maxProbInteraction.r2.to = energy.getAccessibility2().getReversedIndex(j2+i2offset);
	}
}


////////////////////////////////////////////////////////////////////////////

