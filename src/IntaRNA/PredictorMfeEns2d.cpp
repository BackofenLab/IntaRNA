
#include "IntaRNA/PredictorMfeEns2d.h"

#include <stdexcept>

namespace IntaRNA {

////////////////////////////////////////////////////////////////////////////

PredictorMfeEns2d::
PredictorMfeEns2d(
		const InteractionEnergy & energy
		, OutputHandler & output
		, PredictionTracker * predTracker )
 : PredictorMfeEns(energy,output,predTracker)
	, hybridZ( 0,0 )
{
}


////////////////////////////////////////////////////////////////////////////

PredictorMfeEns2d::
~PredictorMfeEns2d()
{
	// clean up
}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2d::
predict( const IndexRange & r1
		, const IndexRange & r2 )
{
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) <<"predicting ensemble mfe interactions in O(n^2) space..."; }
	// measure timing
	TIMED_FUNC_IF(timerObj,VLOG_IS_ON(9));

#if INTARNA_IN_DEBUG_MODE
	// check indices
	if (!(r1.isAscending() && r2.isAscending()) )
		throw std::runtime_error("PredictorMfeEns2d::predict("+toString(r1)+","+toString(r2)+") is not sane");
#endif

	// set index offset
	energy.setOffset1(r1.from);
	energy.setOffset2(r2.from);

	// resize matrix
	hybridZ.resize( std::min( energy.size1()
						, (r1.to==RnaSequence::lastPos?energy.size1()-1:r1.to)-r1.from+1 )
				, std::min( energy.size2()
						, (r2.to==RnaSequence::lastPos?energy.size2()-1:r2.to)-r2.from+1 ) );

	// initialize mfe interaction for updates
	initOptima();
	// initialize overall partition function for updates
	initZ();

	// for all right ends j1
	for (size_t j1 = hybridZ.size1(); j1-- > 0; ) {
		// check if j1 is accessible
		if (!energy.isAccessible1(j1))
			continue;
		// iterate over all right ends j2
		for (size_t j2 = hybridZ.size2(); j2-- > 0; ) {
			// check if j2 is accessible
			if (!energy.isAccessible2(j2))
				continue;
			// check if base pair (j1,j2) possible
			if (!energy.areComplementary( j1, j2 ))
				continue;

			// fill matrix and store best interaction
			fillHybridZ( j1, j2, 0, 0, true );

		}
	}

	// report mfe interaction
	reportOptima();
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2d::
fillHybridZ( const size_t j1, const size_t j2
			, const size_t i1init, const size_t i2init
			, const bool callUpdateZ )
{
	// temporary access
	const OutputConstraint & outConstraint = output.getOutputConstraint();
#if INTARNA_IN_DEBUG_MODE
	if (i1init > j1)
		throw std::runtime_error("PredictorMfeEns2d::fillHybridZ() : i1init > j1 : "+toString(i1init)+" > "+toString(j1));
	if (i2init > j2)
		throw std::runtime_error("PredictorMfeEns2d::fillHybridZ() : i2init > j2 : "+toString(i2init)+" > "+toString(j2));
#endif

	// get minimal start indices heeding max interaction length
	const size_t i1start = std::max(i1init,j1-std::min(j1,energy.getAccessibility1().getMaxLength()+1));
	const size_t i2start = std::max(i2init,j2-std::min(j2,energy.getAccessibility2().getMaxLength()+1));

	// global vars to avoid reallocation
	size_t i1,i2,w1,w2,k1,k2;

	// determine whether or not lonely base pairs are allowed or if we have to
	// ensure a stacking to the right of the left boundary (i1,i2)
	const size_t noLpShift = outConstraint.noLP ? 1 : 0;
	Z_type iStackZ = Z_type(1);

	//////////  COMPUTE HYBRIDIZATION ENERGIES  ////////////

	// iterate over all window starts i1 (seq1) and i2 (seq2)
	for (i1=j1+1; i1-- > i1start; ) {
		// w1 = interaction width in seq1
		w1 = j1-i1+1;
		// screen for left boundaries in seq2
		for (i2=j2+1; i2-- > i2start; ) {

			// init: mark as invalid boundary
			hybridZ(i1,i2) = Z_type(0.0);

			// check if this cell is to be computed (!=E_INF)
			if( energy.isAccessible1(i1)
				&& energy.isAccessible2(i2)
				&& energy.areComplementary(i1,i2)
			)
			{
				// w2 = interaction width in seq2
				w2 = j2-i2+1;

				// reference access to cell value
				Z_type &curZ = hybridZ(i1,i2);

				// either interaction initiation
				if ( i1==j1 && i2==j2)  {
					if (noLpShift == 0) {
						// single base pair
						curZ = energy.getBoltzmannWeight(energy.getE_init());
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
						curZ = energy.getBoltzmannWeight(energy.getE_interLeft(i1,j1,i2,j2))
								* hybridZ(j1,j2);
					} else {
						// no lp allowed
						// check if right-side stacking of (i1,i2) is possible
						if (energy.isAccessible1(i1+noLpShift)
							&& energy.isAccessible2(i2+noLpShift)
							&& energy.areComplementary(i1+noLpShift,i2+noLpShift))
						{
							// get stacking term to avoid recomputation
							iStackZ = energy.getBoltzmannWeight(energy.getE_interLeft(i1,i1+noLpShift,i2,i2+noLpShift));

							// init with stacking only
							curZ = iStackZ * ((w1==2&&w2==2) ? energy.getBoltzmannWeight(energy.getE_init()) : hybridZ(i1+noLpShift, i2+noLpShift) );
						} else {
							//
							iStackZ = Z_INF;
						}
					}
					// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
					// ensure stacking is possible if no LP allowed
					if (w1 > 2 && w2 > 2 && Z_isNotINF(iStackZ)) {
						for (k1=std::min(j1-1,i1+energy.getMaxInternalLoopSize1()+1+noLpShift); k1>i1+noLpShift; k1--) {
						for (k2=std::min(j2-1,i2+energy.getMaxInternalLoopSize2()+1+noLpShift); k2>i2+noLpShift; k2--) {
							// ensure at least one unpaired base in connecting loop for noLP predictions
							if (outConstraint.noLP && k1-1==i1+noLpShift && k2-1==i2+noLpShift) {
								continue;
							}
							// check if (k1,k2) are valid left boundary
							if ( ! Z_equal( hybridZ(k1,k2), 0.0 ) ) {
								// update minimal value
								curZ += iStackZ * energy.getBoltzmannWeight(energy.getE_interLeft(i1+noLpShift,k1,i2+noLpShift,k2)) * hybridZ(k1,k2);
							}
						}
						}
					}
				}

				// update mfe if needed
				if (callUpdateZ) {
					updateZ(i1, j1, i2, j2, curZ, true);
				}

			} // complementary base pair
		}
	}

}

////////////////////////////////////////////////////////////////////////////


} // namespace
