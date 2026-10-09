
#include "IntaRNA/PredictorMfeEns2d.h"
#include "IntaRNA/PartitionArithmetic.h"

#include <stdexcept>

namespace IntaRNA {

////////////////////////////////////////////////////////////////////////////

PredictorMfeEns2d::
PredictorMfeEns2d(
		const InteractionEnergy & energy
		, OutputHandler & output
		, PredictionTracker * predTracker
		, BasePairProbabilities * pairProbabilities )
 : PredictorMfeEns(energy,output,predTracker)
	, pairProbabilities(pairProbabilities)
	, hybridZ( 0,0 )
{
	if (pairProbabilities && !output.getOutputConstraint().needZall)
		throw std::invalid_argument("unseeded pair probabilities require needZall");
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
	try { predictRegion(r1,r2); }
	catch (...) { if (pairProbabilities) pairProbabilities->fail(); throw; }
}

void
PredictorMfeEns2d::predictRegion(const IndexRange & r1,const IndexRange & r2)
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

	if (pairProbabilities && (!r1.isAscending() || !r2.isAscending()
			|| r1.from>=energy.getAccessibility1().getSequence().size()
			|| r2.from>=energy.getAccessibility2().getSequence().size()))
		throw std::invalid_argument("unseeded pair probabilities: invalid region");
	if (pairProbabilities && pairProbabilities->status()!=BasePairProbabilities::Status::pending)
		throw std::logic_error("unseeded pair probabilities: accumulator is not pending");

	// set index offset
	energy.setOffset1(r1.from);
	energy.setOffset2(r2.from);

	// resize matrix
	hybridZ.resize( std::min( energy.size1()
						, r1.to==RnaSequence::lastPos?energy.size1():r1.to-r1.from+1 )
				, std::min( energy.size2()
						, r2.to==RnaSequence::lastPos?energy.size2():r2.to-r2.from+1 ) );

	// initialize mfe interaction for updates
	initOptima();
	// initialize overall partition function for updates
	initZ();
	// No probability buffers or reverse work are needed for ordinary prediction.
	Z2dMatrix outside, masses;
	Z_type denominator=0;
	if (pairProbabilities) {
		outside.resize(hybridZ.size1(),hybridZ.size2());
		masses.resize(hybridZ.size1(),hybridZ.size2());
	}

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
			if (pairProbabilities)
				accumulatePairMasses(j1,j2,outside,masses,denominator);

		}
	}
	if (pairProbabilities && denominator!=PartitionArithmetic::check(Zall))
		throw std::logic_error("unseeded pair probabilities: updateZ override changed the partition objective");

	// report mfe interaction
	reportOptima();
	if (pairProbabilities) {
		const size_t n=hybridZ.size1(), m=hybridZ.size2();
		const auto first=energy.getBasePair(0,0), last=energy.getBasePair(n-1,m-1);
		// An unseeded prediction has no seed annotations, even for an API sink
		// that requested a mask for use by a shared writer.
		const bool annotate=pairProbabilities->collectsSeedPairs();
		Matrix<unsigned char> seeds(annotate?n:0,annotate?m:0,0);
		pairProbabilities->addRegion(IndexRange(first.first,last.first),
				IndexRange(last.second,first.second),denominator,masses,annotate?&seeds:nullptr);
	}
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
			if( energy.areComplementary(i1,i2)
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
						if (energy.areComplementary(i1+noLpShift,i2+noLpShift))
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
					updateCompleteZ(i1, j1, i2, j2, curZ, true);
				}

			} // complementary base pair
		}
	}

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2d::accumulatePairMasses(size_t j1,size_t j2,Z2dMatrix & outside,
		Z2dMatrix & masses,Z_type & denominator) const
{
	using A=PartitionArithmetic;
	const bool noLP=output.getOutputConstraint().noLP;
	const size_t shift=noLP?1:0;
	// Use exactly the rectangle and epsilon gates of fillHybridZ/updateZ: this
	// reverses the represented ensemble, including its existing zero cutoff.
	const size_t begin1=j1-std::min(j1,energy.getAccessibility1().getMaxLength()+1);
	const size_t begin2=j2-std::min(j2,energy.getAccessibility2().getMaxLength()+1);
	for(size_t i1=j1+1;i1-- > begin1;) for(size_t i2=j2+1;i2-- > begin2;) {
		outside(i1,i2)=0;
		const Z_type h=A::check(hybridZ(i1,i2));
		if (!Z_equal(h,0) && isValidOutputSite(i1,j1,i2,j2)) {
			const Z_type b=A::check(energy.getBoltzmannWeight(energy.getE(i1,j1,i2,j2,E_type(0))));
			outside(i1,i2)=b;
			denominator=A::add(denominator,A::multiply(h,b));
		}
	}

	const auto first=energy.getBasePair(0,0);
	const auto last=energy.getBasePair(masses.size1()-1,masses.size2()-1);
	auto addMass=[&](size_t i1,size_t i2,Z_type value) {
		const auto bp=energy.getBasePair(i1,i2);
		auto & mass=masses(bp.first-first.first,bp.second-last.second);
		mass=A::add(mass,value);
	};
	const Z_type initiation=A::check(energy.getBoltzmannWeight(energy.getE_init()));
	// Every suffix state contains its left pair once. Reverse topological
	// traversal gathers all enclosing left boundaries before visiting a suffix.
	for(size_t i1=begin1;i1<=j1;++i1) for(size_t i2=begin2;i2<=j2;++i2) {
		const Z_type a=outside(i1,i2), h=hybridZ(i1,i2);
		if (a==0 || h==0) continue;
		addMass(i1,i2,A::multiply(a,h));
		const size_t w1=j1-i1+1, w2=j2-i2+1;
		if (w1<=1 || w2<=1) continue;

		Z_type stack=1;
		if (!noLP) {
			const Z_type edge=A::multiply(a,
					energy.getBoltzmannWeight(energy.getE_interLeft(i1,j1,i2,j2)));
			outside(j1,j2)=A::add(outside(j1,j2),edge);
		} else {
			if (!energy.areComplementary(i1+1,i2+1)) continue;
			stack=A::check(energy.getBoltzmannWeight(energy.getE_interLeft(i1,i1+1,i2,i2+1)));
			const Z_type edge=A::multiply(a,stack);
			if (w1==2 && w2==2) {
				// The noLP terminal has no singleton H state: its second pair
				// is implicit in stack * initiation and must be counted here.
				addMass(j1,j2,A::multiply(edge,initiation));
			} else {
				outside(i1+1,i2+1)=A::add(outside(i1+1,i2+1),edge);
			}
		}
		if (w1<=2 || w2<=2) continue;
		for(size_t k1=std::min(j1-1,i1+energy.getMaxInternalLoopSize1()+1+shift);k1>i1+shift;--k1) {
			for(size_t k2=std::min(j2-1,i2+energy.getMaxInternalLoopSize2()+1+shift);k2>i2+shift;--k2) {
				if (noLP && k1-1==i1+shift && k2-1==i2+shift) continue;
				if (Z_equal(hybridZ(k1,k2),0)) continue;
				const Z_type loop=A::check(energy.getBoltzmannWeight(
						energy.getE_interLeft(i1+shift,k1,i2+shift,k2)));
				const Z_type edge=A::multiply(a,A::multiply(stack,loop));
				outside(k1,k2)=A::add(outside(k1,k2),edge);
				if (noLP) {
					// This decomposition consumes (i,i+) before jumping over
					// the loop to k; only i and k have explicit suffix states.
					addMass(i1+1,i2+1,A::multiply(edge,hybridZ(k1,k2)));
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////


} // namespace
