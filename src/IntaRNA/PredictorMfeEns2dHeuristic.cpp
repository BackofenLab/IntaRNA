
#include "IntaRNA/PredictorMfeEns2dHeuristic.h"
#include "IntaRNA/PartitionArithmetic.h"

#include <stdexcept>

namespace IntaRNA {

////////////////////////////////////////////////////////////////////////////

PredictorMfeEns2dHeuristic::
PredictorMfeEns2dHeuristic(
		const InteractionEnergy & energy
		, OutputHandler & output
		, PredictionTracker * predTracker
		, BasePairProbabilities * pairProbabilities )
 : PredictorMfeEns2d(energy,output,predTracker,pairProbabilities)
{
}


////////////////////////////////////////////////////////////////////////////

PredictorMfeEns2dHeuristic::
~PredictorMfeEns2dHeuristic()
{
	// clean up
}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2dHeuristic::
predict( const IndexRange & r1
		, const IndexRange & r2
		)
{
	try { predictRegionHeuristic(r1,r2); }
	catch (...) { if (pairProbabilities) pairProbabilities->fail(); throw; }
}

void
PredictorMfeEns2dHeuristic::
predictRegionHeuristic(const IndexRange & r1,const IndexRange & r2)
{
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) <<"predicting ensemble mfe interactions heuristically in O(n^2) space and time..."; }
	// measure timing
	TIMED_FUNC_IF(timerObj,VLOG_IS_ON(9));

#if INTARNA_IN_DEBUG_MODE
	// check indices
	if (!(r1.isAscending() && r2.isAscending()) )
		throw std::runtime_error("PredictorMfeEns2dHeuristic::predict("+toString(r1)+","+toString(r2)+") is not sane");
#endif


	if (!r1.isAscending() || !r2.isAscending()
			|| r1.from>=energy.getAccessibility1().getSequence().size()
			|| r2.from>=energy.getAccessibility2().getSequence().size())
		throw std::invalid_argument("heuristic ensemble prediction: invalid region");
	if (pairProbabilities && pairProbabilities->status()!=BasePairProbabilities::Status::pending)
		throw std::logic_error("heuristic pair probabilities: accumulator is not pending");
	if (pairProbabilities) pairProbabilities->markApproximate();
	// set index offset
	energy.setOffset1(r1.from);
	energy.setOffset2(r2.from);

	// resize matrix
	hybridZ.resize( std::min( energy.size1()
						, r1.to==RnaSequence::lastPos?energy.size1():r1.to-r1.from+1 )
				, std::min( energy.size2()
						, r2.to==RnaSequence::lastPos?energy.size2():r2.to-r2.from+1 ) );
	if (pairProbabilities) {
		probabilityContinuation = Matrix<ProbabilityContinuation>(hybridZ.size1(),hybridZ.size2());
		probabilityFlow = Matrix<Z_type>(hybridZ.size1(),hybridZ.size2(),0);
		probabilityMass = Matrix<Z_type>(hybridZ.size1(),hybridZ.size2(),0);
		probabilityDenominator = 0;
	}

	// init mfe for later updates
	initOptima();
	// initialize overall partition function for updates
	initZ();

	// compute table and update mfeInteraction
	fillHybridZ();

	// trace back and output handler update
	reportOptima();
	if (pairProbabilities) commitProbabilityRegion();

}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2dHeuristic::
fillHybridZ()
{
	// temporary access
	const OutputConstraint & outConstraint = output.getOutputConstraint();
	// compute entries
	// current minimal value
	Z_type curZ = Z_INF;
	E_type curEtotal = E_INF, curCellEtotal = E_INF;
	size_t i1,i2,w1,w2;

	// determine whether or not lonely base pairs are allowed or if we have to
	// ensure a stacking to the right of the left boundary (i1,i2)
	const size_t noLpShift = outConstraint.noLP ? 1 : 0;
	Z_type iStackZ = Z_type(1);

	BestInteractionZ * curCell = NULL;
	const BestInteractionZ * rightExt = NULL;
	// iterate (decreasingly) over all left interaction starts
	for (i1=hybridZ.size1(); i1-- > 0;) {
		for (i2=hybridZ.size2(); i2-- > 0;) {
			// direct cell access
			curCell = &(hybridZ(i1,i2));

			// init as invalid boundary
			*curCell = BestInteractionZ(0.0, RnaSequence::lastPos, RnaSequence::lastPos);
			curCellEtotal = E_INF;

			// check if positions can form interaction
			if (energy.areComplementary(i1,i2) )
			{
				// no lp allowed
				if (noLpShift != 0) {
					// check if right-side stacking of (i1,i2) is possible
					if ( i1+noLpShift < hybridZ.size1()
						&& i2+noLpShift < hybridZ.size2()
						&& energy.areComplementary(i1+noLpShift,i2+noLpShift))
					{
						// get stacking term to avoid recomputation
						iStackZ = energy.getBoltzmannWeight(energy.getE_interLeft(i1,i1+noLpShift,i2,i2+noLpShift));
					} else {
						// skip further processing, since no stacking possible
						continue;
					}
				}

				// if valid right boundary
				if (!outConstraint.noGUend || !energy.isGU(i1+noLpShift,i2+noLpShift))
				{
					// set to interaction initiation with according boundary
					*curCell = BestInteractionZ(iStackZ * energy.getBoltzmannWeight(energy.getE_init()), i1+noLpShift, i2+noLpShift);
					// current best total energy value (covers to far E_init only)
					curCellEtotal = energy.getE(i1,i1+noLpShift, i2,i2+noLpShift ,energy.getE(curCell->val));
					// update overall partition function information for initial bps only
					updateZ( i1,curCell->j1, i2,curCell->j2, curCell->val, true );
					if (pairProbabilities) {
						probabilityContinuation(i1,i2)={RnaSequence::lastPos,RnaSequence::lastPos,noLpShift!=0};
						recordProbabilityCandidate(i1,curCell->j1,i2,curCell->j2,curCell->val,probabilityContinuation(i1,i2));
					}

				}

				if(outConstraint.noLP) {
					/////////////////////////////////////////
					// check direct extension to the right of the noLP stacking
					/////////////////////////////////////////

					// direct cell access (const)
					rightExt = &(hybridZ(i1+noLpShift,i2+noLpShift));
					// check if right side can pair and interaction length is within boundary
					if (!Z_equal(rightExt->val, 0.0)
						&& (rightExt->j1 +1 -i1) <= energy.getAccessibility1().getMaxLength()
						&& (rightExt->j2 +1 -i2) <= energy.getAccessibility2().getMaxLength() )
					{
						// compute Z for direct extension with stacking
						curZ = iStackZ * rightExt->val;

						// update overall partition function information for current right extension
						updateZ( i1,rightExt->j1, i2,rightExt->j2, curZ, true );
						if (pairProbabilities)
							recordProbabilityCandidate(i1,rightExt->j1,i2,rightExt->j2,curZ,{i1+noLpShift,i2+noLpShift,false});

						// check if this combination yields better energy
						curEtotal = energy.getE(i1,rightExt->j1, i2,rightExt->j2, energy.getE(curZ));

						// update best right extension for (i1,i2) in curCell
						if ( curEtotal < curCellEtotal )
						{
							// update current best for this left boundary
							// copy right boundary
							*curCell = *rightExt;
							// set new partition function
							curCell->val = curZ;
							// store total energy to avoid recomputation
							curCellEtotal = curEtotal;
							if (pairProbabilities)
								probabilityContinuation(i1,i2)={i1+noLpShift,i2+noLpShift,false};
						}
					}
				}


				// iterate over all loop sizes w1 (seq1) and w2 (seq2) (minus 1)
				for (w1=1; w1-1 <= energy.getMaxInternalLoopSize1() && i1+w1+noLpShift<hybridZ.size1(); w1++) {
				for (w2=1; w2-1 <= energy.getMaxInternalLoopSize2() && i2+w2+noLpShift<hybridZ.size2(); w2++) {
					// For noLP, the adjacent continuation is already represented
					// by the direct extension above. Counting it again as the
					// (w1,w2)=(1,1) loop duplicates the same interaction paths.
					if (noLpShift != 0 && w1 == 1 && w2 == 1) {
						continue;
					}
					// direct cell access (const)
					rightExt = &(hybridZ(i1+noLpShift+w1,i2+noLpShift+w2));
					// check if right side can pair
					if (Z_equal(rightExt->val, 0.0)) {
						continue;
					}
					// check if interaction length is within boundary
					if ( (rightExt->j1 +1 -i1) > energy.getAccessibility1().getMaxLength()
						|| (rightExt->j2 +1 -i2) > energy.getAccessibility2().getMaxLength() )
					{
						continue;
					}

					// compute Z for this loop sizes
					curZ = iStackZ * energy.getBoltzmannWeight(energy.getE_interLeft(i1+noLpShift,i1+noLpShift+w1,i2+noLpShift,i2+noLpShift+w2)) * rightExt->val;

					// update overall partition function information for current right extension
					updateZ( i1,rightExt->j1, i2,rightExt->j2, curZ, true );
					if (pairProbabilities)
						recordProbabilityCandidate(i1,rightExt->j1,i2,rightExt->j2,curZ,{i1+noLpShift+w1,i2+noLpShift+w2,noLpShift!=0});

					// check if this combination yields better energy
					curEtotal = energy.getE(i1,rightExt->j1, i2,rightExt->j2, energy.getE(curZ));

					// update best right extension for (i1,i2) in curCell
					if ( curEtotal < curCellEtotal )
					{
						// update current best for this left boundary
						// copy right boundary
						*curCell = *rightExt;
						// set new partition function
						curCell->val = curZ;
						// store total energy to avoid recomputation
						curCellEtotal = curEtotal;
						if (pairProbabilities)
							probabilityContinuation(i1,i2)={i1+noLpShift+w1,i2+noLpShift+w2,noLpShift!=0};
					}

				} // w2
				} // w1

			} // valid base pair

		} // i2
	} // i1

}

////////////////////////////////////////////////////////////////////////////

void PredictorMfeEns2dHeuristic::recordProbabilityCandidate(size_t i1,size_t j1,
		size_t i2,size_t j2,Z_type hybrid,const ProbabilityContinuation & continuation)
{
	// Match the legacy heuristic's numerical-zero and output-site admission.
	if (Z_equal(hybrid,0) || !isValidOutputSite(i1,j1,i2,j2)) return;
	const E_type boundaryEnergy=energy.getE(i1,j1,i2,j2,E_type(0));
	const Z_type boundary=energy.getBoltzmannWeight(boundaryEnergy);
	if (boundary==0 && !E_isINF(boundaryEnergy))
		throw std::range_error("heuristic pair probabilities: boundary weight underflow");
	const Z_type weight=PartitionArithmetic::multiply(hybrid,boundary);
	probabilityDenominator=PartitionArithmetic::add(probabilityDenominator,weight);
	probabilityMass(i1,i2)=PartitionArithmetic::add(probabilityMass(i1,i2),weight);
	if (continuation.extraPair)
		probabilityMass(i1+1,i2+1)=PartitionArithmetic::add(probabilityMass(i1+1,i2+1),weight);
	if (continuation.next1!=RnaSequence::lastPos)
		probabilityFlow(continuation.next1,continuation.next2)=PartitionArithmetic::add(
				probabilityFlow(continuation.next1,continuation.next2),weight);
}

void PredictorMfeEns2dHeuristic::commitProbabilityRegion()
{
	if (probabilityDenominator!=Zall)
		throw std::logic_error("heuristic updateZ override did not preserve the pair-probability objective");
	const size_t n=hybridZ.size1(),m=hybridZ.size2();
	// Every retained cell represents exactly one chain. Its incoming candidate
	// mass can therefore pass unchanged to the selected child, without division.
	for(size_t i=0;i<n;++i) for(size_t j=0;j<m;++j) {
		const Z_type weight=probabilityFlow(i,j);
		if (weight==0) continue;
		const auto & continuation=probabilityContinuation(i,j);
		probabilityMass(i,j)=PartitionArithmetic::add(probabilityMass(i,j),weight);
		if (continuation.extraPair)
			probabilityMass(i+1,j+1)=PartitionArithmetic::add(probabilityMass(i+1,j+1),weight);
		if (continuation.next1!=RnaSequence::lastPos)
			probabilityFlow(continuation.next1,continuation.next2)=PartitionArithmetic::add(
					probabilityFlow(continuation.next1,continuation.next2),weight);
	}
	const auto first=energy.getBasePair(0,0),last=energy.getBasePair(n-1,m-1);
	const IndexRange target(first.first,last.first),query(last.second,first.second);
	Matrix<Z_type> original(n,m,0);
	for(size_t i=0;i<n;++i) for(size_t j=0;j<m;++j) {
		const auto bp=energy.getBasePair(i,j);
		original(bp.first-target.from,bp.second-query.from)=probabilityMass(i,j);
	}
	const bool annotate=pairProbabilities->collectsSeedPairs();
	Matrix<unsigned char> seeds(annotate?n:0,annotate?m:0,0);
	pairProbabilities->addRegion(target,query,Zall,original,annotate?&seeds:nullptr);
}

////////////////////////////////////////////////////////////////////////////

} // namespace
