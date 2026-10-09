
#include "IntaRNA/PredictorMfeEns2dSeedExtension.h"
#include "IntaRNA/SeededPartitionFunction.h"
#include "IntaRNA/PartitionArithmetic.h"
#include <cmath>
#include <limits>

namespace IntaRNA {

// The heuristic objective contains signed seed-overlap corrections. Trace the
// implemented arithmetic directly, retaining pair ownership on every term.
// This is allocated only for requested probability output; ordinary prediction
// follows its original arithmetic and storage path.
struct PredictorMfeEns2dSeedExtension::ExtensionProbabilityTrace {
	using Pair=Interaction::BasePair;
	struct Edge { size_t child; Z_type weight; std::vector<Pair> pairs; };
	struct Tape {
		size_t columns=0;
		std::vector<std::vector<Edge>> edges;
		std::vector<Z_type> adjoint;
		void reset(size_t n,size_t m);
		void add(size_t p1,size_t p2,size_t c1,size_t c2,Z_type weight,std::vector<Pair> pairs);
		void reverse(const Matrix<Z_type> & values,ExtensionProbabilityTrace & trace);
	};
	Matrix<Z_type> mass,absoluteMass;
	Matrix<unsigned char> seeds;
	Tape left,right;
	std::vector<Pair> anchor;
	size_t si1=0,si2=0,sj1=0,sj2=0;
	Z_type seedZ=0,z=0;
	explicit ExtensionProbabilityTrace(size_t n,size_t m);
	static Z_type finite(Z_type x);
	static Z_type product(Z_type a,Z_type b);
	void addMass(const Pair & bp,Z_type value);
};

void
PredictorMfeEns2dSeedExtension::ExtensionProbabilityTrace::Tape::reset(size_t n,size_t m)
{
	columns=m;
	edges.assign(n*m,{});
	adjoint.assign(n*m,0);
}

void
PredictorMfeEns2dSeedExtension::ExtensionProbabilityTrace::Tape::add(
		size_t p1,size_t p2,size_t c1,size_t c2,Z_type weight,std::vector<Pair> pairs)
{
	finite(weight);
	edges.at(p1*columns+p2).push_back({c1*columns+c2,weight,std::move(pairs)});
}

void
PredictorMfeEns2dSeedExtension::ExtensionProbabilityTrace::Tape::reverse(
		const Matrix<Z_type> & values,ExtensionProbabilityTrace & trace)
{
	for(size_t parent=edges.size();parent-- >0;) {
		finite(values(parent/columns,parent%columns));
		for(const auto & edge:edges[parent]) {
			const Z_type a=product(adjoint[parent],edge.weight);
			adjoint[edge.child]=finite(adjoint[edge.child]+a);
			const Z_type mass=product(a,values(edge.child/columns,edge.child%columns));
			for(const auto & bp:edge.pairs) trace.addMass(bp,mass);
		}
	}
}

PredictorMfeEns2dSeedExtension::ExtensionProbabilityTrace::ExtensionProbabilityTrace(size_t n,size_t m)
	: mass(n,m,0),absoluteMass(n,m,0),seeds(n,m,0)
{}

Z_type
PredictorMfeEns2dSeedExtension::ExtensionProbabilityTrace::finite(Z_type x)
{
	if (!std::isfinite(x)) throw std::range_error("heuristic seeded probabilities: nonfinite arithmetic");
	return x;
}

Z_type
PredictorMfeEns2dSeedExtension::ExtensionProbabilityTrace::product(Z_type a,Z_type b)
{
	finite(a); finite(b);
	const Z_type value=finite(a*b);
	if(a!=0 && b!=0 && (value==0 || std::abs(value)<std::numeric_limits<Z_type>::min()))
		throw std::range_error("heuristic seeded probabilities: product underflow");
	return value;
}

void
PredictorMfeEns2dSeedExtension::ExtensionProbabilityTrace::addMass(const Pair & bp,Z_type value)
{
	mass(bp.first,bp.second)=finite(mass(bp.first,bp.second)+value);
	absoluteMass(bp.first,bp.second)=finite(absoluteMass(bp.first,bp.second)+std::abs(value));
}

void PredictorMfeEns2dSeedExtension::initExtensionProbabilities(size_t n,size_t m)
{
	if(pairProbabilities) extensionProbabilityTrace=std::make_unique<ExtensionProbabilityTrace>(n,m);
}

std::vector<Interaction::BasePair>
PredictorMfeEns2dSeedExtension::extensionSeedPrefix(size_t i1,size_t i2,size_t before) const
{
	Interaction seed(energy.getAccessibility1().getSequence(),energy.getAccessibility2().getAccessibilityOrigin().getSequence());
	seed.basePairs.push_back(energy.getBasePair(i1,i2));
	seedHandler.traceBackSeed(seed,i1,i2);
	seed.basePairs.push_back(energy.getBasePair(i1+seedHandler.getSeedLength1(i1,i2)-1,
			i2+seedHandler.getSeedLength2(i1,i2)-1));
	std::vector<Interaction::BasePair> pairs;
	for(const auto & bp:seed.basePairs) {
		const size_t p=energy.getIndex1(bp),q=energy.getIndex2(bp);
		if(p<before) pairs.emplace_back(p,q);
	}
	return pairs;
}

void PredictorMfeEns2dSeedExtension::beginExtensionProbabilitySeed(size_t i1,size_t i2)
{
	if(!extensionProbabilityTrace) return;
	auto & trace=*extensionProbabilityTrace;
	trace.si1=i1; trace.si2=i2;
	trace.sj1=i1+seedHandler.getSeedLength1(i1,i2)-1;
	trace.sj2=i2+seedHandler.getSeedLength2(i1,i2)-1;
	trace.seedZ=PartitionArithmetic::check(energy.getBoltzmannWeight(seedHandler.getSeedE(i1,i2)));
	trace.anchor=extensionSeedPrefix(i1,i2,trace.sj1+1);
	for(const auto & bp:trace.anchor) trace.seeds(bp.first,bp.second)=1;
	trace.left.reset(0,0); trace.right.reset(0,0);
}

void PredictorMfeEns2dSeedExtension::addExtensionProbabilityRoot(size_t i1,size_t j1,size_t i2,size_t j2,
		Z_type partZ,bool left,bool right)
{
	if(!extensionProbabilityTrace || Z_equal(partZ,0) || !isValidOutputSite(i1,j1,i2,j2)) return;
	auto & trace=*extensionProbabilityTrace;
	using A=PartitionArithmetic;
	const E_type boundaryE=energy.getE(i1,j1,i2,j2,E_type(0));
	const Z_type boundary=A::check(energy.getBoltzmannWeight(boundaryE));
	if(boundary==0 && !E_isINF(boundaryE)) throw std::range_error("heuristic seeded probabilities: boundary underflow");
	const Z_type weight=A::multiply(partZ,boundary);
	trace.z=A::add(trace.z,weight);
	for(const auto & bp:trace.anchor) trace.addMass(bp,weight);
	const Z_type l=left?hybridZ_left(trace.si1-i1,trace.si2-i2):energy.getBoltzmannWeight(energy.getE_init());
	const Z_type r=right?hybridZ_right(j1-trace.sj1,j2-trace.sj2):Z_type(1);
	if(left) {
		Z_type & a=trace.left.adjoint.at((trace.si1-i1)*trace.left.columns+trace.si2-i2);
		a=A::add(a,A::multiply(A::multiply(boundary,trace.seedZ),r));
	}
	if(right) {
		Z_type & a=trace.right.adjoint.at((j1-trace.sj1)*trace.right.columns+j2-trace.sj2);
		a=A::add(a,A::multiply(A::multiply(boundary,trace.seedZ),l));
	}
}

void PredictorMfeEns2dSeedExtension::finishExtensionProbabilitySeed()
{
	if(!extensionProbabilityTrace) return;
	auto & trace=*extensionProbabilityTrace;
	trace.left.reverse(hybridZ_left,trace);
	trace.right.reverse(hybridZ_right,trace);
}

void PredictorMfeEns2dSeedExtension::commitExtensionProbabilities()
{
	if(!extensionProbabilityTrace) return;
	auto & trace=*extensionProbabilityTrace;
	PartitionArithmetic::check(Zall);
	const size_t n=trace.mass.size1(),m=trace.mass.size2();
	const Z_type tolerance=512*std::numeric_limits<Z_type>::epsilon()*Z_type(n+m+1);
	if(std::abs(Zall-trace.z)>tolerance*std::max(std::abs(trace.z),std::abs(Zall)))
		throw std::logic_error("heuristic seeded probabilities: outside objective differs from reported partition");
	const auto first=energy.getBasePair(0,0),last=energy.getBasePair(n-1,m-1);
	const IndexRange target(first.first,last.first),query(last.second,first.second);
	Matrix<Z_type> original(n,m,0);
	Matrix<unsigned char> annotations(n,m,0);
	for(size_t i=0;i<n;++i) for(size_t j=0;j<m;++j) {
		Z_type mass=trace.mass(i,j);
		if(mass<0 && -mass<=tolerance*trace.absoluteMass(i,j)) mass=0;
		PartitionArithmetic::check(mass);
		const auto bp=energy.getBasePair(i,j);
		original(bp.first-target.from,bp.second-query.from)=mass;
		annotations(bp.first-target.from,bp.second-query.from)=trace.seeds(i,j);
	}
	pairProbabilities->addRegion(target,query,trace.z,original,&annotations);
	extensionProbabilityTrace.reset();
}

//////////////////////////////////////////////////////////////////////////

PredictorMfeEns2dSeedExtension::
PredictorMfeEns2dSeedExtension(
		const InteractionEnergy & energy
		, OutputHandler & output
		, PredictionTracker * predTracker
		, SeedHandler * seedHandlerInstance
		, BasePairProbabilities * pairProbabilities )
 :
	PredictorMfeEns(energy,output,predTracker)
	, pairProbabilities(pairProbabilities)
	, seedHandler(seedHandlerInstance)
	, hybridZ_left( 0,0 )
	, hybridZ_right( 0,0 )
{
	if (pairProbabilities && !output.getOutputConstraint().needZall)
		throw std::invalid_argument("pair probabilities require needZall");
}

//////////////////////////////////////////////////////////////////////////

PredictorMfeEns2dSeedExtension::
~PredictorMfeEns2dSeedExtension()
{
}

//////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2dSeedExtension::
predict( const IndexRange & r1, const IndexRange & r2 )
{
	if (seedHandler.guaranteesStackOnlySeeds()) {
		try { predictStackSeeds(r1,r2); }
		catch (...) { if (pairProbabilities) pairProbabilities->fail(); throw; }
		return;
	}
	if (pairProbabilities) {
		pairProbabilities->fail();
		throw std::invalid_argument("exact bulged-seed pair probabilities are unsupported");
	}
	if (seedHandler.getConstraint().getBasePairs()<2)
		throw std::invalid_argument("legacy bulged seed extension requires at least two seed pairs");
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) <<"predicting ensemble mfe interactions with seed in O(n^2) space and O(n^4) time..."; }
	// measure timing
	TIMED_FUNC_IF(timerObj,VLOG_IS_ON(9));

#if INTARNA_IN_DEBUG_MODE
	// check indices
	if (!(r1.isAscending() && r2.isAscending()) )
		throw std::runtime_error("PredictorMfeEns2dSeedExtension::predict("+toString(r1)+","+toString(r2)+") is not sane");
#endif

	// setup index offset
	energy.setOffset1(r1.from);
	energy.setOffset2(r2.from);
	seedHandler.setOffset1(r1.from);
	seedHandler.setOffset2(r2.from);

	const size_t range_size1 = std::min( energy.size1()
			, (r1.to==RnaSequence::lastPos?energy.size1()-1:r1.to)-r1.from+1 );
	const size_t range_size2 = std::min( energy.size2()
			, (r2.to==RnaSequence::lastPos?energy.size2()-1:r2.to)-r2.from+1 );

	// compute seed interactions for whole range
	// and check if any seed possible
	if (seedHandler.fillSeed( 0, range_size1-1, 0, range_size2-1 ) == 0) {
		// trigger empty interaction reporting
		initOptima();
		initZ();
		reportOptima();
		// stop computation
		return;
	}

	// initialize mfe interaction for updates
	initOptima();
	// initialize overall partition function for updates
	initZ();

	size_t si1 = RnaSequence::lastPos, si2 = RnaSequence::lastPos;
	while( seedHandler.updateToNextSeed(si1,si2
			, 0, range_size1+1-seedHandler.getConstraint().getBasePairs()
			, 0, range_size2+1-seedHandler.getConstraint().getBasePairs()) )
	{
		// get Z and boundaries of seed
		const Z_type seedZ = energy.getBoltzmannWeight( seedHandler.getSeedE(si1, si2) );
		const size_t sl1 = seedHandler.getSeedLength1(si1, si2);
		const size_t sl2 = seedHandler.getSeedLength2(si1, si2);
		const size_t sj1 = si1+sl1-1;
		const size_t sj2 = si2+sl2-1;
		// check if seed fits into interaction range
		if (sj1 > range_size1 || sj2 > range_size2)
			continue;
		const size_t maxMatrixLen1 = energy.getAccessibility1().getMaxLength()-sl1+1;
		const size_t maxMatrixLen2 = energy.getAccessibility2().getMaxLength()-sl2+1;

		// ER
		hybridZ_right.resize( std::min(range_size1-sj1, maxMatrixLen1), std::min(range_size2-sj2, maxMatrixLen2) );
		fillHybridZ_right(sj1, sj2);

		// EL
		hybridZ_left.resize( std::min(si1+1, maxMatrixLen1), std::min(si2+1, maxMatrixLen2) );
		fillHybridZ_left(si1, si2);

		// updateZ for all boundary combinations
		for (size_t l1 = 0; l1<hybridZ_left.size1(); l1++) {
			for (size_t l2 = 0; l2< hybridZ_left.size2(); l2++) {
				// check complementarity of boundary
				if ( Z_equal(hybridZ_left(l1,l2), 0.0) ) continue;
				// iterate extension right of seed in seq 1
				for (size_t r1 = 0; r1 < hybridZ_right.size1() ; r1++) {
					// ensure max interaction length in seq 1
					if (sj1+r1-si1+l1 >= energy.getAccessibility1().getMaxLength()) break;
					// iterate extension right of seed in seq 2
					for (size_t r2 = 0; r2 < hybridZ_right.size2() ; r2++) {
						// ensure max interaction length in seq 2
						if (sj2+r2-si2+l2 >= energy.getAccessibility2().getMaxLength()) break;
						// check complementarity of boundary
						if (Z_equal(hybridZ_right(r1,r2),0.0)) continue;
						// compute partition function given the current seed
						updateZ(si1-l1, sj1+r1, si2-l2, sj2+r2, hybridZ_left(l1,l2) * seedZ * hybridZ_right(r1,r2), true);

					} // r2
				} // l2
			} // r1
		} // l1

	} // si1 / si2

	// report mfe interaction
	reportOptima();

}

//////////////////////////////////////////////////////////////////////////

Z_type
PredictorMfeEns2dSeedExtension::exactBoundaryWeight(size_t i1,size_t j1,size_t i2,size_t j2) const
{
	if (!isValidOutputSite(i1,j1,i2,j2)) return 0;
	const E_type e=energy.getE(i1,j1,i2,j2,E_type(0));
	return E_isINF(e)?Z_type(0):PartitionArithmetic::exp(-E_2_Z(e)/energy.getRT());
}

void
PredictorMfeEns2dSeedExtension::predictStackSeeds(const IndexRange & r1,const IndexRange & r2)
{
	const size_t size1=energy.getAccessibility1().getSequence().size();
	const size_t size2=energy.getAccessibility2().getSequence().size();
	if (!r1.isAscending() || !r2.isAscending() || r1.from>=size1 || r2.from>=size2)
		throw std::invalid_argument("exact seeded prediction: invalid region");
	energy.setOffset1(r1.from); energy.setOffset2(r2.from);
	seedHandler.setOffset1(r1.from); seedHandler.setOffset2(r2.from);
	const size_t n=r1.to==RnaSequence::lastPos?energy.size1():std::min(energy.size1(),r1.to-r1.from+1);
	const size_t m=r2.to==RnaSequence::lastPos?energy.size2():std::min(energy.size2(),r2.to-r2.from+1);
	initOptima(); initZ(); Zall=0;
	const size_t seedCount=seedHandler.fillSeed(0,n-1,0,m-1);
	const size_t span1=std::min(n,energy.getAccessibility1().getMaxLength());
	const size_t span2=std::min(m,energy.getAccessibility2().getMaxLength());
	StackSeedDomain seeds=seedCount ? StackSeedDomain(seedHandler,n,m,span1,span2)
			: StackSeedDomain(std::vector<StackSeedDomain::Occurrence>(),n,m,span1,span2);
	using K=SeededPartitionFunction;
	const K::Domain domain{n,m,span1,span2,energy.getMaxInternalLoopSize1(),energy.getMaxInternalLoopSize2(),output.getOutputConstraint().noLP};
	auto boltzmann=[&](E_type e) { return E_isINF(e)?Z_type(0):PartitionArithmetic::exp(-E_2_Z(e)/energy.getRT()); };
	K::Weights weights{
		[&](K::Pair p){return energy.isAccessible1(p[0]) && energy.isAccessible2(p[1]) && energy.areComplementary(p[0],p[1]);},
		[&](K::Pair){return boltzmann(energy.getE_init());},
		[&](K::Pair p,K::Pair q){return boltzmann(energy.getE_interLeft(p[0],q[0],p[1],q[1]));},
		[&](K::Pair p,K::Pair q){return exactBoundaryWeight(p[0],q[0],p[1],q[1]);},
		[&](K::Pair p,K::Pair q,Z_type h,Z_type b){updateExactCompleteZ(p[0],q[0],p[1],q[1],h,b);}
	};
	const auto result=K::compute(domain,seeds,weights,pairProbabilities!=nullptr);
	if (result.z!=Zall) throw std::logic_error("exact seeded updateZ override did not preserve the partition objective");
	output.setExactPartition(true);
	reportOptima();
	if (pairProbabilities) {
		const auto first=energy.getBasePair(0,0), last=energy.getBasePair(n-1,m-1);
		const IndexRange target(first.first,last.first), query(last.second,first.second);
		Matrix<Z_type> original(n,m,0);
		const bool annotate=pairProbabilities->collectsSeedPairs();
		Matrix<unsigned char> seedPairs(annotate?n:0,annotate?m:0,0);
		for(size_t i=0;i<n;++i) for(size_t j=0;j<m;++j) {
			const auto bp=energy.getBasePair(i,j);
			original(bp.first-target.from,bp.second-query.from)=result.mass(i,j);
			// Mark the same admitted, contained occurrences used by the kernel.
			if (annotate) for(size_t k=0;k<seeds.seedLength(i,j);++k) {
				const auto seedPair=energy.getBasePair(i+k,j+k);
				seedPairs(seedPair.first-target.from,seedPair.second-query.from)=1;
			}
		}
		pairProbabilities->addRegion(target,query,result.z,original,annotate?&seedPairs:nullptr);
	}
}

E_type
PredictorMfeEns2dSeedExtension::
getNonOverlappingEnergy( const size_t si1, const size_t si2, const size_t si1p, const size_t si2p ) {

#if INTARNA_IN_DEBUG_MODE
	// check indices
	if( !seedHandler.isSeedBound(si1,si2) )
		throw std::runtime_error("PredictorMfeEns2dSeedExtension::getNonOverlappingEnergy( si "+toString(si1)+","+toString(si2)+",..) is no seed bound");
	if( !seedHandler.isSeedBound(si1p,si2p) )
		throw std::runtime_error("PredictorMfeEns2dSeedExtension::getNonOverlappingEnergy( sip "+toString(si1p)+","+toString(si2p)+",..) is no seed bound");
	if( si1 > si1p )
		throw std::runtime_error("PredictorMfeEns2dSeedExtension::getNonOverlappingEnergy( si "+toString(si1)+","+toString(si2)+", sip "+toString(si1p)+","+toString(si2p)+",..) si1 > sj1 !");
	// check if loop-overlapping (i.e. share at least one loop)
	if( !seedHandler.areLoopOverlapping(si1,si2,si1p,si2p) ) {
		throw std::runtime_error("PredictorMfeEns2dSeedExtension::getNonOverlappingEnergy( si "+toString(si1)+","+toString(si2)+", sip "+toString(si1p)+","+toString(si2p)+",..) are not loop overlapping");
	}
#endif

	// identity check
	if( si1 == si1p ) {
		return E_type(0);
	}

	// trace seed at (si1,si2)
	Interaction interaction = Interaction(energy.getAccessibility1().getSequence(), energy.getAccessibility2().getAccessibilityOrigin().getSequence());
	interaction.basePairs.push_back( energy.getBasePair(si1, si2) );
	seedHandler.traceBackSeed( interaction, si1, si2 );

	E_type fullE = 0;
	size_t k1old = si1, k2old = si2;
	for (size_t i = 1; i < interaction.basePairs.size(); i++) {
		// get index of current base pair
		size_t k1 = energy.getIndex1(interaction.basePairs[i]);
		// check if overlap done
		if (k1 > si1p) break;
		size_t k2 = energy.getIndex2(interaction.basePairs[i]);
		// add hybridization energy
		fullE += energy.getE_interLeft(k1old,k1,k2old,k2);

		// store
		k1old = k1;
		k2old = k2;
	}
	return fullE;
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2dSeedExtension::
fillHybridZ_left( const size_t si1, const size_t si2 )
{
	if(extensionProbabilityTrace) extensionProbabilityTrace->left.reset(hybridZ_left.size1(),hybridZ_left.size2());
	// temporary access
	const OutputConstraint & outConstraint = output.getOutputConstraint();
#if INTARNA_IN_DEBUG_MODE
	// check indices
	if (!energy.areComplementary(si1,si2) )
		throw std::runtime_error("PredictorMfeEns2dSeedExtension::fillHybridZ_left("+toString(si1)+","+toString(si2)+",..) are not complementary");
#endif

	// global vars to avoid reallocation
	size_t i1,i2,k1,k2;

	// determine whether or not lonely base pairs are allowed or if we have to
	// ensure a stacking to the right of the left boundary (i1,i2)
	const size_t noLpShift = outConstraint.noLP ? 1 : 0;
	Z_type iStackZ = Z_type(1);

	// iterate over all window starts i1 (seq1) and i2 (seq2)
	for (size_t l1=0; l1 < hybridZ_left.size1(); l1++) {
		for (size_t l2=0; l2 < hybridZ_left.size2(); l2++) {
			i1 = si1-l1;
			i2 = si2-l2;

			// referencing cell access
			Z_type & curZ = hybridZ_left(si1-i1,si2-i2);

			// init current cell (0 if not just right-most (j1,j2) base pair)
			curZ = (i1==si1 && i2==si2) ? energy.getBoltzmannWeight(energy.getE_init()) : 0.0;

			// check if complementary (use global sequence indexing)
			if( i1<si1
				&& i2<si2
				&& energy.areComplementary(i1,i2) )
			{

				// right-stacking of i if no-LP
				if (outConstraint.noLP) {
					// skip if no stacking possible
					if (!energy.areComplementary(i1+noLpShift,i2+noLpShift))
					{
						continue;
					}
					// get stacking energy to avoid recomputation in recursion below
					iStackZ = energy.getBoltzmannWeight(energy.getE_interLeft(i1,i1+noLpShift,i2,i2+noLpShift));
					// check just stacked
					curZ += iStackZ * hybridZ_left(l1-noLpShift,l2-noLpShift);
					if(extensionProbabilityTrace) extensionProbabilityTrace->left.add(l1,l2,l1-noLpShift,l2-noLpShift,iStackZ,{{i1,i2}});
				}

				// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
				for (k1=i1+noLpShift; k1++ < si1; ) {
					// ensure maximal loop length
					if (k1-i1-noLpShift > energy.getMaxInternalLoopSize1()+1) break;
					for (k2=i2+noLpShift; k2++ < si2; ) {
						// ensure maximal loop length
						if (k2-i2-noLpShift > energy.getMaxInternalLoopSize2()+1) break;
						// check if (k1,k2) are valid left boundary
						if ( ! Z_equal(hybridZ_left(si1-k1,si2-k2), 0.0) ) {
							curZ += (iStackZ
									* energy.getBoltzmannWeight(energy.getE_interLeft(i1+noLpShift,k1,i2+noLpShift,k2))
									* hybridZ_left(si1-k1,si2-k2));
							if(extensionProbabilityTrace) {
								std::vector<Interaction::BasePair> pairs{{i1,i2}};
								if(noLpShift) pairs.emplace_back(i1+1,i2+1);
								extensionProbabilityTrace->left.add(l1,l2,si1-k1,si2-k2,
										iStackZ*energy.getBoltzmannWeight(energy.getE_interLeft(i1+noLpShift,k1,i2+noLpShift,k2)),std::move(pairs));
							}
						}
					} // k2
				} // k1

				// correction for left seeds
				if ( i1<si1 && i2<si2 && seedHandler.isSeedBound(i1, i2) ) {

					// check if seed is to be processed:
					bool substractThisSeed =
							// check if left of anchor seed
										( i1+seedHandler.getSeedLength1(i1,i2)-1 <= si1
										&& i2+seedHandler.getSeedLength2(i1,i2)-1 <= si2 )
							// check if overlapping with anchor seed
									||	seedHandler.areLoopOverlapping(i1,i2,si1,si2);

					if (substractThisSeed) {

						// iterate seeds in S region
						size_t sj1 = RnaSequence::lastPos, sj2 = RnaSequence::lastPos;
						size_t si1overlap = si1+1, si2overlap = si2+1;
						// find left-most loop-overlapping seed
						while( seedHandler.updateToNextSeed(sj1,sj2
								, i1, std::min(si1,i1+seedHandler.getSeedLength1(i1, i2)-2)
								, i2, std::min(si2,i2+seedHandler.getSeedLength2(i1, i2)-2)) )
						{
							// check if right of i1,i2 and overlapping
							if (sj1 > i1 && seedHandler.areLoopOverlapping(i1, i2, sj1, sj2)) {
								// update left-most loop-overlapping seed
								if (sj1 < si1overlap) {
									si1overlap = sj1;
									si2overlap = sj2;
								}
							}
						}

						// if we found an overlapping seed
						if (si1overlap <= si1) {
							// check if right side is non-empty
							if ( ! E_equal(hybridZ_left( si1-si1overlap, si2-si2overlap),0) ) {
								// compute Energy of loop S \ S'
								E_type nonOverlapE = getNonOverlappingEnergy(i1, i2, si1overlap, si2overlap);
								// subtract energy.getBoltzmannWeight( nonOverlapE ) * hybridZ_left( si1overlap, si2overlap up to anchor seed [==1 if equal])
								Z_type correctionTerm = energy.getBoltzmannWeight( nonOverlapE )
														* hybridZ_left( si1-si1overlap, si2-si2overlap);
								curZ -= correctionTerm;
								if(extensionProbabilityTrace) extensionProbabilityTrace->left.add(l1,l2,si1-si1overlap,si2-si2overlap,
										-energy.getBoltzmannWeight(nonOverlapE),extensionSeedPrefix(i1,i2,si1overlap));
							// sanity insurance
								if (curZ < 0) {
									curZ = Z_type(0.0);
								if(extensionProbabilityTrace) extensionProbabilityTrace->left.edges[l1*extensionProbabilityTrace->left.columns+l2].clear();
								}
							}
						} else {
							// get data for seed to be removed
							const Z_type seedZ_rm = energy.getBoltzmannWeight(seedHandler.getSeedE(i1, i2));
							const size_t sj1_rm = i1+seedHandler.getSeedLength1(i1,i2)-1;
							const size_t sj2_rm = i2+seedHandler.getSeedLength2(i1,i2)-1;
							// if no S'
							// substract seedZ * hybridZ_left(right end seed up to anchor seed)
							Z_type correctionTerm = seedZ_rm
									* hybridZ_left( si1-sj1_rm
												  , si2-sj2_rm );
							if(extensionProbabilityTrace) extensionProbabilityTrace->left.add(l1,l2,si1-sj1_rm,si2-sj2_rm,
									-seedZ_rm,extensionSeedPrefix(i1,i2,sj1_rm));
							// if noLP : handle explicit loop right of current seed
							if (outConstraint.noLP) {
								for (k1=sj1_rm; k1++ < si1; ) {
									// ensure maximal loop length
									if (k1-sj1_rm > energy.getMaxInternalLoopSize1()+1) break;
									for (k2=sj2_rm; k2++ < si2; ) {
										// ensure at least one unpaired base in interior loop following the seed to be removed
										if (sj1_rm-k1 + sj2_rm-k2 == 2) {continue;}
										// ensure maximal loop length
										if (k2-sj2_rm > energy.getMaxInternalLoopSize2()+1) break;
										// check if (k1,k2) are valid left boundary
										if ( ! Z_equal(hybridZ_left(si1-k1,si2-k2), 0.0) ) {
											correctionTerm += (seedZ_rm
													* energy.getBoltzmannWeight(energy.getE_interLeft(sj1_rm,k1,sj2_rm,k2))
													* hybridZ_left(si1-k1,si2-k2) );
											if(extensionProbabilityTrace) extensionProbabilityTrace->left.add(l1,l2,si1-k1,si2-k2,
													-seedZ_rm*energy.getBoltzmannWeight(energy.getE_interLeft(sj1_rm,k1,sj2_rm,k2)),extensionSeedPrefix(i1,i2,k1));
										}
									} // k2
								} // k1

							}
							curZ -= correctionTerm;
							// sanity insurance
							if (curZ < 0) {
								curZ = Z_type(0.0);
								if(extensionProbabilityTrace) extensionProbabilityTrace->left.edges[l1*extensionProbabilityTrace->left.columns+l2].clear();
							}
						}
					} // substractThisSeed
				}

			} // complementary

		} // i2
	} // i1

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2dSeedExtension::
fillHybridZ_right( const size_t sj1, const size_t sj2 )
{
	if(extensionProbabilityTrace) extensionProbabilityTrace->right.reset(hybridZ_right.size1(),hybridZ_right.size2());
	// temporary access
	const OutputConstraint & outConstraint = output.getOutputConstraint();
#if INTARNA_IN_DEBUG_MODE
	// check indices
	if (!energy.areComplementary(sj1,sj2) )
		throw std::runtime_error("PredictorMfeEns2dSeedExtension::fillHybridZ_right("+toString(sj1)+","+toString(sj2)+",..) are not complementary");
#endif

	// global vars to avoid reallocation
	size_t j1,j2,k1,k2;

	// determine whether or not lonely base pairs are allowed or if we have to
	// ensure a stacking to the right of the left boundary (i1,i2)
	const size_t noLpShift = outConstraint.noLP ? 1 : 0;
	Z_type iStackZ = Z_type(1);

	// iterate over all window ends j1 (seq1) and j2 (seq2)
	for (j1=sj1; j1-sj1 < hybridZ_right.size1(); j1++ ) {
		for (j2=sj2; j2-sj2 < hybridZ_right.size2(); j2++ ) {

			// referencing cell access
			Z_type & curZ = hybridZ_right(j1-sj1,j2-sj2);

			// init partition function for current cell -> (i1,i2) are complementary per definition
			curZ = sj1==j1 && sj2==j2 ? energy.getBoltzmannWeight(0.0) : 0.0;

			// check if complementary free base pair
			if( sj1<j1
				&& sj2<j2
				&& energy.areComplementary(j1,j2) )
			{

				// left-stacking of j if no-LP
				if (outConstraint.noLP) {
					// skip if no stacking possible
					if (!energy.areComplementary(j1-noLpShift,j2-noLpShift))
					{
						continue;
					}
					// get stacking energy to avoid recomputation in recursion below
					iStackZ = energy.getBoltzmannWeight(energy.getE_interLeft(j1-noLpShift,j1,j2-noLpShift,j2));
					// check just stacked seed extension
					if (j1-noLpShift==sj1 && j2-noLpShift==sj2) {
						curZ += iStackZ * hybridZ_right(0,0);
						if(extensionProbabilityTrace) extensionProbabilityTrace->right.add(j1-sj1,j2-sj2,0,0,iStackZ,{{j1,j2}});
					}
				}

				// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
				for (k1=j1-noLpShift; k1-- > sj1; ) {
					// ensure maximal loop length
					if (j1-noLpShift-k1 > energy.getMaxInternalLoopSize1()+1) break;
					for (k2=j2-noLpShift; k2-- > sj2; ) {
						// ensure maximal loop length
						if (j2-noLpShift-k2 > energy.getMaxInternalLoopSize2()+1) break;
						// check if (k1,k2) are valid left boundary
						if ( ! Z_equal(hybridZ_right(k1-sj1,k2-sj2), 0.0) ) {
							// update partition function
							curZ += ( hybridZ_right(k1-sj1,k2-sj2)
									* energy.getBoltzmannWeight(energy.getE_interLeft(k1,j1-noLpShift,k2,j2-noLpShift))
									* iStackZ );
							if(extensionProbabilityTrace) {
								std::vector<Interaction::BasePair> pairs{{j1,j2}};
								if(noLpShift) pairs.emplace_back(j1-1,j2-1);
								extensionProbabilityTrace->right.add(j1-sj1,j2-sj2,k1-sj1,k2-sj2,
										energy.getBoltzmannWeight(energy.getE_interLeft(k1,j1-noLpShift,k2,j2-noLpShift))*iStackZ,std::move(pairs));
							}
						}
					} // k2
				} // k1
			}
		}
	}

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns2dSeedExtension::
traceBack( Interaction & interaction )
{

	// forward tracing
	PredictorMfeEns::traceBack(interaction);

	// add seeds in region
	seedHandler.addSeeds( interaction );

}

////////////////////////////////////////////////////////////////////////////


} // namespace
