#include "catch.hpp"

#include "IntaRNA/AccessibilityConstraint.h"
#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/BasePairProbabilities.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/OutputHandlerInteractionList.h"
#include "IntaRNA/PredictorMfeEns2d.h"
#include "IntaRNA/ReverseAccessibility.h"

#include <algorithm>
#include <functional>
#include <iterator>
#include <string>
#include <vector>

using namespace IntaRNA;

namespace {

struct UnseededOracle {
	Z_type z=0;
	Matrix<Z_type> mass;
	UnseededOracle(size_t n,size_t m) : mass(n,m,0) {}
};

// Enumerate increasing pair chains directly, without the predictor's suffix
// recurrence. Filter lonely pairs by inspecting both neighbours in the chain.
UnseededOracle enumerateUnseeded(const InteractionEnergy & energy,
		const OutputConstraint & constraint,const IndexRange & r1,const IndexRange & r2)
{
	UnseededOracle result(energy.size1(),energy.size2());
	const size_t end1=r1.to==RnaSequence::lastPos?energy.size1()-1:r1.to;
	const size_t end2=r2.to==RnaSequence::lastPos?energy.size2()-1:r2.to;
	std::vector<std::pair<size_t,size_t>> chain;
	std::function<void(E_type)> visit=[&](E_type hybrid) {
		const auto [i1,i2]=chain.front();
		const auto [j1,j2]=chain.back();
		bool accepted=true;
		if (constraint.noLP) {
			for(size_t k=0;k<chain.size();++k) {
				const bool left=k>0 && chain[k-1].first+1==chain[k].first
						&& chain[k-1].second+1==chain[k].second;
				const bool right=k+1<chain.size() && chain[k].first+1==chain[k+1].first
						&& chain[k].second+1==chain[k+1].second;
				if (!left && !right) accepted=false;
			}
		}
		if (constraint.noGUend && (energy.isGU(i1,i2) || energy.isGU(j1,j2))) accepted=false;
		if (energy.getED1(i1,j1)>constraint.maxED || energy.getED2(i2,j2)>constraint.maxED) accepted=false;
		const E_type complete=energy.getE(i1,j1,i2,j2,hybrid);
		if (accepted && E_isNotINF(complete)) {
			const Z_type weight=energy.getBoltzmannWeight(complete);
			result.z+=weight;
			for (const auto & p:chain) {
				const auto bp=energy.getBasePair(p.first,p.second);
				result.mass(bp.first,bp.second)+=weight;
			}
		}
		for(size_t k1=j1+1;k1<=end1;++k1) for(size_t k2=j2+1;k2<=end2;++k2) {
			if (!energy.isAccessible1(k1) || !energy.isAccessible2(k2)
					|| !energy.areComplementary(k1,k2)) continue;
			const E_type loop=energy.getE_interLeft(j1,k1,j2,k2);
			if (E_isINF(loop)) continue;
			chain.emplace_back(k1,k2);visit(hybrid+loop);chain.pop_back();
		}
	};
	for(size_t i1=r1.from;i1<=end1;++i1) for(size_t i2=r2.from;i2<=end2;++i2) {
		if (!energy.isAccessible1(i1) || !energy.isAccessible2(i2)
				|| !energy.areComplementary(i1,i2)) continue;
		chain.emplace_back(i1,i2);visit(energy.getE_init());chain.pop_back();
	}
	return result;
}

void compareUnseeded(const InteractionEnergy & energy,const OutputConstraint & constraint,
		const IndexRange & r1,const IndexRange & r2,bool compareOracle=true)
{
	OutputHandlerInteractionList ordinary(constraint,10), annotated(constraint,10);
	PredictorMfeEns2d baseline(energy,ordinary,nullptr);
	BasePairProbabilities result(energy.size1(),energy.size2(),true);
	PredictorMfeEns2d predictor(energy,annotated,nullptr,&result);
	baseline.predict(r1,r2);
	predictor.predict(r1,r2);
	result.finalize();
	REQUIRE(predictor.getZall()==baseline.getZall());
	REQUIRE(result.getZ()==baseline.getZall());
	REQUIRE(std::distance(ordinary.begin(),ordinary.end())==std::distance(annotated.begin(),annotated.end()));
	auto a=ordinary.begin(), b=annotated.begin();
	for(;a!=ordinary.end();++a,++b) {
		REQUIRE((*a)->energy==(*b)->energy);
		REQUIRE((*a)->basePairs==(*b)->basePairs);
	}
	for(size_t i=0;i<energy.size1();++i) for(size_t j=0;j<energy.size2();++j)
		REQUIRE(result.seedPairs()(i,j)==0);
	if (!compareOracle) return;
	const auto oracle=enumerateUnseeded(energy,constraint,r1,r2);
	REQUIRE(result.getZ()==Approx(oracle.z).epsilon(2e-12).margin(1e-14));
	for(size_t i=0;i<energy.size1();++i) for(size_t j=0;j<energy.size2();++j) {
		CAPTURE(i,j);
		REQUIRE(result.rawMasses()(i,j)==Approx(oracle.mass(i,j)).epsilon(2e-12).margin(1e-14));
	}
	if (oracle.z==0) REQUIRE(result.status()==BasePairProbabilities::Status::empty);
	else {
		const auto probabilities=result.probabilities();
		for(size_t i=0;i<energy.size1();++i) for(size_t j=0;j<energy.size2();++j)
			REQUIRE(probabilities(i,j)==Approx(oracle.mass(i,j)/oracle.z).epsilon(2e-12).margin(1e-14));
	}
}

class ProbabilityWidthAccessibility : public Accessibility {
public:
	ProbabilityWidthAccessibility(const RnaSequence & sequence,size_t maximum,E_type perNt)
		: Accessibility(sequence,maximum,nullptr),perNt(perNt) {}
	E_type getED(size_t from,size_t to) const override {
		checkIndices(from,to);
		return to-from+1>getMaxLength()+1?ED_UPPER_BOUND:E_type((to-from+1)*perNt);
	}
private:
	E_type perNt;
};

class TinyInitiationEnergy : public InteractionEnergyBasePair {
public:
	using InteractionEnergyBasePair::InteractionEnergyBasePair;
	E_type getE_init() const override { return Ekcal_2_E(60); }
};

class ChangedObjectivePredictor : public PredictorMfeEns2d {
public:
	using PredictorMfeEns2d::PredictorMfeEns2d;
protected:
	void updateZ(size_t i1,size_t j1,size_t i2,size_t j2,Z_type z,bool hybrid) override {
		PredictorMfeEns::updateZ(i1,j1,i2,j2,2*z,hybrid);
	}
};

}

TEST_CASE("exact unseeded pair probabilities enumerate physical pair chains", "[UnseededProbability][PredictorTinyOracle]") {
	#include "testEasyLoggingSetup.icc"
	for (const auto & sequence:std::vector<std::pair<std::string,std::string>>{
			{"GGGGGG","CCCCCC"},{"GUGCA","UGCGC"},{"AAGGGAA","ACCCUAA"},
			{"AAAA","AAAA"},{"G","C"}}) {
		for(bool noLP:{false,true}) for(bool noGU:{false,true}) for(size_t loop:{0,2}) {
			CAPTURE(sequence.first,sequence.second,noLP,noGU,loop);
			RnaSequence target("target",sequence.first), query("query",sequence.second);
			AccessibilityDisabled targetAcc(target,0,nullptr), queryAcc(query,0,nullptr);
			ReverseAccessibility reverse(queryAcc);
			InteractionEnergyBasePair energy(targetAcc,reverse,loop,loop,false,0.7,Ekcal_2_E(-0.8),3,Ekcal_2_E(0.13));
			OutputConstraint constraint(4,OutputConstraint::OVERLAP_BOTH,E_INF,E_INF,false,noLP,noGU,true,true);
			compareUnseeded(energy,constraint,{0,RnaSequence::lastPos},{0,RnaSequence::lastPos});
		}
	}
}

TEST_CASE("unseeded pair probabilities preserve filtered offset ensembles", "[UnseededProbability]") {
	#include "testEasyLoggingSetup.icc"
	RnaSequence target("target","GUGGGCG"),query("query","CUCCCCC");
	for(bool noLP:{false,true}) for(bool noGU:{false,true}) {
		ProbabilityWidthAccessibility targetAcc(target,3,Ekcal_2_E(0.10));
		ProbabilityWidthAccessibility queryAcc(query,4,Ekcal_2_E(0.15));
		ReverseAccessibility reverse(queryAcc);
		InteractionEnergyBasePair energy(targetAcc,reverse,1,2,false,1,Ekcal_2_E(-1),3,0,true,false);
		OutputConstraint constraint(5,OutputConstraint::OVERLAP_NONE,E_INF,E_INF,false,noLP,noGU,true,true,Ekcal_2_E(0.45));
		compareUnseeded(energy,constraint,{1,5},{2,6});
		compareUnseeded(energy,constraint,{2,RnaSequence::lastPos},{1,RnaSequence::lastPos});
	}
	AccessibilityConstraint blocked(target,"...b...",0,"","","");
	AccessibilityDisabled targetAcc(target,0,&blocked),queryAcc(query,0,nullptr);
	ReverseAccessibility reverse(queryAcc);
	InteractionEnergyBasePair energy(targetAcc,reverse,2,2);
	OutputConstraint constraint(1,OutputConstraint::OVERLAP_BOTH,E_INF,E_INF,false,false,false,true,true);
	compareUnseeded(energy,constraint,{1,5},{1,5});

	BasePairProbabilities result(target.size(),query.size());
	OutputHandlerInteractionList output(constraint,4);
	PredictorMfeEns2d predictor(energy,output,nullptr,&result);
	predictor.predict({0,2},{0,2});
	predictor.predict({4,6},{3,6});
	result.finalize();
	const auto first=enumerateUnseeded(energy,constraint,{0,2},{0,2});
	const auto second=enumerateUnseeded(energy,constraint,{4,6},{3,6});
	REQUIRE(result.getZ()==Approx(first.z+second.z).epsilon(2e-12));
	for(size_t i=0;i<target.size();++i) for(size_t j=0;j<query.size();++j)
		REQUIRE(result.rawMasses()(i,j)==Approx(first.mass(i,j)+second.mass(i,j)).epsilon(2e-12).margin(1e-14));
}

TEST_CASE("unseeded probability collection follows legacy numerical gates and fails atomically", "[UnseededProbability]") {
	#include "testEasyLoggingSetup.icc"
	RnaSequence target("target","GGGGG"),query("query","CCCCC");
	AccessibilityDisabled targetAcc(target,0,nullptr), queryAcc(query,0,nullptr);
	ReverseAccessibility reverse(queryAcc);
	for(bool noLP:{false,true}) {
		OutputConstraint constraint(2,OutputConstraint::OVERLAP_BOTH,E_INF,E_INF,false,noLP,false,true,true);
		TinyInitiationEnergy energy(targetAcc,reverse,2,2,false,1,Ekcal_2_E(-40));
		compareUnseeded(energy,constraint,{0,4},{0,4},false);
		InteractionEnergyBasePair vanishing(targetAcc,reverse,2,2,false,1,Ekcal_2_E(50));
		compareUnseeded(vanishing,constraint,{0,4},{0,4},false);
	}
	OutputConstraint constraint(0,OutputConstraint::OVERLAP_BOTH,E_INF,E_INF,false,false,false,true,true);
	OutputHandlerInteractionList output(constraint,1);
	InteractionEnergyBasePair energy(targetAcc,reverse);
	BasePairProbabilities unsupported(5,5);
	ChangedObjectivePredictor changed(energy,output,nullptr,&unsupported);
	REQUIRE_THROWS(changed.predict());
	REQUIRE(unsupported.status()==BasePairProbabilities::Status::failed);
	BasePairProbabilities invalid(5,5);
	PredictorMfeEns2d predictor(energy,output,nullptr,&invalid);
	REQUIRE_THROWS(predictor.predict({5,6},{0,4}));
	REQUIRE(invalid.status()==BasePairProbabilities::Status::failed);
	BasePairProbabilities overflow(5,5);
	InteractionEnergyBasePair extreme(targetAcc,reverse,2,2,false,1,Ekcal_2_E(-1000));
	PredictorMfeEns2d huge(extreme,output,nullptr,&overflow);
	REQUIRE_THROWS(huge.predict());
	REQUIRE(overflow.status()==BasePairProbabilities::Status::failed);
}
