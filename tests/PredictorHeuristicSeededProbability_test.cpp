#include "catch.hpp"
#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/OutputHandlerInteractionList.h"
#include "IntaRNA/PredictorMfeEns2dHeuristicSeedExtension.h"
#include "IntaRNA/SeedHandlerExplicit.h"
#include "IntaRNA/SeedHandlerMfe.h"
#include <array>
#include <cmath>
#include <map>
#include <limits>

using namespace IntaRNA;

namespace {
// Introduce an independent fugacity x at one actual pair: every loop pays for
// its left pair, and the complete boundary pays for the last pair. Every
// physical chain is multilinear in x, including signed overlap corrections.
// Consequently Z(x)=Z_without+x*M_pair and a finite difference gives M exactly,
// without consulting the reverse tape or its pair-ownership bookkeeping.
class PairFugacityEnergy: public InteractionEnergyBasePair {
public:
	using InteractionEnergyBasePair::InteractionEnergyBasePair;
	size_t marked1=RnaSequence::lastPos,marked2=RnaSequence::lastPos;
	E_type penalty=0;
	E_type getE_interLeft(size_t i1,size_t j1,size_t i2,size_t j2) const override {
		const E_type e=InteractionEnergyBasePair::getE_interLeft(i1,j1,i2,j2);
		return E_isINF(e)?e:e+((i1==marked1 && i2==marked2)?penalty:0);
	}
	E_type getE(size_t i1,size_t j1,size_t i2,size_t j2,E_type hybrid) const override {
		const E_type e=InteractionEnergyBasePair::getE(i1,j1,i2,j2,hybrid);
		return E_isINF(e)?e:e+((j1==marked1 && j2==marked2)?penalty:0);
	}
};

class FixedChoiceSeededHeuristic: public PredictorMfeEns2dHeuristicSeedExtension {
public:
	using PredictorMfeEns2dHeuristicSeedExtension::PredictorMfeEns2dHeuristicSeedExtension;
	using Key=std::array<size_t,2>;
	struct Choice { size_t j1,j2; Z_type e; };
	using Choices=std::map<Key,Choice>;
	Choices choices;
	const Choices * fixed=nullptr;
protected:
	void fillHybridZ_right(size_t sj1,size_t sj2,size_t si1,size_t si2) override {
		PredictorMfeEns2dHeuristicSeedExtension::fillHybridZ_right(sj1,sj2,si1,si2);
		// Condition on the original heuristic's chosen domain; a fugacity is a
		// counting variable and must not change the domain being differentiated.
		if(fixed) {
			const auto & c=fixed->at({si1,si2}); j1opt=c.j1; j2opt=c.j2; E_right_opt=c.e;
		} else choices[{si1,si2}]={j1opt,j2opt,E_right_opt};
	}
};

class CorruptSeededHeuristicPartition: public PredictorMfeEns2dHeuristicSeedExtension {
public:
	using PredictorMfeEns2dHeuristicSeedExtension::PredictorMfeEns2dHeuristicSeedExtension;
	Z_type injected=0;
protected:
	void updateZ(size_t i1,size_t j1,size_t i2,size_t j2,Z_type partZ,bool isHybridZ) override;
};

void CorruptSeededHeuristicPartition::updateZ(size_t i1,size_t j1,size_t i2,size_t j2,
		Z_type partZ,bool isHybridZ)
{
	PredictorMfeEns2dHeuristicSeedExtension::updateZ(i1,j1,i2,j2,partZ,isHybridZ);
	Zall=injected;
}

SeedConstraint explicitSeeds(const std::string & patterns,size_t nominal=2) {
	return SeedConstraint(nominal,6,3,3,E_INF,E_INF,E_INF,IndexRangeList(),IndexRangeList(),patterns,true,true,true);
}
OutputConstraint constraints(bool noLP,bool noGU=false) {
	return OutputConstraint(0,OutputConstraint::OVERLAP_BOTH,E_INF,E_INF,false,noLP,noGU,true,false);
}
}

TEST_CASE("seeded heuristic outside agrees with independent pair fugacities", "[BasePairProbabilities][HeuristicSeededProbability]") {
	#include "testEasyLoggingSetup.icc"
	RnaSequence t("t","GGGGGG"),q("q","CCCCCC");
	AccessibilityDisabled at(t,0,nullptr),aq(q,0,nullptr); ReverseAccessibility ar(aq);
	PairFugacityEnergy energy(at,ar,2,2);
	// Isolated and competing seeds exercise left/right/both objectives,
	// overlapping and nonoverlapping subtraction, and bulged actual pairs.
	for(const auto & patterns:{"3||&3||","2||&4||,3||&3||","1||&5||,4||&2||",
			"2|.|&3||","1|.||&3|||,3||&3||","2|||&3|||,4||&2||"}) {
		const auto sc=explicitSeeds(patterns);
		for(bool noLP:{false,true}) {
			CAPTURE(patterns,noLP);
			const auto oc=constraints(noLP);
			energy.penalty=0;
			BasePairProbabilities result(t.size(),q.size(),true);
			OutputHandlerInteractionList output(oc,0);
			FixedChoiceSeededHeuristic predictor(energy,output,nullptr,new SeedHandlerExplicit(energy,sc),&result);
			predictor.predict(); result.finalize();
			REQUIRE(result.isApproximate());
			REQUIRE(result.getZ()==Approx(predictor.getZall()).epsilon(2e-12));
			const Z_type baseline=result.getZ();
			for(size_t i=0;i<t.size();++i) for(size_t j=0;j<q.size();++j) {
				CAPTURE(i,j);
				energy.marked1=i;energy.marked2=j;energy.penalty=100;
				OutputHandlerInteractionList changedOutput(oc,0);
				FixedChoiceSeededHeuristic changed(energy,changedOutput,nullptr,new SeedHandlerExplicit(energy,sc));
				changed.fixed=&predictor.choices;
				changed.predict();
				const Z_type expected=(baseline-changed.getZall())/(-std::expm1(-1.0));
				const auto bp=energy.getBasePair(i,j);
				REQUIRE(result.rawMasses()(bp.first,bp.second)==Approx(expected).epsilon(2e-10).margin(2e-10));
			}
		}
	}
}

TEST_CASE("seeded heuristic probabilities handle original coordinates regions and empty domains", "[BasePairProbabilities][HeuristicSeededProbability]") {
	#include "testEasyLoggingSetup.icc"
	RnaSequence t("t","GGGGGG"),q("q","CCCCCC");
	AccessibilityDisabled at(t,0,nullptr),aq(q,0,nullptr);ReverseAccessibility ar(aq);
	InteractionEnergyBasePair energy(at,ar,2,2);
	const auto oc=constraints(false);
	const auto sc=explicitSeeds("1||&5||,4||&2||",5);
	BasePairProbabilities result(6,6,true);OutputHandlerInteractionList output(oc,0);
	FixedChoiceSeededHeuristic predictor(energy,output,nullptr,new SeedHandlerExplicit(energy,sc),&result);
	predictor.predict({0,2},{0,2});const Z_type first=result.getZ();
	REQUIRE(first>0);
	predictor.predict({3,RnaSequence::lastPos},{3,RnaSequence::lastPos});result.finalize();
	REQUIRE(result.getZ()==Approx(2*first).epsilon(2e-12));
	REQUIRE(result.seedPairs()(0,5)==1);REQUIRE(result.seedPairs()(1,4)==1);
	REQUIRE(result.seedPairs()(3,2)==1);REQUIRE(result.seedPairs()(4,1)==1);
	REQUIRE(result.probabilities()(0,0)==0);
	for(const std::string pattern:{"2||&4||","1|||&4|||"}) {
		const auto shortSc=explicitSeeds(pattern);
		BasePairProbabilities empty(6,6,true);OutputHandlerInteractionList emptyOutput(oc,0);
		FixedChoiceSeededHeuristic p(energy,emptyOutput,nullptr,new SeedHandlerExplicit(energy,shortSc),&empty);
		p.predict({0,1},{0,1});empty.finalize();
		REQUIRE(empty.status()==BasePairProbabilities::Status::empty);
	}
	const auto singleton=explicitSeeds("2|&4|",2);
	BasePairProbabilities failed(6,6);OutputHandlerInteractionList failedOutput(oc,0);
	FixedChoiceSeededHeuristic bad(energy,failedOutput,nullptr,new SeedHandlerExplicit(energy,singleton),&failed);
	REQUIRE_THROWS(bad.predict());REQUIRE(failed.status()==BasePairProbabilities::Status::failed);
}

TEST_CASE("computed bulged seeded heuristic marginals condition on the selected seed family", "[BasePairProbabilities][HeuristicSeededProbability]") {
	#include "testEasyLoggingSetup.icc"
	RnaSequence t("t","GGAGGCG"),q("q","CCGCACC");
	AccessibilityDisabled at(t,0,nullptr),aq(q,0,nullptr);ReverseAccessibility ar(aq);
	PairFugacityEnergy energy(at,ar,2,2);
	const SeedConstraint sc(3,2,1,1,E_INF,E_INF,E_INF,IndexRangeList(),IndexRangeList(),"",true,true,true);
	// Freeze the concrete MFE seed selected at every start by encoding its
	// actual pairs explicitly. Perturbations then cannot select another seed.
	SeedHandlerMfe selected(energy,sc);selected.fillSeed(0,t.size()-1,0,q.size()-1);
	std::string patterns;
	size_t i=RnaSequence::lastPos,j=RnaSequence::lastPos;
	while(selected.updateToNextSeed(i,j,0,t.size()-1,0,q.size()-1)) {
		const size_t sl1=selected.getSeedLength1(i,j),sl2=selected.getSeedLength2(i,j);
		Interaction seed(t,q);seed.basePairs.push_back(energy.getBasePair(i,j));
		selected.traceBackSeed(seed,i,j);seed.basePairs.push_back(energy.getBasePair(i+sl1-1,j+sl2-1));
		const size_t qStart=q.size()-j-sl2;
		std::string ts(sl1,'.'),qs(sl2,'.');
		for(const auto & bp:seed.basePairs) { ts[bp.first-i]='|';qs[bp.second-qStart]='|'; }
		if(!patterns.empty()) patterns+=',';
		patterns+=std::to_string(i+1)+ts+'&'+std::to_string(qStart+1)+qs;
	}
	REQUIRE_FALSE(patterns.empty());
	const auto fixedSeeds=explicitSeeds(patterns,3);
	for(bool noLP:{false,true}) {
		const auto oc=constraints(noLP);energy.penalty=0;
		BasePairProbabilities result(t.size(),q.size());OutputHandlerInteractionList output(oc,0);
		FixedChoiceSeededHeuristic predictor(energy,output,nullptr,new SeedHandlerMfe(energy,sc),&result);
		predictor.predict();result.finalize();
		for(size_t a=0;a<t.size();++a) for(size_t b=0;b<q.size();++b) {
			CAPTURE(noLP,a,b);
			energy.marked1=a;energy.marked2=b;energy.penalty=100;
			OutputHandlerInteractionList changedOutput(oc,0);
			FixedChoiceSeededHeuristic changed(energy,changedOutput,nullptr,new SeedHandlerExplicit(energy,fixedSeeds));
			changed.fixed=&predictor.choices;changed.predict();
			const Z_type expected=(result.getZ()-changed.getZall())/(-std::expm1(-1.0));
			const auto bp=energy.getBasePair(a,b);
			REQUIRE(result.rawMasses()(bp.first,bp.second)==Approx(expected).epsilon(2e-10).margin(2e-10));
		}
	}
}

TEST_CASE("seeded heuristic rejects nonfinite reported partitions before committing masses", "[BasePairProbabilities][HeuristicSeededProbability]") {
	#include "testEasyLoggingSetup.icc"
	RnaSequence t("t","GGGG"),q("q","CCCC");
	AccessibilityDisabled at(t,0,nullptr),aq(q,0,nullptr);ReverseAccessibility ar(aq);
	InteractionEnergyBasePair energy(at,ar,2,2);
	const auto oc=constraints(false);
	const auto sc=explicitSeeds("1||&3||");
	for(Z_type injected:{std::numeric_limits<Z_type>::infinity(),std::numeric_limits<Z_type>::quiet_NaN()}) {
		CAPTURE(injected);
		BasePairProbabilities result(t.size(),q.size());OutputHandlerInteractionList output(oc,0);
		CorruptSeededHeuristicPartition predictor(energy,output,nullptr,new SeedHandlerExplicit(energy,sc),&result);
		predictor.injected=injected;
		REQUIRE_THROWS_AS(predictor.predict(),std::range_error);
		REQUIRE(result.status()==BasePairProbabilities::Status::failed);
		REQUIRE(result.getZ()==0);
		REQUIRE_THROWS(result.finalize());
	}
}
