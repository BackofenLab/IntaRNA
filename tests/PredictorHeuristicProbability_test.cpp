#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/AccessibilityBasePair.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/InteractionEnergyIdxOffset.h"
#include "IntaRNA/OutputHandlerInteractionList.h"
#include "IntaRNA/PredictorMfeEns2dHeuristic.h"
#include "IntaRNA/ReverseAccessibility.h"

#include <algorithm>
#include <array>
#include <functional>
#include <iterator>
#include <map>
#include <stdexcept>
#include <vector>

using namespace IntaRNA;

namespace {

using Pair=std::array<size_t,2>;
using Chain=std::vector<Pair>;

struct ChainCandidate {
	Chain chain;
	E_type hybrid;
};

struct HeuristicReference {
	Z_type z=0;
	std::map<Pair,Z_type> mass;
	size_t candidates=0;
};

bool adjacent(Pair p,Pair q) {
	return p[0]+1==q[0] && p[1]+1==q[1];
}

/** Enumerate complete chains first, then filter by the documented selected
 * suffix rule. Unlike the implementation, this stores complete pair lists,
 * sums integer chain energies and directly counts pairs in admitted chains.
 */
HeuristicReference enumerateHeuristic(const InteractionEnergy & energy,
		const OutputConstraint & constraint,size_t n,size_t m)
{
	std::map<Pair,std::vector<ChainCandidate>> chains;
	Chain chain;
	std::function<void(E_type)> enumerate=[&](E_type hybrid) {
		const Pair left=chain.front(),right=chain.back();
		if (right[0]-left[0]+1>energy.getAccessibility1().getMaxLength()
				|| right[1]-left[1]+1>energy.getAccessibility2().getMaxLength()) return;
		chains[left].push_back({chain,hybrid});
		for(size_t i=right[0]+1;i<n;++i) for(size_t j=right[1]+1;j<m;++j) {
			if (!energy.areComplementary(i,j)) continue;
			const E_type edge=energy.getE_interLeft(right[0],i,right[1],j);
			if (E_isINF(edge)) continue;
			chain.push_back({i,j});enumerate(hybrid+edge);chain.pop_back();
		}
	};
	for(size_t i=0;i<n;++i) for(size_t j=0;j<m;++j) {
		if (!energy.areComplementary(i,j)) continue;
		chain={{i,j}};enumerate(energy.getE_init());
	}
	HeuristicReference result;
	std::map<Pair,Chain> selected;
	for(size_t i=n;i-- >0;) for(size_t j=m;j-- >0;) {
		const Pair left{i,j};
		std::vector<std::pair<std::array<size_t,3>,const ChainCandidate *>> admitted;
		for(const auto & candidate:chains[left]) {
			const auto & c=candidate.chain;
			if (constraint.noGUend && energy.isGU(c.back()[0],c.back()[1])) continue;
			size_t suffix=1;
			std::array<size_t,3> order{0,0,0};
			if (constraint.noLP) {
				if (c.size()<2 || !adjacent(c[0],c[1])) continue;
				if (c.size()==2) suffix=c.size();
				else if (selected[c[1]]==Chain(c.begin()+1,c.end())) order={1,0,0};
				else {
					if (adjacent(c[1],c[2])) continue;
					suffix=2;order={2,c[2][0],c[2][1]};
				}
			} else if (c.size()>1) order={1,c[1][0],c[1][1]};
			if (suffix<c.size() && selected[c[suffix]]!=Chain(c.begin()+suffix,c.end())) continue;
			admitted.push_back({order,&candidate});
		}
		std::sort(admitted.begin(),admitted.end(),[](const auto & a,const auto & b){return a.first<b.first;});
		E_type best=E_INF;
		for(const auto & entry:admitted) {
			const auto & candidate=*entry.second;
			const Pair right=candidate.chain.back();
			const E_type full=energy.getE(i,right[0],j,right[1],candidate.hybrid);
			if (full<best) {best=full;selected[left]=candidate.chain;}
			if (constraint.noGUend && energy.isGU(i,j)) continue;
			if (constraint.maxED<Accessibility::ED_UPPER_BOUND
					&& (energy.getED1(i,right[0])>constraint.maxED || energy.getED2(j,right[1])>constraint.maxED)) continue;
			const Z_type weight=energy.getBoltzmannWeight(full);
			result.z+=weight;++result.candidates;
			for(const Pair pair:candidate.chain) result.mass[pair]+=weight;
		}
	}
	return result;
}

class FailingHeuristic : public PredictorMfeEns2dHeuristic {
public:
	using PredictorMfeEns2dHeuristic::PredictorMfeEns2dHeuristic;
	bool reject=false;
protected:
	void updateZ(size_t i1,size_t j1,size_t i2,size_t j2,Z_type z,bool hybrid) override {
		if (reject) throw std::runtime_error("test probability failure");
		PredictorMfeEns::updateZ(i1,j1,i2,j2,z,hybrid);
	}
};

class HeuristicTinyInitiation : public InteractionEnergyBasePair {
public:
	using InteractionEnergyBasePair::InteractionEnergyBasePair;
	E_type getE_init() const override { return Ekcal_2_E(60); }
};

}

TEST_CASE("heuristic pair probabilities count the retained candidate chains", "[BasePairProbabilities][PredictorHeuristicProbability]") {
	#include "testEasyLoggingSetup.icc"
	for(const auto & sequences:std::vector<std::pair<std::string,std::string>>{
			{"GGGGG","CCCCC"},{"GGAGG","CCUCC"},{"GUGGC","GCCUU"}}) {
		RnaSequence target("target",sequences.first),query("query",sequences.second);
		AccessibilityBasePair at(target,4,nullptr),aq(query,4,nullptr);ReverseAccessibility reverse(aq);
		InteractionEnergyBasePair energy(at,reverse,2,1);
		for(bool noLP:{false,true}) for(bool noGU:{false,true}) for(unsigned regionKind:{0u,1u,2u}) for(E_type maxED:{Accessibility::ED_UPPER_BOUND,E_type(50)}) {
			const IndexRange region=regionKind==0?IndexRange(0,4):regionKind==1?IndexRange(1,3):IndexRange(2,RnaSequence::lastPos);
			OutputConstraint constraint(1,OutputConstraint::OVERLAP_BOTH,E_INF,E_INF,false,noLP,noGU,true,true,maxED);
			InteractionEnergyIdxOffset shifted(energy,region.from,region.from);
			const size_t length=region.to==RnaSequence::lastPos?target.size()-region.from:region.to-region.from+1;
			const auto reference=enumerateHeuristic(shifted,constraint,length,length);
			OutputHandlerInteractionList enabled(constraint,1),disabled(constraint,1);
			BasePairProbabilities result(target.size(),query.size(),true);
			PredictorMfeEns2dHeuristic on(energy,enabled,nullptr,&result),off(energy,disabled,nullptr);
			on.predict(region,region);off.predict(region,region);result.finalize();
			REQUIRE(result.isApproximate());
			CAPTURE(sequences.first,sequences.second,noLP,noGU,regionKind,maxED);
			REQUIRE(on.getZall()==off.getZall());
			REQUIRE(result.getZ()==on.getZall());
			REQUIRE(result.getZ()==Approx(reference.z).epsilon(3e-12));
			Matrix<Z_type> expected(target.size(),query.size(),0);
			for(const auto & entry:reference.mass) {
				const auto bp=shifted.getBasePair(entry.first[0],entry.first[1]);
				expected(bp.first,bp.second)=entry.second;
			}
			for(size_t i=0;i<target.size();++i) for(size_t j=0;j<query.size();++j) {
				REQUIRE(result.rawMasses()(i,j)==Approx(expected(i,j)).epsilon(3e-12));
				REQUIRE(result.seedPairs()(i,j)==0);
			}
			REQUIRE(std::distance(enabled.begin(),enabled.end())==std::distance(disabled.begin(),disabled.end()));
			auto a=enabled.begin(),b=disabled.begin();
			for(;a!=enabled.end();++a,++b) {
				REQUIRE((*a)->energy==(*b)->energy);
				REQUIRE((*a)->basePairs==(*b)->basePairs);
			}
		}
	}
}

TEST_CASE("heuristic pair probabilities handle empty and singleton domains", "[BasePairProbabilities][PredictorHeuristicProbability]") {
	#include "testEasyLoggingSetup.icc"
	for(const std::string & queryText:{std::string("C"),std::string("A")}) for(bool noLP:{false,true}) {
		RnaSequence target("target","G"),query("query",queryText);
		AccessibilityDisabled at(target,0,nullptr),aq(query,0,nullptr);ReverseAccessibility reverse(aq);
		InteractionEnergyBasePair energy(at,reverse);
		OutputConstraint constraint(1,OutputConstraint::OVERLAP_BOTH,E_INF,E_INF,false,noLP,false,true,false);
		OutputHandlerInteractionList output(constraint,1);
		BasePairProbabilities result(1,1);
		PredictorMfeEns2dHeuristic predictor(energy,output,nullptr,&result);
		predictor.predict();result.finalize();
		REQUIRE(result.isApproximate());
		const bool empty=noLP || queryText=="A";
		REQUIRE((result.status()==BasePairProbabilities::Status::empty)==empty);
		if (!empty) REQUIRE(result.probabilities()(0,0)==1);
	}
}

TEST_CASE("heuristic pair owner merges disjoint regions and rejects failed results", "[BasePairProbabilities][PredictorHeuristicProbability]") {
	#include "testEasyLoggingSetup.icc"
	RnaSequence target("target","GGGGGG"),query("query","CCCCCC");
	AccessibilityDisabled at(target,0,nullptr),aq(query,0,nullptr);ReverseAccessibility reverse(aq);
	InteractionEnergyBasePair energy(at,reverse,2,2);
	OutputConstraint constraint(1,OutputConstraint::OVERLAP_BOTH,E_INF,E_INF,false,true,false,true,false);
	OutputHandlerInteractionList output(constraint,1);
	BasePairProbabilities result(6,6);
	FailingHeuristic predictor(energy,output,nullptr,&result);
	predictor.predict({0,2},{0,2});const Z_type first=predictor.getZall();
	predictor.predict({3,5},{3,5});result.finalize();
	REQUIRE(result.getZ()==2*first);
	REQUIRE(result.probabilities()(0,0)==0);
	BasePairProbabilities failed(6,6);OutputHandlerInteractionList failedOutput(constraint,1);
	FailingHeuristic failing(energy,failedOutput,nullptr,&failed);
	failing.predict({0,2},{0,2});failing.reject=true;
	REQUIRE_THROWS(failing.predict({3,5},{3,5}));
	REQUIRE(failed.status()==BasePairProbabilities::Status::failed);
	REQUIRE_THROWS(failed.finalize());
	BasePairProbabilities invalid(6,6);OutputHandlerInteractionList invalidOutput(constraint,1);
	FailingHeuristic badRange(energy,invalidOutput,nullptr,&invalid);
	REQUIRE_THROWS(badRange.predict({8,9},{0,2}));
	REQUIRE(invalid.status()==BasePairProbabilities::Status::failed);
}

TEST_CASE("heuristic probability collection preserves numerical-zero admission", "[BasePairProbabilities][PredictorHeuristicProbability]") {
	#include "testEasyLoggingSetup.icc"
	RnaSequence target("target","GGGGG"),query("query","CCCCC");
	AccessibilityDisabled at(target,0,nullptr),aq(query,0,nullptr);ReverseAccessibility reverse(aq);
	HeuristicTinyInitiation tiny(at,reverse,2,2,false,1,Ekcal_2_E(-40));
	InteractionEnergyBasePair vanishing(at,reverse,2,2,false,1,Ekcal_2_E(50));
	for(const InteractionEnergy * energy:{static_cast<InteractionEnergy *>(&tiny),static_cast<InteractionEnergy *>(&vanishing)}) for(bool noLP:{false,true}) {
		OutputConstraint constraint(1,OutputConstraint::OVERLAP_BOTH,E_INF,E_INF,false,noLP,false,true,false);
		OutputHandlerInteractionList baselineOut(constraint,1),probabilityOut(constraint,1);
		BasePairProbabilities result(5,5);
		PredictorMfeEns2dHeuristic baseline(*energy,baselineOut,nullptr),predictor(*energy,probabilityOut,nullptr,&result);
		baseline.predict();predictor.predict();result.finalize();
		REQUIRE(result.getZ()==baseline.getZall());
		REQUIRE(predictor.getZall()==baseline.getZall());
		if (energy==&vanishing) REQUIRE(result.status()==BasePairProbabilities::Status::empty);
	}
	OutputConstraint constraint(0,OutputConstraint::OVERLAP_BOTH,E_INF,E_INF,false,false,false,true,false);
	OutputHandlerInteractionList output(constraint,1);
	BasePairProbabilities failed(5,5);
	InteractionEnergyBasePair underflow(at,reverse,2,2,false,0.01,Ekcal_2_E(-1),3,Ekcal_2_E(1000));
	PredictorMfeEns2dHeuristic predictor(underflow,output,nullptr,&failed);
	REQUIRE_THROWS(predictor.predict());
	REQUIRE(failed.status()==BasePairProbabilities::Status::failed);
}
