#include "catch.hpp"
#include "IntaRNA/StackSeedDomain.h"
#include "IntaRNA/SeedHandlerExplicit.h"
#include "IntaRNA/SeedHandlerNoBulge.h"
#include "IntaRNA/SeedHandlerMfe.h"
#include "IntaRNA/SeedHandlerIdxOffset.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/AccessibilityDisabled.h"
#include <set>
using namespace IntaRNA;

TEST_CASE("stack capability and active occurrence domain", "[StackSeedDomain]") {
	#include "testEasyLoggingSetup.icc"
	RnaSequence t("t","GGGGGG"), q("q","CCCCCC");
	AccessibilityDisabled at(t,0,nullptr), aq(q,0,nullptr);
	ReverseAccessibility ar(aq);
	InteractionEnergyBasePair energy(at,ar);
	SeedConstraint computed(2,0,4,3,E_INF,Accessibility::ED_UPPER_BOUND,E_INF,IndexRangeList(),IndexRangeList(),"",false,false,false);
	SeedHandlerMfe mfe(energy,computed);
	SeedHandlerNoBulge stack(energy,computed);
	REQUIRE(mfe.guaranteesStackOnlySeeds());
	REQUIRE(stack.guaranteesStackOnlySeeds());
#if INTARNA_IN_DEBUG_MODE
	// Unsupported handlers are rejected only by the debug precondition check.
	SeedConstraint bulged(2,1,1,0,E_INF,Accessibility::ED_UPPER_BOUND,E_INF,IndexRangeList(),IndexRangeList(),"",false,false,false);
	SeedHandlerMfe bulge(energy,bulged);
	REQUIRE_FALSE(bulge.guaranteesStackOnlySeeds());
	REQUIRE_THROWS_AS(StackSeedDomain(bulge,6,6,6,6),std::invalid_argument);
#endif
	SeedConstraint explicitC(1,0,0,0,E_INF,Accessibility::ED_UPPER_BOUND,E_INF,
		IndexRangeList(),IndexRangeList(),"1|&6|,2|||&3|||,5||&1||",false,false,false);
	SeedHandlerIdxOffset explicitH(new SeedHandlerExplicit(energy,explicitC));
	REQUIRE(explicitH.guaranteesStackOnlySeeds());
	explicitH.setOffset1(1); explicitH.setOffset2(1);
	explicitH.fillSeed(0,3,0,3);
	StackSeedDomain domain(explicitH,4,4,4,4);
	REQUIRE(domain.seedLength(0,0)==3);
	REQUIRE(domain.seedLength(3,3)==0); // end outside the active range
	REQUIRE(domain.maxSeedLength()==3);
	StackSeedDomain mixed({{0,0,1},{2,2,2},{4,4,3}},6,6,4,3);
	REQUIRE(mixed.maxSeedLength()==2);
	std::set<Interaction::BasePair> starts;
	size_t calls=0;
	mixed.forEachStart([&](size_t i,size_t j){ ++calls; starts.emplace(i,j); });
	REQUIRE(calls==starts.size());
	for (size_t i=0;i<6;++i) for(size_t j=0;j<6;++j) {
		const bool expected=(i==0 && j==0) || (i<=2 && j>=1 && j<=2);
		REQUIRE((starts.count({i,j})!=0)==expected);
	}
	SeedConstraint explicitB(2,0,0,0,E_INF,Accessibility::ED_UPPER_BOUND,E_INF,
		IndexRangeList(),IndexRangeList(),"1|.|&4|.|",false,false,false);
	SeedHandlerExplicit eb(energy,explicitB);
	REQUIRE_FALSE(eb.guaranteesStackOnlySeeds());
}

TEST_CASE("stack-domain lookups and traversal preserve overlapping seed geometry", "[StackSeedDomain]") {
	// Different starts can share an ending pair; the shortest suffix must win
	// independently of insertion order, while duplicate seeds add no new starts.
	std::vector<StackSeedDomain::Occurrence> seeds{{3,2,1},{1,0,3},{2,1,2},{3,2,1},{0,4,2}};
	const std::vector<Interaction::BasePair> expected{{0,0},{0,1},{0,2},{0,3},{0,4},
		{1,0},{1,1},{1,2},{2,0},{2,1},{2,2},{3,0},{3,1},{3,2}};
	for (bool reverse:{false,true}) {
		if (reverse) std::reverse(seeds.begin(),seeds.end());
		StackSeedDomain domain(seeds,6,6,4,4);
		REQUIRE(domain.seedLength(1,0)==3);
		REQUIRE(domain.seedLength(2,1)==2);
		REQUIRE(domain.seedLength(3,2)==1);
		REQUIRE(domain.seedLength(0,1)==0); // swapped coordinates are a different pair
		REQUIRE(domain.seedEndingLength(3,2)==1);
		REQUIRE(domain.seedEndingLength(1,5)==2);
		REQUIRE(domain.seedEndingLength(2,3)==0);
		REQUIRE(domain.maxSeedLength()==3);
		std::vector<Interaction::BasePair> starts;
		domain.forEachStart([&](size_t i,size_t j){ starts.emplace_back(i,j); });
		REQUIRE(starts==expected);
	}
#if INTARNA_IN_DEBUG_MODE
	// Conflicting lengths violate a debug-checked caller precondition.
	REQUIRE_THROWS_AS((StackSeedDomain({{1,1,1},{1,1,2}},4,4,4,4)),std::invalid_argument);
#endif
}

TEST_CASE("stack-domain admission preserves handler-specific threshold equality", "[StackSeedDomain]") {
	#include "testEasyLoggingSetup.icc"
	RnaSequence t("t","GGG"),q("q","CCC");
	AccessibilityDisabled at(t,0,nullptr),aq(q,0,nullptr);ReverseAccessibility ar(aq);
	InteractionEnergyBasePair e(at,ar);
	for(bool hybrid:{false,true}) {
		SeedConstraint sc(2,0,0,0,hybrid?E_INF:-200,Accessibility::ED_UPPER_BOUND,hybrid?-200:E_INF,IndexRangeList(),IndexRangeList(),"",false,false,false);
		SeedHandlerMfe mfe(e,sc);SeedHandlerNoBulge nb(e,sc);
		REQUIRE(mfe.fillSeed(0,2,0,2)==4);
		REQUIRE(nb.fillSeed(0,2,0,2)==0);
		REQUIRE(StackSeedDomain(mfe,3,3,3,3).maxSeedLength()==2);
		REQUIRE(StackSeedDomain(nb,3,3,3,3).maxSeedLength()==0);
	}
	SeedConstraint exactED(2,0,0,0,E_INF,0,E_INF,IndexRangeList(),IndexRangeList(),"",false,false,false);
	SeedHandlerMfe mfe(e,exactED);SeedHandlerNoBulge nb(e,exactED);
	REQUIRE(mfe.fillSeed(0,2,0,2)==0);REQUIRE(nb.fillSeed(0,2,0,2)==0);
	// Explicit seeds ignore these computed-seed admission thresholds.
	SeedConstraint explicitC(1,0,0,0,-10000,0,-10000,IndexRangeList(),IndexRangeList(),"1|&3|,2||&1||",true,true,true);
	SeedHandlerExplicit ex(e,explicitC);REQUIRE(ex.fillSeed(0,2,0,2)==2);
	REQUIRE(StackSeedDomain(ex,3,3,3,3).maxSeedLength()==2);
}

TEST_CASE("explicit stack capability remains conservative for mixed encodings", "[StackSeedDomain]") {
	#include "testEasyLoggingSetup.icc"
	RnaSequence t("t","GGGG"),q("q","CCCC");
	AccessibilityDisabled at(t,0,nullptr),aq(q,0,nullptr);ReverseAccessibility ar(aq);
	InteractionEnergyBasePair e(at,ar);
	SeedConstraint sc(2,0,0,0,E_INF,Accessibility::ED_UPPER_BOUND,E_INF,
			IndexRangeList(),IndexRangeList(),"1|.|&2|.|,1|||&2|||",false,false,false);
	SeedHandlerExplicit handler(e,sc);
	// The retained minimum-energy pattern is a stack, but mixed input encodings
	// conservatively disable the stack-only handler contract before dispatch.
	REQUIRE_FALSE(handler.guaranteesStackOnlySeeds());
	REQUIRE(handler.fillSeed(0,3,0,3)==1);
	REQUIRE(handler.getSeedLength1(0,0)==3);
	REQUIRE(handler.getSeedLength2(0,0)==3);
#if INTARNA_IN_DEBUG_MODE
	REQUIRE_THROWS_AS(StackSeedDomain(handler,4,4,4,4),std::invalid_argument);
#endif
}
