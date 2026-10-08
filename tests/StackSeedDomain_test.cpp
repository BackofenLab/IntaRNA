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
	SeedConstraint bulged(2,1,1,0,E_INF,Accessibility::ED_UPPER_BOUND,E_INF,IndexRangeList(),IndexRangeList(),"",false,false,false);
	SeedHandlerMfe bulge(energy,bulged);
	REQUIRE_FALSE(bulge.guaranteesStackOnlySeeds());
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
	std::set<std::pair<size_t,size_t>> starts;
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
