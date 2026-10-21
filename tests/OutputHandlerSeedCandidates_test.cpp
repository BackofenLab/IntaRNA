#include "catch.hpp"
#include "IntaRNA/OutputHandlerSeedCandidates.h"
#include "IntaRNA/OutputHandlerInteractionList.h"
using namespace IntaRNA;

TEST_CASE("seed quotas precede site deduplication and final output restrictions", "[OutputPerSeed]") {
	#include "testEasyLoggingSetup.icc"
	RnaSequence target("target","GGGGGGGGGG"),query("query","CCCCCCCCCC");
	auto candidate=[&](size_t i,size_t j,size_t q,size_t r,E_type e) {
		Interaction x(target,query);x.basePairs={{i,q},{j,r}};x.energy=e;return x;
	};
	const auto a=candidate(0,3,9,6,-1000),b=candidate(1,4,8,5,-900),c=candidate(6,7,3,2,-700),d=candidate(8,9,1,0,-600);
	for (size_t quota : {1,2}) for (auto overlap : {OutputConstraint::OVERLAP_BOTH,OutputConstraint::OVERLAP_NONE,OutputConstraint::OVERLAP_SEQ1,OutputConstraint::OVERLAP_SEQ2}) {
		CAPTURE(quota,overlap);
		OutputConstraint constraint(10,overlap,E_INF,E_INF,false,false,false,false,true,Accessibility::ED_UPPER_BOUND,quota);
		OutputHandlerSeedCandidates collector(constraint);
		// Equal boundaries from another seed retain the eligible better trace.
		auto alternative=a;alternative.energy=-950;alternative.basePairs.insert(alternative.basePairs.begin()+1,{2,7});
		collector.addSeed({0,9},b);collector.addSeed({0,9},a);
		collector.addSeed({1,8},alternative);collector.addSeed({1,8},d);
		// A second window rediscovers the same seed/site; it spends no new quota.
		collector.addSeed({0,9},a);collector.addSeed({6,3},c);
		OutputHandlerInteractionList output(constraint,10);collector.report(output);
		REQUIRE(output.begin()!=output.end());REQUIRE((**output.begin()).basePairs==a.basePairs);
		REQUIRE(output.reported()==(quota==2?(overlap==OutputConstraint::OVERLAP_BOTH?4:3):2));
	}
	SECTION("a rejected best site does not refill its seed quota") {
		OutputConstraint constraint(10,OutputConstraint::OVERLAP_NONE,E_INF,E_INF,false,false,false,false,true,Accessibility::ED_UPPER_BOUND,1);
		OutputHandlerSeedCandidates collector(constraint);collector.addSeed({0,9},a);
		collector.addSeed({1,8},b);collector.addSeed({1,8},d);
		OutputHandlerInteractionList output(constraint,10);collector.report(output);
		REQUIRE(output.reported()==1);
	}
	SECTION("energy window and final count are applied after the seed union") {
		OutputConstraint constraint(1,OutputConstraint::OVERLAP_BOTH,-800,50,false,false,false,false,true,Accessibility::ED_UPPER_BOUND,2);
		OutputHandlerSeedCandidates collector(constraint);collector.addSeed({0,9},a);collector.addSeed({6,3},c);
		OutputHandlerInteractionList output(constraint,10);collector.report(output);REQUIRE(output.reported()==1);
		collector.clear();OutputHandlerInteractionList empty(constraint,10);collector.report(empty);REQUIRE(empty.empty());
	}
}
