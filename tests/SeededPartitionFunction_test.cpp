#include "catch.hpp"
#include "SeededChainOracle.h"
#include <cmath>

TEST_CASE("seeded chain oracle counts structures instead of seed witnesses", "[SeededPartitionFunction]") {
	using namespace seeded_oracle;
	Model m{3,3,3,3,false, {{{0,0},{1,1}},{{1,1},{2,2}}},
		[](Pair p){return p[0]==p[1];}, [](Pair){return 2.L;},
		[](Pair a, Pair b){return stack(a,b) ? 3.L : 0.L;},
		[](Pair,Pair){return 5.L;}};
	auto r=enumerate(m);
	REQUIRE(r.z == 150.L); // two two-pair chains (30 each), one three-pair (90)
	REQUIRE(r.mass[Pair{1,1}] == 150.L);
	REQUIRE(r.pairCount == 390.L);
	m.seeds={{{1,1}}};
	r=enumerate(m);
	REQUIRE(r.z == 160.L);
	m.noLP=true;
	REQUIRE(enumerate(m).z == 150.L);
	m.span1=1;
	REQUIRE(enumerate(m).z == 0.L);
}

#include "SeededTheoryFixtures.h"
TEST_CASE("unrounded fixed-seed theory gold agrees with exhaustive chains", "[SeededPartitionFunction]") {
	using namespace seeded_oracle;
	auto fixtures=theoryFixtures();
	REQUIRE(fixtures.size()==10);
	for (const auto & f:fixtures) {
		INFO(f.name);
		const auto got=enumerate(f.model);
		REQUIRE(double(got.z)==Approx(double(f.expected.z)).epsilon(2e-12));
		for (const auto & [pair, mass]:f.expected.mass) {
			auto it=got.mass.find(pair);
			REQUIRE(double(it==got.mass.end()?0:it->second)==Approx(double(mass)).epsilon(2e-12));
		}
		if (f.name=="random03") {
			REQUIRE(double(got.mass.at({9,4})/got.z)==Approx(double(std::exp(1.L)/(1+2*std::exp(1.L)))).epsilon(2e-12));
		}
	}
}
