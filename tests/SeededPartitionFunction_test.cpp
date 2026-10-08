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
		CAPTURE(f.name);
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

#include "IntaRNA/SeededPartitionFunction.h"
#include "IntaRNA/PartitionArithmetic.h"
#include "SeededStackKernel.h"
#include <random>
#include <chrono>
#include <iostream>

namespace {
using Kernel=IntaRNA::SeededPartitionFunction;
using IntaRNA::Z_type;
IntaRNA::StackSeedDomain domainFor(const seeded_oracle::Model & m) {
	std::vector<IntaRNA::StackSeedDomain::Occurrence> seeds;
	for(const auto & s:m.seeds) seeds.push_back({s[0][0],s[0][1],s.size()});
	return {seeds,m.n,m.m,m.span1,m.span2};
}
Kernel::Weights weightsFor(const seeded_oracle::Model & m) {
	return {m.valid,[&](auto p){return Z_type(m.init(p));},
		[&](auto p,auto q){return Z_type(m.edge(p,q));},
		[&](auto p,auto q){return Z_type(m.boundary(p,q));},{}};
}
void compareKernels(const seeded_oracle::Model & m) {
	const auto oracle=seeded_oracle::enumerate(m);
	const auto seeds=domainFor(m);
	Kernel::Domain d{m.n,m.m,m.span1,m.span2,m.n,m.m,m.noLP};
	for(bool suffix:{false,true}) {
		CAPTURE(suffix);
		auto w=weightsFor(m);
		std::map<seeded_oracle::Boundary,Z_type> boundaries;
		w.complete=[&](auto p,auto q,auto h,auto b){boundaries[{p[0],p[1],q[0],q[1]}]+=h*b;};
		auto got=suffix?Kernel::compute(d,seeds,w,true):seededStack(d,seeds,w,true);
		REQUIRE(double(got.z)==Approx(double(oracle.z)).epsilon(2e-12));
		for(auto [b,z]:oracle.boundaries) REQUIRE(double(boundaries[b])==Approx(double(z)).epsilon(2e-12));
		Z_type expectedCount=0;
		for(size_t i=0;i<m.n;++i) for(size_t j=0;j<m.m;++j) {
			auto it=oracle.mass.find({i,j});
			REQUIRE(double(got.mass(i,j))==Approx(double(it==oracle.mass.end()?0:it->second)).epsilon(2e-12));
			expectedCount+=got.mass(i,j);
		}
		REQUIRE(double(expectedCount)==Approx(double(oracle.pairCount)).epsilon(2e-12));
		auto only=suffix?Kernel::compute(d,seeds,w,false):seededStack(d,seeds,w,false);
		REQUIRE(only.z==got.z);
		REQUIRE(only.mass.storageSize()==0);
	}
}
}
TEST_CASE("both numerical kernels match all theory fixtures", "[SeededPartitionFunction]") {
	for(const auto & f:seeded_oracle::theoryFixtures()) { CAPTURE(f.name);compareKernels(f.model); }
}
TEST_CASE("heterogeneous weighted kernels match independent enumeration", "[SeededPartitionFunction]") {
	using namespace seeded_oracle;
	std::mt19937 rng(257);
	for(size_t trial=0;trial<100;++trial) {
		CAPTURE(trial);
		Model m{6,5,2+rng()%5,2+rng()%4,bool(trial%2),{},
			[trial](Pair p){return (p[0]+3*p[1]+trial)%7!=0;},
			[](Pair p){return .7L+.13L*p[0];},
			[trial](Pair p,Pair q){return (p[0]+q[1]+trial)%11==0?0.L:.1L+.17L*(q[0]-p[0])+.21L*(q[1]-p[1]);},
			[trial](Pair p,Pair q){return (q[0]+p[1]+trial)%5==0?0.L:.2L+.07L*(q[0]+p[1]);}};
		for(size_t i=0;i<m.n;++i) for(size_t j=0;j<m.m;++j) if(rng()%3==0) {
			size_t len=1+rng()%3;
			Chain seed;
			for(size_t k=0;k<len && i+k<m.n && j+k<m.m;++k) seed.push_back({i+k,j+k});
			m.seeds.push_back(seed);
		}
		compareKernels(m);
	}
}
TEST_CASE("checked partition arithmetic distinguishes empty and range failure", "[SeededPartitionFunction]") {
	using A=IntaRNA::PartitionArithmetic;
	REQUIRE(A::add(0,1e-100)==Z_type(1e-100));
	REQUIRE(A::multiply(0,1)==0);
	REQUIRE_THROWS_AS(A::exp(-1e6),std::range_error);
	REQUIRE_THROWS_AS(A::exp(1e6),std::range_error);
	REQUIRE_THROWS_AS(A::multiply(1e-200,1e-200),std::range_error);
	REQUIRE_THROWS_AS(A::divide(1e-200,1e200),std::range_error);
	REQUIRE_THROWS_AS(A::check(std::numeric_limits<Z_type>::quiet_NaN()),std::range_error);
	REQUIRE_THROWS_AS(A::add(std::numeric_limits<Z_type>::max(),std::numeric_limits<Z_type>::max()),std::range_error);
}
TEST_CASE("numerical seeded kernel comparison timing", "[.][SeededBenchmark]") {
	using namespace seeded_oracle;
	for (std::string kind:{"none","sparse","overlap","dense","long-stack","mixed","asymmetric"}) {
		const size_t n=kind=="long-stack"?90:45;
		Model m{n,n,kind=="asymmetric"?12u:30u,kind=="asymmetric"?24u:30u,false,{},
			[kind](Pair p){return kind=="long-stack"?p[0]==p[1]:(7*p[0]+11*p[1])%17!=0;},
			[](Pair){return .8L;},[](Pair p,Pair q){return q[0]-p[0]<=3 && q[1]-p[1]<=3?.7L:0.L;},
			[](Pair,Pair){return .6L;}};
		if(kind!="none") for(size_t i=0;i+4<n;++i) for(size_t j=0;j+4<n;++j) {
			if(kind=="sparse" && (i!=n/2 || j!=n/2)) continue;
			if((kind=="overlap" || kind=="long-stack") && i!=j) continue;
			size_t len=kind=="mixed"?1+(i+j)%4:4;Chain s;
			for(size_t k=0;k<len;++k) s.push_back({i+k,j+k});m.seeds.push_back(s);
		}
		const auto seeds=domainFor(m);const auto w=weightsFor(m);
		Kernel::Domain d{n,n,m.span1,m.span2,2,2,false};
		Z_type reference=0;
		for(bool suffix:{false,true}) for(bool outside:{false,true}) {
			const auto start=std::chrono::steady_clock::now();
			auto r=suffix?Kernel::compute(d,seeds,w,outside):seededStack(d,seeds,w,outside);
			const auto ms=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
			if(!suffix && !outside) reference=r.z; else REQUIRE(double(r.z)==Approx(double(reference)).epsilon(2e-12));
			std::cout<<kind<<","<<(suffix?"suffix":"stack")<<","<<outside<<","<<ms<<","<<r.z<<"\n";
		}
	}
}
