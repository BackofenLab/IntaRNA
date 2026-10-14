#include "catch.hpp"
#include "SeededChainOracle.h"
#include <cmath>

TEST_CASE("seeded chain oracle counts structures instead of seed witnesses", "[SeededPartitionFunction]") {
	using namespace seeded_oracle;
	Model m{3,3,3,3,false, {{{0,0},{1,1}},{{1,1},{2,2}}},
		[](Pair p){return p.first==p.second;}, [](Pair){return 2.L;},
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
	for(const auto & s:m.seeds) seeds.push_back({s[0].first,s[0].second,s.size()});
	return {seeds,m.n,m.m,m.span1,m.span2};
}
Kernel::Weights weightsFor(const seeded_oracle::Model & m) {
	return {m.valid,[&](auto p){return Z_type(m.init(p));},
		[&](auto p,auto q){return Z_type(m.edge(p,q));},
		[&](auto p,auto q){return Z_type(m.boundary(p,q));},{},m.selectRight};
}
void compareKernels(const seeded_oracle::Model & m) {
	const auto oracle=seeded_oracle::enumerate(m);
	const auto seeds=domainFor(m);
	Kernel::Domain d{m.n,m.m,m.span1,m.span2,m.n,m.m,m.noLP};
	// Every allowed vertex as a singleton seed denotes the unconstrained domain.
	// Check the partition decomposition on all starts, including those pruned in
	// the original sparse-seed calculation.
	auto all=m;all.seeds.clear();all.selectRight={};
	for(size_t i=0;i<m.n;++i) for(size_t j=0;j<m.m;++j)
		if(m.valid({i,j})) all.seeds.push_back({{i,j}});
	const auto unconstrained=Kernel::compute(d,domainFor(all),weightsFor(all),false);
	auto exact=m;exact.selectRight={};
	const auto exactOracle=seeded_oracle::enumerate(exact);
	REQUIRE(double(unconstrained.z)==Approx(double(exactOracle.z+exactOracle.unseeded)).epsilon(2e-12));
	for(bool suffix:{false,true}) {
		if(!suffix && m.selectRight) continue; // maximal-stack gold is exact only
		CAPTURE(suffix);
		auto w=weightsFor(m);
		std::map<seeded_oracle::Boundary,Z_type> boundaries;
		std::map<seeded_oracle::Boundary,size_t> calls;
		w.complete=[&](auto p,auto q,auto h,auto b){
			const seeded_oracle::Boundary boundary{p.first,p.second,q.first,q.second};
			boundaries[boundary]+=h*b;
			++calls[boundary];
		};
		auto got=suffix?Kernel::compute(d,seeds,w,true):seededStack(d,seeds,w,true);
		REQUIRE(double(got.z)==Approx(double(oracle.z)).epsilon(2e-12));
		for(auto [b,z]:oracle.boundaries) {
			const auto found=boundaries.find(b);
			REQUIRE(double(found==boundaries.end()?0:found->second)==Approx(double(z)).epsilon(2e-12));
		}
		for(auto [b,z]:boundaries) {
			const auto expected=oracle.boundaries.find(b);
			REQUIRE(double(z)==Approx(double(expected==oracle.boundaries.end()?0:expected->second)).epsilon(2e-12));
			REQUIRE(calls[b]==1);
		}
		Z_type expectedCount=0;
		std::vector<Z_type> rows(m.n,0),columns(m.m,0);
		for(size_t i=0;i<m.n;++i) for(size_t j=0;j<m.m;++j) {
			auto it=oracle.mass.find({i,j});
			REQUIRE(double(got.mass(i,j))==Approx(double(it==oracle.mass.end()?0:it->second)).epsilon(2e-12));
			Z_type spot=0;
			for(auto [b,z]:boundaries)
				if(b[0]<=i && i<=b[2] && b[1]<=j && j<=b[3]) spot+=z;
			const auto expectedSpot=oracle.spots.find({i,j});
			REQUIRE(double(spot)==Approx(double(expectedSpot==oracle.spots.end()?0:expectedSpot->second)).epsilon(2e-12));
			REQUIRE(got.mass(i,j)>=0);
			REQUIRE(got.mass(i,j)<=spot+std::abs(got.z)*2e-12);
			REQUIRE(spot<=got.z+std::abs(got.z)*2e-12);
			rows[i]+=got.mass(i,j);
			columns[j]+=got.mass(i,j);
			expectedCount+=got.mass(i,j);
		}
		for(const auto mass:rows) REQUIRE(mass<=got.z+std::abs(got.z)*2e-12);
		for(const auto mass:columns) REQUIRE(mass<=got.z+std::abs(got.z)*2e-12);
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
			[trial](Pair p){return (p.first+3*p.second+trial)%7!=0;},
			[](Pair p){return .7L+.13L*p.first;},
			[trial](Pair p,Pair q){return (p.first+q.second+trial)%11==0?0.L:.1L+.17L*(q.first-p.first)+.21L*(q.second-p.second);},
			[trial](Pair p,Pair q){return (q.first+p.second+trial)%5==0?0.L:.2L+.07L*(q.first+p.second);}};
		for(size_t i=0;i<m.n;++i) for(size_t j=0;j<m.m;++j) if(rng()%3==0) {
			size_t len=1+rng()%3;
			Chain seed;
			for(size_t k=0;k<len && i+k<m.n && j+k<m.m;++k) seed.push_back({i+k,j+k});
			m.seeds.push_back(seed);
		}
		compareKernels(m);
		if(trial<40) {
			m.selectRight=[trial](Pair p,Pair q){return Pair{
				std::min(size_t(5),q.first+(p.first+trial)%3),
				std::min(size_t(4),q.second+(p.second+trial)%2)};};
			compareKernels(m);
		}
	}
}
TEST_CASE("first seed acceptance handles source hits overlapping and separated anchors", "[SeededPartitionFunction]") {
	using namespace seeded_oracle;
	// The unseeded prefixes and unrestricted suffixes include skipped diagonal
	// pairs, so an interior spot need not be an actual interaction pair.
	Model m{8,9,8,9,false,{},[](Pair p){return p.first==p.second;},
		[](Pair p){return 1.1L+.1L*p.first;},
		[](Pair p,Pair q){return stack(p,q)?.7L:.2L;},
		[](Pair p,Pair q){return .4L+.03L*p.first+.07L*q.second;}};
	const std::vector<std::pair<std::string,std::vector<Chain>>> cases{
		{"empty family",{}},
		{"source singleton",{{{0,0}}}},
		{"later singleton",{{{3,3}}}},
		{"two-pair anchor",{{{2,2},{3,3}}}},
		{"overlapping anchors",{{{1,1},{2,2}},{{2,2},{3,3}}}},
		{"adjacent anchors",{{{1,1},{2,2}},{{3,3},{4,4}}}},
		{"separated anchors",{{{0,0},{1,1}},{{5,5},{6,6}}}},
		{"filtered earlier stacks",{{{4,4},{5,5}}}},
		{"shared accepting pair",{{{1,1},{2,2},{3,3}},{{2,2},{3,3}},{{3,3}}}},
		{"length-seven overlap",{{{0,0},{1,1},{2,2},{3,3},{4,4},{5,5},{6,6}},
			{{1,1},{2,2},{3,3},{4,4},{5,5},{6,6},{7,7}}}}
	};
	for(const auto & [name,seeds]:cases) for(bool noLP:{false,true}) {
		CAPTURE(name,noLP);
		m.seeds=seeds;m.noLP=noLP;
		compareKernels(m);
	}
	// A singleton source is the entire interaction, and noLP rejects it.
	m.n=1;m.m=1;m.span1=1;m.span2=1;m.seeds={{{0,0}}};
	for(bool noLP:{false,true}) { m.noLP=noLP;compareKernels(m); }
	m.valid=[](Pair){return false;};
	compareKernels(m);
}

TEST_CASE("first-hit pair masses are derivatives of one frozen declared seed family", "[SeededPartitionFunction]") {
	using namespace seeded_oracle;
	Model m{6,7,5,6,false,{{{1,1},{2,2}},{{2,2},{3,3}},{{4,4}}},
		[](Pair p){return p.first==p.second || p.first+1==p.second;},
		[](Pair p){return .8L+.07L*p.first;},
		[](Pair p,Pair q){return .2L+.1L*(q.first-p.first)+.03L*(q.second-p.second);},
		[](Pair p,Pair q){return .4L+.02L*p.second+.05L*q.first;}};
	const auto seeds=domainFor(m);
	for(bool noLP:{false,true}) {
		const Kernel::Domain d{m.n,m.m,m.span1,m.span2,m.n,m.m,noLP};
		const auto baseline=Kernel::compute(d,seeds,weightsFor(m),true);
		for(size_t i=0;i<m.n;++i) for(size_t j=0;j<m.m;++j) {
			CAPTURE(noLP,i,j);
			const Pair marked{i,j};
			auto weights=weightsFor(m);
			// Every incoming edge owns its destination pair; initiation owns the
			// source pair. Doubling one pair gives Z(2)-Z(1)=its raw mass.
			weights.initiation=[&](Pair p){return Z_type(m.init(p)*(p==marked?2:1));};
			weights.transition=[&](Pair p,Pair q){return Z_type(m.edge(p,q)*(q==marked?2:1));};
			const auto markedZ=Kernel::compute(d,seeds,weights,false).z;
			REQUIRE(double(markedZ-baseline.z)==Approx(double(baseline.mass(i,j))).epsilon(2e-12).margin(2e-12));
		}
	}
}

TEST_CASE("frozen right extensions retain complete paths by their first seed owner", "[SeededPartitionFunction]") {
	using namespace seeded_oracle;
	Model m{7,8,7,8,false,{{{1,1},{2,2}},{{2,2},{3,3}},{{3,3}},{{4,4},{5,5}}},
		[](Pair p){return p.first==p.second || p.first+1==p.second;},
		[](Pair p){return .7L+.04L*p.first;},
		[](Pair p,Pair q){return stack(p,q)?.8L:.3L;},
		[](Pair p,Pair q){return .5L+.03L*p.first+.02L*q.second;}};
	for(bool noLP:{false,true}) for(size_t extension:{0,1}) {
		CAPTURE(noLP,extension);
		m.noLP=noLP;m.selectRight={};
		const auto exact=enumerate(m);
		m.selectRight=[extension](Pair,Pair end){return Pair{end.first+extension,end.second+extension};};
		const auto retained=enumerate(m);
		REQUIRE(retained.z>0);
		REQUIRE(retained.z<exact.z);
		compareKernels(m);
		const Kernel::Domain d{m.n,m.m,m.span1,m.span2,m.n,m.m,noLP};
		const auto baseline=Kernel::compute(d,domainFor(m),weightsFor(m),true);
		for(const Pair marked:{Pair{0,0},Pair{1,1},Pair{2,2},Pair{3,3},Pair{6,6}}) {
			// The retained family is fixed before inserting a pair fugacity.
			auto w=weightsFor(m);
			w.initiation=[&](Pair p){return Z_type(m.init(p)*(p==marked?2:1));};
			w.transition=[&](Pair p,Pair q){return Z_type(m.edge(p,q)*(q==marked?2:1));};
			const auto markedZ=Kernel::compute(d,domainFor(m),w,false).z;
			REQUIRE(double(markedZ-baseline.z)==Approx(double(baseline.mass(marked.first,marked.second))).epsilon(2e-12).margin(2e-12));
		}
	}
	// The full five-pair stack hits (1,1)-(2,2) first and is discarded by
	// that anchor. Its later seed must not admit the same path a second time.
	m={5,5,5,5,false,{{{1,1},{2,2}},{{2,2},{3,3}}},
		[](Pair p){return p.first==p.second;},[](Pair){return 2.L;},
		[](Pair p,Pair q){return stack(p,q)?3.L:0.L;},[](Pair,Pair){return 5.L;},
		[](Pair start,Pair end){return start.first==2?Pair{4,4}:end;}};
	const auto retained=enumerate(m);
	REQUIRE(retained.boundaries.count({0,0,4,4})==0);
	compareKernels(m);
}

TEST_CASE("seeded kernels preserve tiny mass and reject invalid weight arithmetic", "[SeededPartitionFunction]") {
	using namespace seeded_oracle;
	Model m{2,2,2,2,false,{{{0,0},{1,1}}},[](Pair p){return p.first==p.second;},
		[](Pair){return 1e-200L;},[](Pair,Pair){return .5L;},
		[](Pair,Pair){return 1e-40L;}};
	const auto seeds=domainFor(m);
	const Kernel::Domain d{2,2,2,2,0,0,false};
	const auto tiny=Kernel::compute(d,seeds,weightsFor(m),true);
	const Z_type expected=Z_type(1e-200)*Z_type(.5)*Z_type(1e-40);
	REQUIRE(tiny.z>0);
	REQUIRE(double(tiny.z/expected)==Approx(1).epsilon(2e-12));
	REQUIRE(tiny.mass(0,0)==tiny.z);
	REQUIRE(tiny.mass(1,1)==tiny.z);
	for(Z_type invalid:{Z_type(-1),std::numeric_limits<Z_type>::infinity(),std::numeric_limits<Z_type>::quiet_NaN()}) {
		CAPTURE(invalid);
		auto weights=weightsFor(m);
		weights.initiation=[invalid](Pair){return invalid;};
		REQUIRE_THROWS_AS(Kernel::compute(d,seeds,weights,true),std::range_error);
		weights=weightsFor(m);
		weights.transition=[invalid](Pair,Pair){return invalid;};
		REQUIRE_THROWS_AS(Kernel::compute(d,seeds,weights,true),std::range_error);
		weights=weightsFor(m);
		weights.boundary=[invalid](Pair,Pair){return invalid;};
		REQUIRE_THROWS_AS(Kernel::compute(d,seeds,weights,true),std::range_error);
	}
	for(bool overflow:{false,true}) {
		auto weights=weightsFor(m);
		weights.initiation=[overflow](Pair){return overflow?std::numeric_limits<Z_type>::max():std::numeric_limits<Z_type>::min();};
		weights.transition=[overflow](Pair,Pair){return overflow?Z_type(2):Z_type(.1);};
		REQUIRE_THROWS_AS(Kernel::compute(d,seeds,weights,false),std::range_error);
	}
	// A forbidden accepting edge is exact zero, never a tiny positive mass.
	auto zero=weightsFor(m);
	zero.transition=[](Pair,Pair){return Z_type(0);};
	zero.boundary=[](Pair,Pair)->Z_type{throw std::logic_error("unreachable boundary");};
	const auto empty=Kernel::compute(d,seeds,zero,true);
	REQUIRE(empty.z==0);
	REQUIRE(empty.mass(0,0)==0);
	REQUIRE(empty.mass(1,1)==0);
}

TEST_CASE("checked partition arithmetic distinguishes empty and range failure", "[SeededPartitionFunction]") {
	using A=IntaRNA::PartitionArithmetic;
	REQUIRE(A::add(0,1e-100)==Z_type(1e-100));
	REQUIRE(A::multiply(0,1)==0);
	REQUIRE_THROWS_AS(A::exp(-1e6),std::range_error);
	REQUIRE_THROWS_AS(A::exp(1e6),std::range_error);
	REQUIRE_THROWS_AS(A::multiply(std::numeric_limits<Z_type>::min(),Z_type(.1)),std::range_error);
	REQUIRE_THROWS_AS(A::divide(std::numeric_limits<Z_type>::min(),std::numeric_limits<Z_type>::max()),std::range_error);
	REQUIRE_THROWS_AS(A::check(std::numeric_limits<Z_type>::quiet_NaN()),std::range_error);
	REQUIRE_THROWS_AS(A::add(std::numeric_limits<Z_type>::max(),std::numeric_limits<Z_type>::max()),std::range_error);
}
TEST_CASE("numerical seeded kernel comparison timing", "[.][SeededBenchmark]") {
	using namespace seeded_oracle;
	for (std::string kind:{"none","sparse","overlap","dense","long-stack","mixed","asymmetric"}) {
		const size_t n=kind=="long-stack"?90:45;
		Model m{n,n,kind=="asymmetric"?12u:30u,kind=="asymmetric"?24u:30u,false,{},
			[kind](Pair p){return kind=="long-stack"?p.first==p.second:(7*p.first+11*p.second)%17!=0;},
			[](Pair){return .8L;},[](Pair p,Pair q){return q.first-p.first<=3 && q.second-p.second<=3?.7L:0.L;},
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
