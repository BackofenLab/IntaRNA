#include "catch.hpp"
#include "IntaRNA/BasePairProbabilities.h"
using namespace IntaRNA;
TEST_CASE("raw pair regions commit once and finalize only on success", "[BasePairProbabilities]") {
	BasePairProbabilities r(4,5);
	Matrix<Z_type> m(2,2,0);m(0,1)=2;m(1,0)=2;
	r.addRegion({0,1},{1,2},2,m);
	REQUIRE_THROWS(r.probabilities());
	r.addRegion({2,3},{3,4},2,m);
	r.finalize();
	REQUIRE(r.getZ()==4);
	REQUIRE(r.status()==BasePairProbabilities::Status::nonempty);
	auto p=r.probabilities();
	REQUIRE(p(0,2)==.5);REQUIRE(p(3,3)==.5);REQUIRE(p(0,4)==0);
	BasePairProbabilities empty(1,1);empty.finalize();
	REQUIRE(empty.status()==BasePairProbabilities::Status::empty);
	BasePairProbabilities overlap(4,5);overlap.addRegion({0,1},{1,2},2,m);
	REQUIRE_THROWS(overlap.addRegion({1,2},{2,3},2,m));
	REQUIRE(overlap.status()==BasePairProbabilities::Status::failed);
	REQUIRE(overlap.getZ()==2);REQUIRE_THROWS(overlap.finalize());
	BasePairProbabilities cancelled(4,5);cancelled.addRegion({0,1},{1,2},2,m);cancelled.fail();
	REQUIRE_THROWS(cancelled.finalize());REQUIRE_THROWS(cancelled.probabilities());
}
TEST_CASE("pair result rejects inconsistent and out-of-range numerical results", "[BasePairProbabilities]") {
	Matrix<Z_type> m(1,1,1e-100);BasePairProbabilities tiny(1,1);
	tiny.addRegion({0,0},{0,0},1e-100,m);tiny.finalize();REQUIRE(tiny.probabilities()(0,0)==1);
	m(0,0)=2;BasePairProbabilities excess(1,1);excess.addRegion({0,0},{0,0},1,m);
	REQUIRE_THROWS(excess.finalize());REQUIRE(excess.status()==BasePairProbabilities::Status::failed);
	m(0,0)=std::numeric_limits<Z_type>::max();BasePairProbabilities overflow(2,1);
	overflow.addRegion({0,0},{0,0},m(0,0),m);
	REQUIRE_THROWS(overflow.addRegion({1,1},{0,0},m(0,0),m));
	REQUIRE(overflow.status()==BasePairProbabilities::Status::failed);
}
