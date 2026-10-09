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

#include "IntaRNA/BasePairProbabilityWriter.h"
#include <sstream>
TEST_CASE("pair writer validates before emitting and checks stream failure", "[BasePairProbabilities]") {
	RnaSequence t("target","GA"),q("query","UC");
	BasePairProbabilities r(2,2);std::ostringstream out;
	REQUIRE_THROWS(BasePairProbabilityWriter::write(out,r,t,q));REQUIRE(out.str().empty());
	Matrix<Z_type> m(2,2,0);m(0,1)=1e-100;
	r.addRegion({0,1},{0,1},1,m);r.finalize();
	BasePairProbabilityWriter::write(out,r,t,q);
	REQUIRE(out.str().find("bpProb;U_1;C_2\nG_1;0;1e-100\nA_2;0;0\n")==0);
	std::ostringstream broken;broken.setstate(std::ios::badbit);
	REQUIRE_THROWS(BasePairProbabilityWriter::write(broken,r,t,q));
	BasePairProbabilities empty(2,2);empty.finalize();std::ostringstream na;
	BasePairProbabilityWriter::write(na,empty,t,q);
	REQUIRE(na.str()=="bpProb;U_1;C_2\nG_1;NA;NA\nA_2;NA;NA\n");
}

#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
TEST_CASE("seed annotations commit with regional masses and SVG requires success", "[BasePairProbabilities]") {
	#include "testEasyLoggingSetup.icc"
	RnaSequence t("target<&", "GGGG"),q("query", "CCCCC");
	AccessibilityDisabled at(t,0,nullptr),aq(q,0,nullptr);
	ReverseAccessibility ar(aq); InteractionEnergyBasePair energy(at,ar);
	BasePairProbabilities r(4,5,true);
	Matrix<Z_type> mass(2,2,0);mass(0,1)=2;mass(1,0)=2;
	Matrix<unsigned char> seeds(2,2,0);seeds(0,1)=1;
	r.addRegion({0,1},{1,2},2,mass,&seeds);
	r.addRegion({2,3},{3,4},2,mass,&seeds);
	REQUIRE(r.collectsSeedPairs());REQUIRE(r.seedPairs()(0,2)==1);
	REQUIRE(r.seedPairs()(2,4)==1);REQUIRE(r.seedPairs()(1,1)==0);
	std::ostringstream pending;
	REQUIRE_THROWS(BasePairProbabilityWriter::writeSvg(pending,r,energy));
	REQUIRE(pending.str().empty());
	r.finalize();std::ostringstream svg;
	BasePairProbabilityWriter::writeSvg(svg,r,energy);
	REQUIRE(svg.str().find("target&lt;&amp;")!=std::string::npos);
	REQUIRE(svg.str().find("Orange outline:")!=std::string::npos);
	std::ostringstream broken;broken.setstate(std::ios::badbit);
	REQUIRE_THROWS(BasePairProbabilityWriter::writeSvg(broken,r,energy));
	r.fail();std::ostringstream failed;
	REQUIRE_THROWS(BasePairProbabilityWriter::writeSvg(failed,r,energy));
	REQUIRE(failed.str().empty());
	BasePairProbabilities missing(4,5,true);
	REQUIRE_THROWS(missing.addRegion({0,1},{1,2},2,mass));
	REQUIRE(missing.status()==BasePairProbabilities::Status::failed);
	REQUIRE(missing.getZ()==0);REQUIRE(missing.seedPairs()(0,2)==0);
	BasePairProbabilities plain(4,5);
	REQUIRE_FALSE(plain.collectsSeedPairs());REQUIRE(plain.seedPairs().size1()==0);
	plain.finalize();std::ostringstream empty;
	BasePairProbabilityWriter::writeSvg(empty,plain,energy);
	REQUIRE(empty.str().find("data-probability=\"NA\"")!=std::string::npos);
	REQUIRE(empty.str().find("Orange outline:")==std::string::npos);
	BasePairProbabilities wrongSize(1,1);wrongSize.finalize();std::ostringstream wrong;
	REQUIRE_THROWS(BasePairProbabilityWriter::writeSvg(wrong,wrongSize,energy));
	REQUIRE(wrong.str().empty());
}
