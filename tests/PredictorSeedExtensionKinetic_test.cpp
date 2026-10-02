#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/InteractionEnergyVrna.h"
#include "IntaRNA/OutputHandler.h"
#include "IntaRNA/PredictorSeedExtensionKinetic.h"
#include "IntaRNA/SeedHandlerExplicit.h"
#include "IntaRNA/VrnaHandler.h"

#include <algorithm>
#include <array>
#include <functional>
#include <map>
#include <string>
#include <tuple>
#include <vector>

using namespace IntaRNA;

namespace {

using Pair = std::pair<size_t, size_t>;
using Chain = std::vector<Pair>;
using Bounds = std::array<size_t, 4>;

class KineticOutput final : public OutputHandler {
public:
	explicit KineticOutput(const OutputConstraint & c);
	void add(const Interaction & i) override;
	std::vector<Interaction> interactions;
};

KineticOutput::KineticOutput(const OutputConstraint & c) : OutputHandler(c) {}

void KineticOutput::add(const Interaction & i) {
	if (!i.basePairs.empty()) interactions.push_back(i);
	++reportedInteractions;
}

// A nonmonotone table deliberately also models imported accessibility data.
class KineticAccessibility final : public AccessibilityDisabled {
public:
	KineticAccessibility(const RnaSequence & s, size_t maxLength = 0);
	E_type getED(size_t from, size_t to) const override;
	std::map<Pair, E_type> values;
};

KineticAccessibility::KineticAccessibility(const RnaSequence & s, size_t maxLength)
 : AccessibilityDisabled(s, maxLength, NULL) {}

E_type KineticAccessibility::getED(size_t from, size_t to) const {
	const E_type base = AccessibilityDisabled::getED(from, to);
	if (base == ED_UPPER_BOUND) return base;
	auto i = values.find({from, to});
	return i == values.end() ? 0 : i->second;
}

class KineticEnergy final : public InteractionEnergyBasePair {
public:
	KineticEnergy(const Accessibility & a, const ReverseAccessibility & b,
			size_t m1 = 3, size_t m2 = 3);
	E_type getE_interLeft(size_t i, size_t j, size_t k, size_t l) const override;
	E_type getE(size_t i, size_t j, size_t k, size_t l, E_type h) const override;
	bool customLoops = false;
	std::map<Bounds, E_type> loops;
	std::map<Bounds, E_type> boundaryTerms;
};

KineticEnergy::KineticEnergy(const Accessibility & a, const ReverseAccessibility & b,
		size_t m1, size_t m2)
 : InteractionEnergyBasePair(a, b, m1, m2, false, 1., -100, 3, 0, false) {}

E_type KineticEnergy::getE_interLeft(size_t i, size_t j, size_t k, size_t l) const {
	if (!isValidInternalLoop(i, j, k, l)) return E_INF;
	if (!customLoops) return InteractionEnergyBasePair::getE_interLeft(i, j, k, l);
	auto p = loops.find({i,j,k,l});
	return p == loops.end() ? E_INF : p->second;
}

E_type KineticEnergy::getE(size_t i, size_t j, size_t k, size_t l, E_type h) const {
	const E_type base = InteractionEnergyBasePair::getE(i,j,k,l,h);
	auto p = boundaryTerms.find({i,j,k,l});
	return E_isINF(base) ? E_INF : base + (p == boundaryTerms.end() ? 0 : p->second);
}

struct KineticFixture {
	RnaSequence first, second;
	KineticAccessibility acc1, acc2;
	ReverseAccessibility reversed;
	KineticEnergy energy;
	KineticFixture(size_t n = 7, size_t m1 = 3, size_t m2 = 3,
			size_t span1 = 0, size_t span2 = 0);
};

KineticFixture::KineticFixture(size_t n, size_t m1, size_t m2, size_t span1, size_t span2)
 : first("target", std::string(n, 'G')), second("query", std::string(n, 'C')),
	acc1(first,span1), acc2(second,span2), reversed(acc2), energy(acc1,reversed,m1,m2) {}

bool stacked(const Pair & a, const Pair & b) {
	return a.first+1 == b.first && a.second+1 == b.second;
}

Bounds bounds(const Chain & chain) {
	return {chain.front().first,chain.back().first,chain.front().second,chain.back().second};
}

// Independent reference: reconstruct every candidate's complete chain, check
// its topology and recompute all loop and boundary contributions from scratch.
E_type chainEnergy(const InteractionEnergy & energy, const Chain & chain,
		const OutputConstraint & out)
{
	const auto b = bounds(chain);
	if (b[1]-b[0]+1 > energy.getAccessibility1().getMaxLength()
			|| b[3]-b[2]+1 > energy.getAccessibility2().getMaxLength()) return E_INF;
	E_type h = energy.getE_init();
	for (size_t p = 0; p < chain.size(); ++p) {
		if (!energy.areComplementary(chain[p].first,chain[p].second)) return E_INF;
		const bool before = p != 0 && stacked(chain[p-1],chain[p]);
		const bool after = p+1 < chain.size() && stacked(chain[p],chain[p+1]);
		if (out.noLP && !before && !after) return E_INF;
		if (p == 0) continue;
		if (out.noGUend && !before && (energy.isGU(chain[p-1].first,chain[p-1].second)
				|| energy.isGU(chain[p].first,chain[p].second))) return E_INF;
		const E_type loop = energy.getE_interLeft(chain[p-1].first,chain[p].first,
				chain[p-1].second,chain[p].second);
		if (E_isINF(loop)) return E_INF;
		h += loop;
	}
	return energy.getE(b[0],b[1],b[2],b[3],h);
}

std::vector<Chain> oracle(const InteractionEnergy & energy, Chain chain,
		const OutputConstraint & out, char score,
		const IndexRange & r1, const IndexRange & r2)
{
	std::vector<Chain> result;
	if (E_isINF(chainEnergy(energy,chain,out))) return result;
	while (true) {
		const auto b = bounds(chain);
		const E_type current = chainEnergy(energy,chain,out);
		if ((!out.noGUend || (!energy.isGU(b[0],b[2]) && !energy.isGU(b[1],b[3])))
				&& energy.getED1(b[0],b[1]) <= out.maxED
				&& energy.getED2(b[2],b[3]) <= out.maxED) result.push_back(chain);
		bool found = false;
		Chain best;
		std::tuple<long double,int,size_t,size_t> bestKey;
		// Enumerate absolute endpoints in deliberately different order from
		// production's side/loop-size enumeration.
		for (size_t x = r1.from; x <= r1.to; ++x) {
			for (size_t y = r2.from; y <= r2.to; ++y) {
				const bool left = x < b[0] && y < b[2];
				const bool right = x > b[1] && y > b[3];
				if (!left && !right) continue;
				const size_t s1 = left ? b[0]-x-1 : x-b[1]-1;
				const size_t s2 = left ? b[2]-y-1 : y-b[3]-1;
				if (s1 > energy.getMaxInternalLoopSize1() || s2 > energy.getMaxInternalLoopSize2()) continue;
				Chain trial = chain;
				trial.push_back({x,y});
				if (out.noLP && s1+s2 != 0) {
					if (left && (x == r1.from || y == r2.from)) continue;
					if (right && (x == r1.to || y == r2.to)) continue;
					trial.push_back(left ? Pair{x-1,y-1} : Pair{x+1,y+1});
				}
				std::sort(trial.begin(),trial.end());
				const E_type total = chainEnergy(energy,trial,out);
				if (E_isINF(total) || total >= current) continue;
				const size_t denominator = score == 'A' ? 1 : score == 'B' ? 1+s1+s2 : 1+2*std::max(s1,s2);
				const auto key = std::make_tuple(static_cast<long double>(total-current)/denominator,
						left ? 0 : 1,s1+s2,s1);
				if (!found || key < bestKey) { found = true; bestKey = key; best = trial; }
			}
		}
		if (!found) return result;
		chain = best;
	}
}

std::string seedEncoding(const Chain & seed, size_t n2) {
	const auto b = bounds(seed);
	std::string a(b[1]-b[0]+1,'.'), bReverse(b[3]-b[2]+1,'.');
	for (const auto & p : seed) { a[p.first-b[0]] = '|'; bReverse[p.second-b[2]] = '|'; }
	std::reverse(bReverse.begin(),bReverse.end());
	return std::to_string(b[0]+1)+a+"&"+std::to_string(n2-b[3])+bReverse;
}

SeedConstraint seedConstraint(const std::string & explicitSeed) {
	return SeedConstraint(2,20,20,20,E_INF,Accessibility::ED_UPPER_BOUND,E_INF,
			IndexRangeList(""),IndexRangeList(""),explicitSeed,false,false,false);
}

std::vector<Interaction> predict(const InteractionEnergy & energy, const Chain & seed,
		const OutputConstraint & out, char score = 'A',
		const IndexRange & r1 = IndexRange(0,RnaSequence::lastPos),
		const IndexRange & r2 = IndexRange(0,RnaSequence::lastPos))
{
	const auto sc = seedConstraint(seedEncoding(seed,energy.size2()));
	KineticOutput output(out);
	PredictorSeedExtensionKinetic predictor(energy,output,NULL,new SeedHandlerExplicit(energy,sc),score);
	predictor.predict(r1,r2);
	return output.interactions;
}

Chain internalChain(const InteractionEnergy & energy, const Interaction & interaction) {
	Chain result;
	for (const auto & p : interaction.basePairs) result.push_back({energy.getIndex1(p),energy.getIndex2(p)});
	return result;
}

void checkOracle(const InteractionEnergy & energy, const Chain & seed,
		const OutputConstraint & out, char score,
		IndexRange r1 = IndexRange(0,RnaSequence::lastPos),
		IndexRange r2 = IndexRange(0,RnaSequence::lastPos))
{
	r1.to = std::min(r1.to,energy.size1()-1); r2.to = std::min(r2.to,energy.size2()-1);
	auto expected = oracle(energy,seed,out,score,r1,r2);
	expected.erase(std::remove_if(expected.begin(),expected.end(),[&](const Chain & c) {
		return chainEnergy(energy,c,out) >= out.maxE;
	}),expected.end());
	std::sort(expected.begin(),expected.end(),[&](const Chain & a,const Chain & b) {
		return chainEnergy(energy,a,out) < chainEnergy(energy,b,out);
	});
	const auto actual = predict(energy,seed,out,score,r1,r2);
	REQUIRE(actual.size() == std::min(expected.size(),out.reportMax));
	for (size_t p = 0; p < actual.size(); ++p) {
		REQUIRE(internalChain(energy,actual[p]) == expected[p]);
		REQUIRE(actual[p].energy == chainEnergy(energy,expected[p],out));
	}
}

} // namespace

TEST_CASE("Kinetic extension agrees with independent whole-chain oracle", "[PredictorSeedExtensionKinetic]") {
	#include "testEasyLoggingSetup.icc"
	const Chain seed{{2,2},{3,3}};
	for (char score : {'A','B','C'}) {
		for (bool noLP : {false,true}) {
			KineticFixture fixture(8,1,3,6,7);
			// Nonmonotone values distinguish exact candidate EDs from a
			// farthest-end penalty and change full boundary energy deltas.
			fixture.acc1.values = {{{0,3},700},{{1,3},30},{{2,5},500},{{2,6},20}};
			fixture.acc2.values = {{{2,5},450},{{1,5},20}};
			OutputConstraint out(100,OutputConstraint::OVERLAP_BOTH,E_INF,E_INF,false,noLP);
			checkOracle(fixture.energy,seed,out,score);
			checkOracle(fixture.energy,seed,out,score,IndexRange(1,6),IndexRange(0,7));
		}
	}
}

TEST_CASE("Kinetic scoring and strict downhill acceptance", "[PredictorSeedExtensionKinetic]") {
	#include "testEasyLoggingSetup.icc"
	const Chain seed{{2,2},{3,3}};
	OutputConstraint out(100,OutputConstraint::OVERLAP_BOTH,E_INF,E_INF);
	SECTION("A B and C select distinct justified moves") {
		KineticFixture f(8);
		f.energy.customLoops = true;
		f.energy.loops = {{{2,3,2,3},-100},{{3,6,3,4},-420},{{3,5,3,5},-390},{{3,4,3,4},-120}};
		const std::map<char,Pair> expected{{'A',{6,4}},{'B',{6,4}},{'C',{5,5}}};
		for (char score : {'A','B','C'}) {
			checkOracle(f.energy,seed,out,score);
			const auto actual = predict(f.energy,seed,out,score);
			REQUIRE(internalChain(f.energy,actual.front()).back() == expected.at(score));
		}
		// B prefers the immediate stack when its per-distance gain wins.
		f.energy.loops[{3,4,3,4}] = -150;
		REQUIRE(internalChain(f.energy,predict(f.energy,seed,out,'B').front()).back() == Pair(4,4));
	}
	SECTION("zero and positive moves stop even when the local loop is favorable") {
		for (E_type delta : {0,1,200}) {
			KineticFixture f(5);
			f.energy.customLoops = true;
			f.energy.loops = {{{2,3,2,3},-100},{{3,4,3,4},-100}};
			f.energy.boundaryTerms[{2,4,2,4}] = 100+delta;
			const auto actual = predict(f.energy,seed,out);
			REQUIRE(actual.size() == 1);
			REQUIRE(internalChain(f.energy,actual.front()) == seed);
		}
	}
	SECTION("ties prefer left then smaller total gap then smaller first gap") {
		KineticFixture f(7);
		f.energy.customLoops = true;
		f.energy.loops = {{{2,3,2,3},-100},{{0,2,1,2},-100},{{1,2,0,2},-100},{{1,2,1,2},-100},{{3,4,3,4},-100}};
		// Limit each span to three: whichever move wins blocks the other side.
		KineticFixture shortF(7,3,3,3,3);
		shortF.energy.customLoops = true; shortF.energy.loops = f.energy.loops;
		REQUIRE(internalChain(shortF.energy,predict(shortF.energy,seed,out).front()).front() == Pair(1,1));
		f.energy.loops.erase({1,2,1,2});
		// Equal total gap: s1=0 (pair 1,0) wins over s1=1 (pair 0,1).
		const auto paths = oracle(f.energy,seed,out,'A',IndexRange(0,6),IndexRange(0,6));
		REQUIRE(paths.at(1).front() == Pair(1,0));
		checkOracle(f.energy,seed,out,'A');
	}
}

TEST_CASE("Kinetic noLP macro-steps and explicit seed validation", "[PredictorSeedExtensionKinetic]") {
	#include "testEasyLoggingSetup.icc"
	KineticFixture f(6,1,1);
	f.energy.customLoops = true;
	f.energy.loops = {{{0,1,0,1},-100},{{1,3,1,3},50},{{3,4,3,4},-330}};
	OutputConstraint out(100,OutputConstraint::OVERLAP_BOTH,E_INF,E_INF,false,true);
	const Chain seed{{0,0},{1,1}};
	checkOracle(f.energy,seed,out,'A');
	const auto actual = predict(f.energy,seed,out);
	REQUIRE(actual.size() == 2);
	REQUIRE(actual.front().basePairs.size() == 4);
	REQUIRE(actual.front().energy == -480);
	SECTION("B and C count one macro move in their denominator") {
		KineticFixture ranked(7,1,1,5,5);
		ranked.energy.customLoops = true;
		ranked.energy.loops = {{{2,3,2,3},-100},{{1,2,1,2},-80},
				{{3,5,3,5},50},{{5,6,5,6},-330}};
		const Chain middle{{2,2},{3,3}};
		for (char score : {'B','C'}) {
			checkOracle(ranked.energy,middle,out,score);
			const auto rankedResult = predict(ranked.energy,middle,out,score);
			// -280/3 beats -80. Incorrectly counting the two added pairs
			// would instead produce -280/4 and choose the left stack.
			REQUIRE(rankedResult.front().basePairs.size() == 4);
			REQUIRE(internalChain(ranked.energy,rankedResult.front()).back() == Pair(6,6));
		}
	}
	SECTION("lonely explicit seeds are skipped") {
		REQUIRE(predict(f.energy,Chain{{1,1},{3,3}},out).empty());
	}
	SECTION("both additional pairs must fit both strand spans") {
		KineticFixture shortF(6,1,1,4,6);
		shortF.energy.customLoops = true; shortF.energy.loops = f.energy.loops;
		const auto limited = predict(shortF.energy,seed,out);
		REQUIRE(limited.size() == 1);
		REQUIRE(limited.front().basePairs.size() == 2);
	}
}

TEST_CASE("Kinetic noGU output retains a valid prefix", "[PredictorSeedExtensionKinetic]") {
	#include "testEasyLoggingSetup.icc"
	RnaSequence first("target","GGG"), second("query","UCC");
	AccessibilityDisabled a(first,0,NULL), b(second,0,NULL);
	ReverseAccessibility reversed(b);
	InteractionEnergyBasePair energy(a,reversed,1,1);
	OutputConstraint out(100,OutputConstraint::OVERLAP_BOTH,E_INF,E_INF,false,false,true);
	const Chain seed{{0,0},{1,1}};
	checkOracle(energy,seed,out,'A');
	const auto actual = predict(energy,seed,out);
	REQUIRE(actual.size() == 1);
	REQUIRE(internalChain(energy,actual.front()) == seed);
	SECTION("a transient GU endpoint can become internal after another stack") {
		RnaSequence t("target","GGGG"), q("query","CUCC");
		AccessibilityDisabled at(t,0,NULL), aq(q,0,NULL);
		ReverseAccessibility rq(aq);
		InteractionEnergyBasePair model(at,rq,1,1);
		checkOracle(model,seed,out,'A');
		const auto recovered = predict(model,seed,out);
		REQUIRE(recovered.size() == 2);
		REQUIRE(recovered.front().basePairs.size() == 4);
		REQUIRE(recovered.front().energy == -400);
		REQUIRE(model.isGU(2,2));
	}
}

TEST_CASE("Kinetic full ViennaRNA energy agrees with whole-chain oracle", "[PredictorSeedExtensionKinetic]") {
	#include "testEasyLoggingSetup.icc"
	RnaSequence first("target","AGCGACGCA"), second("query","UGCGUCGCU");
	KineticAccessibility a(first), b(second);
	a.values = {{{0,4},140},{{1,4},20},{{1,5},350},{{1,6},30},{{2,7},130}};
	b.values = {{{0,4},20},{{1,4},60},{{1,5},20},{{2,6},90}};
	ReverseAccessibility reversed(b);
	VrnaHandler vrna(37,"Turner04",false,false);
	const Chain seed{{2,2},{3,3}};
	for (bool dangles : {false,true}) {
		InteractionEnergyVrna energy(a,reversed,vrna,3,2,false,37,dangles);
		OutputConstraint out(100,OutputConstraint::OVERLAP_BOTH,E_INF,E_INF);
		for (char score : {'A','B','C'}) checkOracle(energy,seed,out,score);
	}
}

TEST_CASE("Kinetic refuses undefined ensemble output and invalid scores", "[PredictorSeedExtensionKinetic]") {
	#include "testEasyLoggingSetup.icc"
	KineticFixture f;
	const auto sc = seedConstraint("3||&4||");
	OutputConstraint ensemble(1,OutputConstraint::OVERLAP_BOTH,0,E_INF,false,false,false,true);
	KineticOutput ensembleOut(ensemble);
	REQUIRE_THROWS_AS(PredictorSeedExtensionKinetic(f.energy,ensembleOut,NULL,
			new SeedHandlerExplicit(f.energy,sc)),std::invalid_argument);
	OutputConstraint ordinary;
	KineticOutput ordinaryOut(ordinary);
	REQUIRE_THROWS_AS(PredictorSeedExtensionKinetic(f.energy,ordinaryOut,NULL,
			new SeedHandlerExplicit(f.energy,sc),'Z'),std::invalid_argument);
}

TEST_CASE("Kinetic Turner loop is rescued by its atomic stack", "[PredictorSeedExtensionKinetic]") {
	#include "testEasyLoggingSetup.icc"
	// A single query-strand A bulge prevents a direct stack after the seed.
	// This is the Turner2004 witness reproduced by the first council:
	// +0.50 kcal/mol loop plus -3.30 kcal/mol outward stack.
	RnaSequence first("target","CCCC"), second("query","GGAGG");
	AccessibilityDisabled a(first,0,NULL), b(second,0,NULL);
	ReverseAccessibility reversed(b);
	VrnaHandler vrna(37,"Turner04",false,false);
	InteractionEnergyVrna energy(a,reversed,vrna,1,1,false,0,false);
	const E_type loop = energy.getE_interLeft(1,2,1,3);
	const E_type stack = energy.getE_interLeft(2,3,3,4);
	REQUIRE(loop == 50);
	REQUIRE(stack == -330);
	OutputConstraint out(100,OutputConstraint::OVERLAP_BOTH,E_INF,E_INF,false,true);
	const Chain seed{{0,0},{1,1}};
	checkOracle(energy,seed,out,'A');
	const auto actual = predict(energy,seed,out);
	REQUIRE(actual.size() == 2);
	REQUIRE((internalChain(energy,actual.front()) == Chain{{0,0},{1,1},{2,3},{3,4}}));
	REQUIRE(actual.front().energy - actual.back().energy == loop+stack);
}

TEST_CASE("Kinetic nonoverlap reporting recovers shorter prefixes", "[PredictorSeedExtensionKinetic]") {
	#include "testEasyLoggingSetup.icc"
	KineticFixture f(8,0,0,5,5);
	f.energy.customLoops = true;
	for (size_t p = 0; p+1 < 8; ++p) f.energy.loops[{p,p+1,p,p+1}] = p < 3 ? -100 : -200;
	const auto sc = seedConstraint("1||&7||,4||&4||");
	for (bool needBPs : {false,true}) {
		OutputConstraint out(3,OutputConstraint::OVERLAP_NONE,E_INF,E_INF,false,false,false,false,needBPs);
		KineticOutput output(out);
		PredictorSeedExtensionKinetic predictor(f.energy,output,NULL,new SeedHandlerExplicit(f.energy,sc));
		predictor.predict();
		REQUIRE(output.interactions.size() == 2);
		REQUIRE(output.interactions[0].energy == -900);
		REQUIRE(output.interactions[1].energy == -300);
		REQUIRE((bounds(internalChain(f.energy,output.interactions[0])) == Bounds{{3,7,3,7}}));
		REQUIRE((bounds(internalChain(f.energy,output.interactions[1])) == Bounds{{0,2,0,2}}));
		if (needBPs) {
			REQUIRE(output.interactions[0].basePairs.size() == 5);
			REQUIRE(output.interactions[1].basePairs.size() == 3);
		}
		// Reuse the instance with a different range and ensure old candidates
		// and seed metadata cannot leak into the second prediction.
		output.interactions.clear();
		predictor.predict(IndexRange(0,2),IndexRange(0,2));
		REQUIRE(output.interactions.size() == 1);
		REQUIRE(output.interactions.front().energy == -300);
		REQUIRE((bounds(internalChain(f.energy,output.interactions.front())) == Bounds{{0,2,0,2}}));
	}
}

TEST_CASE("Kinetic annotations exclude rejected lonely explicit seeds", "[PredictorSeedExtensionKinetic]") {
	#include "testEasyLoggingSetup.icc"
	KineticFixture f(5,0,0);
	const auto sc = seedConstraint("2||&3||,1|&5|");
	OutputConstraint out(1,OutputConstraint::OVERLAP_BOTH,E_INF,E_INF,false,true);
	KineticOutput output(out);
	PredictorSeedExtensionKinetic predictor(f.energy,output,NULL,new SeedHandlerExplicit(f.energy,sc));
	predictor.predict();
	REQUIRE(output.interactions.size() == 1);
	const auto & result = output.interactions.front();
	REQUIRE(result.basePairs.size() == 5);
	REQUIRE(result.seed != NULL);
	REQUIRE(result.seed->size() == 1);
}
