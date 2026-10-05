#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/InteractionEnergyVrna.h"
#include "IntaRNA/OutputHandler.h"
#include "IntaRNA/PredictorSeedExtensionKinetic.h"
#include "IntaRNA/PredictorSeedExtensionKineticPruned.h"
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
	mutable std::map<Bounds, size_t> loopCalls;
};

KineticEnergy::KineticEnergy(const Accessibility & a, const ReverseAccessibility & b,
		size_t m1, size_t m2)
 : InteractionEnergyBasePair(a, b, m1, m2, false, 1., -100, 3, 0, false) {}

E_type KineticEnergy::getE_interLeft(size_t i, size_t j, size_t k, size_t l) const {
	++loopCalls[{i,j,k,l}];
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
		if (p == 0) continue;
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
		std::tuple<long double,int,size_t,size_t,bool> bestKey;
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
				if (s1+s2 != 0 && (out.noGUend || !energy.isInternalLoopGUallowed())
						&& (energy.isGU(left ? b[0] : b[1], left ? b[2] : b[3]) || energy.isGU(x,y))) continue;
				for (bool macro : {false,true}) {
					if (!macro && s1+s2 != 0) continue;
					Chain trial = chain;
					trial.push_back({x,y});
					if (macro) {
						if (left && (x == r1.from || y == r2.from)) continue;
						if (right && (x == r1.to || y == r2.to)) continue;
						trial.push_back(left ? Pair{x-1,y-1} : Pair{x+1,y+1});
					}
					std::sort(trial.begin(),trial.end());
					const E_type total = chainEnergy(energy,trial,out);
					if (E_isINF(total) || total >= current) continue;
					const size_t denominator = score == 'A' ? 1 : score == 'B' ? 1+s1+s2 : 1+2*std::max(s1,s2);
					const auto key = std::make_tuple(static_cast<long double>(total-current)/denominator,
							left ? 0 : 1,s1+s2,s1,macro);
					if (!found || key < bestKey) { found = true; bestKey = key; best = trial; }
				}
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
		f.energy.loops = {{{2,3,2,3},-100},{{3,6,3,4},-420},{{3,5,3,5},-390},{{3,4,3,4},-120},{{6,7,4,5},0},{{5,6,5,6},0}};
		const std::map<char,Pair> expected{{'A',{7,5}},{'B',{7,5}},{'C',{6,6}}};
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
		f.energy.loops = {{{3,4,3,4},-100},{{2,3,1,3},-50},{{1,2,0,1},-50},
				{{1,3,2,3},-50},{{0,1,1,2},-50}};
		const Chain middle{{3,3},{4,4}};
		// Equal total gap: s1=0 (outer pair 1,0) wins over s1=1 (0,1).
		const auto paths = oracle(f.energy,middle,out,'A',IndexRange(0,6),IndexRange(0,6));
		REQUIRE(paths.at(1).front() == Pair(1,0));
		checkOracle(f.energy,middle,out,'A');
	}
}

TEST_CASE("Kinetic always uses stacked extensions and trusts explicit seeds", "[PredictorSeedExtensionKinetic]") {
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
	SECTION("lonely explicit seeds are accepted") {
		REQUIRE_FALSE(predict(f.energy,Chain{{1,1},{3,3}},out).empty());
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
		REQUIRE(output.interactions[1].energy == -200);
		REQUIRE((bounds(internalChain(f.energy,output.interactions[0])) == Bounds{{3,7,3,7}}));
		REQUIRE((bounds(internalChain(f.energy,output.interactions[1])) == Bounds{{0,1,0,1}}));
		if (needBPs) {
			REQUIRE(output.interactions[0].basePairs.size() == 5);
			REQUIRE(output.interactions[1].basePairs.size() == 2);
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

TEST_CASE("Kinetic annotations retain handler-provided lonely explicit seeds", "[PredictorSeedExtensionKinetic]") {
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
	REQUIRE(result.seed->size() == 2);
}

TEST_CASE("Kinetic two-stack moves cross a local barrier atomically", "[PredictorSeedExtensionKinetic]") {
	#include "testEasyLoggingSetup.icc"
	KineticFixture f(4,0,0);
	f.energy.customLoops = true;
	f.energy.loops = {{{0,1,0,1},-100},{{1,2,1,2},50},{{2,3,2,3},-200}};
	const Chain seed{{0,0},{1,1}};
	for (bool noLP : {false,true}) {
		OutputConstraint out(100,OutputConstraint::OVERLAP_BOTH,E_INF,E_INF,false,noLP);
		checkOracle(f.energy,seed,out,'A');
		const auto actual = predict(f.energy,seed,out);
		REQUIRE(actual.size() == 2);
		REQUIRE(actual.front().basePairs.size() == 4);
		REQUIRE(actual.front().energy == -350);
	}
}

TEST_CASE("Kinetic exact move ties prefer the single stack", "[PredictorSeedExtensionKinetic]") {
	#include "testEasyLoggingSetup.icc"
	KineticFixture f(4,0,0);
	f.energy.customLoops = true;
	f.energy.loops = {{{0,1,0,1},-100},{{1,2,1,2},-100},{{2,3,2,3},0}};
	OutputConstraint out(100,OutputConstraint::OVERLAP_BOTH,E_INF,E_INF);
	for (char score : {'A','B','C'}) {
		const auto actual = predict(f.energy,Chain{{0,0},{1,1}},out,score);
		REQUIRE(actual.size() == 2);
		REQUIRE(actual.front().basePairs.size() == 3);
		REQUIRE(actual.front().energy == -300);
	}
}

TEST_CASE("Kinetic caches the unchanged end and trusts the seed energy", "[PredictorSeedExtensionKinetic]") {
	#include "testEasyLoggingSetup.icc"
	KineticFixture f(8,0,0);
	f.energy.customLoops = true;
	for (size_t i = 0; i+1 < 8; ++i) f.energy.loops[{i,i+1,i,i+1}] = i < 2 ? -200 : -100;
	OutputConstraint out(100,OutputConstraint::OVERLAP_BOTH,E_INF,E_INF);
	const auto result = predict(f.energy,Chain{{2,2},{3,3}},out);
	REQUIRE(result.front().basePairs.size() == 8);
	// Seed energy is evaluated by SeedHandlerExplicit just once, never by
	// the predictor. The right first stack is shared by single/two-pair moves
	// and is not reevaluated after the left end grows.
	REQUIRE(f.energy.loopCalls.at(Bounds{{2,3,2,3}}) == 1);
	REQUIRE(f.energy.loopCalls.at(Bounds{{3,4,3,4}}) == 2);
}

TEST_CASE("Kinetic oracle covers varied sequences and endpoint energies", "[PredictorSeedExtensionKinetic]") {
	#include "testEasyLoggingSetup.icc"
	unsigned int random = 254;
	const auto next = [&]() { random = random*1664525u+1013904223u; return random; };
	for (size_t sample = 0; sample < 32; ++sample) {
		std::string t(9,'G'), q(9,'C');
		for (size_t i = 0; i < 9; ++i) { t[i] = "ACGU"[(next() >> 16)%4]; q[i] = "ACGU"[(next() >> 16)%4]; }
		t[3] = t[4] = 'G'; q[4] = q[5] = 'C';
		RnaSequence first("t",t), second("q",q);
		KineticAccessibility a(first), b(second);
		for (size_t i = 0; i < 9; ++i) for (size_t j = i; j < 9; ++j) {
			a.values[{i,j}] = (next() >> 16)%300;
			b.values[{i,j}] = (next() >> 16)%300;
		}
		ReverseAccessibility reversed(b);
		InteractionEnergyBasePair energy(a,reversed,2,3);
		OutputConstraint out(100,OutputConstraint::OVERLAP_BOTH,E_INF,E_INF,false,sample%2,sample%3 == 0);
		for (char score : {'A','B','C'}) checkOracle(energy,Chain{{3,3},{4,4}},out,score);
	}
}

TEST_CASE("Pruned kinetic mode preserves complete energies and resets offsets", "[PredictorSeedExtensionKinetic]") {
	#include "testEasyLoggingSetup.icc"
	RnaSequence t("t","GGGGGGGG"), q("q","CCCCCCCC");
	KineticAccessibility a(t), b(q);
	for (size_t i = 0; i < 8; ++i) for (size_t j = i; j < 8; ++j) {
		a.values[{i,j}] = 20*(j-i); b.values[{i,j}] = 25*(j-i);
	}
	ReverseAccessibility reversed(b);
	VrnaHandler vrna(25,"Turner99",false,false);
	InteractionEnergyVrna energy(a,reversed,vrna,2,3,false,37,false);
	const Chain seed{{3,3},{4,4}};
	const auto sc = seedConstraint(seedEncoding(seed,8));
	for (char score : {'A','B','C'}) {
		OutputConstraint out(100,OutputConstraint::OVERLAP_BOTH,E_INF,E_INF);
		KineticOutput output(out);
		PredictorSeedExtensionKineticPruned predictor(energy,output,NULL,new SeedHandlerExplicit(energy,sc),score);
		for (const IndexRange & range : {IndexRange(0,7),IndexRange(1,6),IndexRange(2,7),IndexRange(0,7)}) {
			output.interactions.clear(); predictor.predict(range,range);
			const auto expected = predict(energy,seed,out,score,range,range);
			REQUIRE(output.interactions.size() == expected.size());
			for (size_t i = 0; i < expected.size(); ++i) {
				REQUIRE(output.interactions[i].basePairs == expected[i].basePairs);
				REQUIRE(output.interactions[i].energy == chainEnergy(energy,internalChain(energy,expected[i]),out));
			}
		}
	}
}

TEST_CASE("Pruned kinetic mode documents its nonmonotone ED tradeoff", "[PredictorSeedExtensionKinetic]") {
	#include "testEasyLoggingSetup.icc"
	KineticFixture f(6,2,2);
	f.acc1.values = {{{0,2},200},{{0,3},300}};
	InteractionEnergyBasePair energy(f.acc1,f.reversed,2,2);
	const Chain seed{{0,0},{1,1}};
	const auto sc = seedConstraint(seedEncoding(seed,6));
	OutputConstraint out(100,OutputConstraint::OVERLAP_BOTH,E_INF,E_INF);
	KineticOutput output(out);
	PredictorSeedExtensionKineticPruned predictor(energy,output,NULL,new SeedHandlerExplicit(energy,sc));
	predictor.predict();
	REQUIRE(output.interactions.size() == 1);
	REQUIRE(output.interactions.front().basePairs.size() == 2);
	REQUIRE(predict(energy,seed,out).front().basePairs.size() > 2);
	// The subclass has no bound for custom energy overrides: safely use K.
	KineticOutput customOutput(out);
	PredictorSeedExtensionKineticPruned custom(f.energy,customOutput,NULL,new SeedHandlerExplicit(f.energy,sc));
	custom.predict();
	REQUIRE(customOutput.interactions.front().basePairs == predict(f.energy,seed,out).front().basePairs);
}

namespace {

class KineticPruningProbe : public PredictorSeedExtensionKineticPruned {
public:
	using PredictorSeedExtensionKineticPruned::PredictorSeedExtensionKineticPruned;
	bool skips(const Bounds & before, const Bounds & after, bool left, size_t s1, size_t s2) const;
};

bool KineticPruningProbe::skips(const Bounds & before, const Bounds & after,
		bool left, size_t s1, size_t s2) const
{
	Candidate c;
	c.bounds = after; c.left = left; c.s1 = s1; c.s2 = s2; c.macro = true;
	return prune(c,before);
}

} // namespace

TEST_CASE("Pruning tables bound local moves at all oriented root types", "[PredictorSeedExtensionKinetic]") {
	#include "testEasyLoggingSetup.icc"
	RnaSequence t("t","CGGUAUCGGUAU");
	std::string reversedQuery = "GCUGUAGCUGUA";
	std::reverse(reversedQuery.begin(),reversedQuery.end());
	RnaSequence q("q",reversedQuery);
	KineticAccessibility a(t), b(q);
	ReverseAccessibility reversed(b);
	OutputConstraint out(100,OutputConstraint::OVERLAP_BOTH,E_INF,E_INF);
	KineticOutput output(out);
	size_t checked = 0;
	for (const std::string & parameters : {"Turner04","Turner99"}) {
		for (double temperature : {20.,37.}) {
			VrnaHandler vrna(temperature,parameters,false,false);
			InteractionEnergyVrna energy(a,reversed,vrna,2,3,false,0,false);
			const auto sc = seedConstraint(seedEncoding(Chain{{0,0},{1,1}},12));
			KineticPruningProbe probe(energy,output,NULL,new SeedHandlerExplicit(energy,sc));
			for (bool left : {false,true}) for (size_t root = 0; root < 12; ++root)
			for (size_t s1 = 0; s1 <= 2; ++s1) for (size_t s2 = 0; s2 <= 3; ++s2) {
				if ((left && root < std::max(s1,s2)+2)
						|| (!left && root+std::max(s1,s2)+2 >= 12)) continue;
				const size_t c1 = left ? root-s1-1 : root+s1+1;
				const size_t c2 = left ? root-s2-1 : root+s2+1;
				const size_t o1 = left ? c1-1 : c1+1, o2 = left ? c2-1 : c2+1;
				const E_type loop = left ? energy.getE_interLeft(c1,root,c2,root)
						: energy.getE_interLeft(root,c1,root,c2);
				const E_type stack = left ? energy.getE_interLeft(o1,c1,o2,c2)
						: energy.getE_interLeft(c1,o1,c2,o2);
				if (E_isINF(loop) || E_isINF(stack) || loop+stack >= 0) continue;
				const Bounds before{root,root,root,root};
				const Bounds after = left ? Bounds{o1,root,o2,root} : Bounds{root,o1,root,o2};
				a.values.clear();
				// A real local move with delta=-1 must never be rejected by
				// a valid local suffix bound (no endpoint terms involved).
				a.values[{after[0],after[1]}] = -(loop+stack)-1;
				REQUIRE_FALSE(probe.skips(before,after,left,s1,s2));
				++checked;
			}
			a.values.clear();
		}
	}
	REQUIRE(checked > 50);
}
