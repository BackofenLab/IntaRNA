#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/OutputHandler.h"
#include "IntaRNA/OutputHandlerInteractionList.h"
#include "IntaRNA/PredictorMfe2dSeed.h"
#include "IntaRNA/PredictorMfe2dSeedExtension.h"
#include "IntaRNA/PredictorMfeEns2dSeedExtension.h"
#include "IntaRNA/PredictorMfeEnsSeedOnly.h"
#include "IntaRNA/PredictorMfeSeedOnly.h"
#include "IntaRNA/ReverseAccessibility.h"
#include "IntaRNA/RnaSequence.h"
#include "IntaRNA/SeedConstraint.h"
#include "IntaRNA/SeedHandlerMfe.h"

#include <algorithm>
#include <array>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace IntaRNA;

namespace {

struct SeedOraclePair {
	size_t i1;
	size_t i2;

	bool operator==(const SeedOraclePair &) const = default;
};

using SeedOracleChain = std::vector<SeedOraclePair>;

struct SeedOccurrence {
	size_t firstPair;
	size_t lastPair;
	E_type hybridEnergy;
	E_type totalEnergy;
};

struct SeededOracleInteraction {
	SeedOracleChain chain;
	std::vector<SeedOccurrence> seeds;
	E_type energy;
};

struct SeedOnlyOracleInteraction {
	SeedOracleChain chain;
	E_type energy;
};

struct SeedOracleResult {
	std::vector<SeededOracleInteraction> interactions;
	std::vector<SeedOnlyOracleInteraction> seedOnlyInteractions;
	E_type mfe = E_INF;
	Z_type unionPartition = 0.0;
};

struct SeedOnlySummary {
	std::vector<SeedOnlyOracleInteraction> optimalInteractions;
	E_type mfe = E_INF;
	Z_type partition = 0.0;
};

class SeedRecordingOutput final : public OutputHandler {
public:
	explicit SeedRecordingOutput(const OutputConstraint & constraint)
	 : OutputHandler(constraint)
	{}

	void add(const Interaction & interaction) override {
		calls.push_back(std::make_unique<Interaction>(interaction));
		++reportedInteractions;
	}

	const Interaction & last() const {
		return *calls.back();
	}

	std::vector<std::unique_ptr<Interaction> > calls;
};

size_t seedOracleLast(const IndexRange & range, const size_t sequenceSize) {
	return range.to == RnaSequence::lastPos ? sequenceSize-1 : range.to;
}

bool isStacked(const SeedOraclePair & left, const SeedOraclePair & right) {
	return left.i1+1 == right.i1 && left.i2+1 == right.i2;
}

bool hasNoLonelyPairs(const SeedOracleChain & chain, const size_t first,
		const size_t last)
{
	for (size_t i = first; i <= last; ++i) {
		const bool stackedLeft = i > first && isStacked(chain[i-1], chain[i]);
		const bool stackedRight = i < last && isStacked(chain[i], chain[i+1]);
		if (!stackedLeft && !stackedRight) {
			return false;
		}
	}
	return first <= last;
}

bool isAllowedSeedPair(const InteractionEnergy & energy,
		const SeedConstraint & constraint, const SeedOraclePair & pair,
		const bool atSeedEnd)
{
	return energy.areComplementary(pair.i1, pair.i2)
			&& energy.getED1(pair.i1, pair.i1) <= constraint.getMaxED()
			&& energy.getED2(pair.i2, pair.i2) <= constraint.getMaxED()
			&& (constraint.isGUallowed() || !energy.isGU(pair.i1, pair.i2))
			&& (!atSeedEnd || constraint.isGUendAllowed()
					|| !energy.isGU(pair.i1, pair.i2))
			&& (constraint.getRanges1().empty()
					|| constraint.getRanges1().covers(pair.i1))
			&& (constraint.getRanges2().empty()
					|| constraint.getRanges2().covers(pair.i2));
}

bool evaluateSeed(const InteractionEnergy & energy,
		const SeedConstraint & constraint, const SeedOracleChain & chain,
		const size_t first, SeedOccurrence & seed)
{
	const size_t seedBP = constraint.getBasePairs();
	if (first+seedBP > chain.size()) {
		return false;
	}
	const size_t last = first+seedBP-1;
	const SeedOraclePair & left = chain[first];
	const SeedOraclePair & right = chain[last];
	const size_t unpaired1 = right.i1-left.i1+1-seedBP;
	const size_t unpaired2 = right.i2-left.i2+1-seedBP;
	if (unpaired1 > constraint.getMaxUnpaired1()
			|| unpaired2 > constraint.getMaxUnpaired2()
			|| unpaired1+unpaired2 > constraint.getMaxUnpairedOverall())
	{
		return false;
	}
	if (!constraint.isLpAllowed() && !hasNoLonelyPairs(chain, first, last)) {
		return false;
	}
	for (size_t i = first; i <= last; ++i) {
		if (!isAllowedSeedPair(energy, constraint, chain[i],
				i == first || i == last))
		{
			return false;
		}
	}

	E_type hybridEnergy = energy.getE_init();
	for (size_t i = first+1; i <= last; ++i) {
		const E_type loopEnergy = energy.getE_interLeft(chain[i-1].i1,
				chain[i].i1, chain[i-1].i2, chain[i].i2);
		if (E_isINF(loopEnergy)) {
			return false;
		}
		hybridEnergy += loopEnergy;
	}
	if (hybridEnergy > constraint.getMaxEhybrid()) {
		return false;
	}
	// SeedHandlerMfe treats the interval ED limit as exclusive.
	if (energy.getED1(left.i1, right.i1) >= constraint.getMaxED()
			|| energy.getED2(left.i2, right.i2) >= constraint.getMaxED())
	{
		return false;
	}
	const E_type totalEnergy = energy.getE(left.i1, right.i1,
			left.i2, right.i2, hybridEnergy);
	if (E_isINF(totalEnergy) || totalEnergy > constraint.getMaxE()) {
		return false;
	}

	seed = SeedOccurrence{first, last, hybridEnergy, totalEnergy};
	return true;
}

std::vector<SeedOccurrence> findSeeds(const InteractionEnergy & energy,
		const SeedConstraint & constraint, const SeedOracleChain & chain)
{
	std::vector<SeedOccurrence> seeds;
	if (constraint.getExplicitSeeds().size() != 0) {
		throw std::logic_error(
				"the computed-seed oracle does not interpret explicit seed strings");
	}
	for (size_t first = 0; first+constraint.getBasePairs() <= chain.size();
			++first)
	{
		SeedOccurrence seed{};
		if (evaluateSeed(energy, constraint, chain, first, seed)) {
			seeds.push_back(seed);
		}
	}
	return seeds;
}

bool passesOutputConstraint(const InteractionEnergy & energy,
		const OutputConstraint & constraint, const SeedOracleChain & chain)
{
	if (constraint.noLP && !hasNoLonelyPairs(chain, 0, chain.size()-1)) {
		return false;
	}
	const SeedOraclePair & left = chain.front();
	const SeedOraclePair & right = chain.back();
	if (right.i1-left.i1+1>energy.getAccessibility1().getMaxLength()
			|| right.i2-left.i2+1>energy.getAccessibility2().getMaxLength()) return false;
	if (constraint.noGUend
			&& (energy.isGU(left.i1, left.i2)
					|| energy.isGU(right.i1, right.i2)))
	{
		return false;
	}
	return energy.getED1(left.i1, right.i1) <= constraint.maxED
			&& energy.getED2(left.i2, right.i2) <= constraint.maxED;
}

void recordSeedOracleChain(const InteractionEnergy & energy,
		const SeedConstraint & seedConstraint,
		const OutputConstraint & outputConstraint,
		const SeedOracleChain & chain, const E_type hybridEnergy,
		SeedOracleResult & result)
{
	const std::vector<SeedOccurrence> seeds = findSeeds(energy,
			seedConstraint, chain);
	if (chain.size() == seedConstraint.getBasePairs() && seeds.size() == 1
			&& passesOutputConstraint(energy, outputConstraint, chain))
	{
		result.seedOnlyInteractions.push_back(
				SeedOnlyOracleInteraction{chain, seeds.front().totalEnergy});
	}
	if (seeds.empty() || !passesOutputConstraint(energy, outputConstraint, chain)) {
		return;
	}
	const SeedOraclePair & left = chain.front();
	const SeedOraclePair & right = chain.back();
	const E_type totalEnergy = energy.getE(left.i1, right.i1,
			left.i2, right.i2, hybridEnergy);
	if (E_isINF(totalEnergy)) {
		return;
	}
	result.interactions.push_back(SeededOracleInteraction{
			chain, seeds, totalEnergy});
	result.mfe = std::min(result.mfe, totalEnergy);
	// A chain belongs to the seeded ensemble once, irrespective of how many
	// overlapping seed witnesses it contains.
	result.unionPartition += energy.getBoltzmannWeight(totalEnergy);
}

void enumerateSeedOracleExtensions(const InteractionEnergy & energy,
		const SeedConstraint & seedConstraint,
		const OutputConstraint & outputConstraint, const size_t to1,
		const size_t to2, SeedOracleChain & chain, const E_type hybridEnergy,
		SeedOracleResult & result)
{
	recordSeedOracleChain(energy, seedConstraint, outputConstraint, chain,
			hybridEnergy, result);
	const SeedOraclePair previous = chain.back();
	for (size_t next1 = previous.i1+1; next1 <= to1; ++next1) {
		for (size_t next2 = previous.i2+1; next2 <= to2; ++next2) {
			if (!energy.areComplementary(next1, next2)) {
				continue;
			}
			const E_type loopEnergy = energy.getE_interLeft(previous.i1,
					next1, previous.i2, next2);
			if (E_isINF(loopEnergy)) {
				continue;
			}
			chain.push_back(SeedOraclePair{next1, next2});
			enumerateSeedOracleExtensions(energy, seedConstraint,
					outputConstraint, to1, to2, chain,
					hybridEnergy+loopEnergy, result);
			chain.pop_back();
		}
	}
}

/**
 * Exhaustive reference calculation for tiny inputs.  It enumerates physical
 * base-pair chains directly and checks their seed witnesses afterwards; it
 * has no predictor matrix, seed-handler state, or dynamic-programming
 * recurrence.
 */
SeedOracleResult enumerateSeededInteractions(const InteractionEnergy & energy,
		const SeedConstraint & seedConstraint,
		const OutputConstraint & outputConstraint, const IndexRange & range1,
		const IndexRange & range2)
{
	SeedOracleResult result;
	const size_t to1 = seedOracleLast(range1, energy.size1());
	const size_t to2 = seedOracleLast(range2, energy.size2());
	for (size_t i1 = range1.from; i1 <= to1; ++i1) {
		for (size_t i2 = range2.from; i2 <= to2; ++i2) {
			if (!energy.areComplementary(i1, i2)) {
				continue;
			}
			SeedOracleChain chain(1, SeedOraclePair{i1, i2});
			enumerateSeedOracleExtensions(energy, seedConstraint,
					outputConstraint, to1, to2, chain, energy.getE_init(),
					result);
		}
	}
	return result;
}

SeedOnlySummary summarizeMfeSeedPerStart(const InteractionEnergy & energy,
		const SeedOracleResult & oracle)
{
	using Start = std::pair<size_t, size_t>;
	std::map<Start, E_type> bestByStart;
	for (const SeedOnlyOracleInteraction & seed : oracle.seedOnlyInteractions) {
		const Start start(seed.chain.front().i1, seed.chain.front().i2);
		auto entry = bestByStart.find(start);
		if (entry == bestByStart.end() || seed.energy < entry->second) {
			bestByStart[start] = seed.energy;
		}
	}

	SeedOnlySummary summary;
	for (const auto & entry : bestByStart) {
		summary.partition += energy.getBoltzmannWeight(entry.second);
		summary.mfe = std::min(summary.mfe, entry.second);
	}
	for (const SeedOnlyOracleInteraction & seed : oracle.seedOnlyInteractions) {
		const Start start(seed.chain.front().i1, seed.chain.front().i2);
		if (seed.energy == bestByStart.at(start) && seed.energy == summary.mfe) {
			summary.optimalInteractions.push_back(seed);
		}
	}
	return summary;
}

Interaction::PairingVec seedOracleBasePairs(const InteractionEnergy & energy,
		const SeedOracleChain & chain)
{
	Interaction::PairingVec basePairs;
	basePairs.reserve(chain.size());
	for (const SeedOraclePair & pair : chain) {
		basePairs.push_back(energy.getBasePair(pair.i1, pair.i2));
	}
	return basePairs;
}

const SeededOracleInteraction * findOptimalTrace(
		const InteractionEnergy & energy, const SeedOracleResult & oracle,
		const Interaction & interaction)
{
	for (const SeededOracleInteraction & candidate : oracle.interactions) {
		if (candidate.energy == oracle.mfe
				&& seedOracleBasePairs(energy, candidate.chain)
						== interaction.basePairs)
		{
			return &candidate;
		}
	}
	return nullptr;
}

bool matchesSeed(const Interaction::Seed & reported,
		const InteractionEnergy & energy, const SeedOracleChain & chain,
		const SeedOccurrence & expected)
{
	return reported.bp_i == energy.getBasePair(chain[expected.firstPair].i1,
				chain[expected.firstPair].i2)
			&& reported.bp_j == energy.getBasePair(chain[expected.lastPair].i1,
					chain[expected.lastPair].i2)
			&& E_equal(reported.energy, expected.totalEnergy);
}

void requireValidSeededMfe(const InteractionEnergy & energy,
		const SeedOracleResult & oracle, const Interaction & interaction)
{
	REQUIRE_FALSE(interaction.isEmpty());
	REQUIRE(interaction.energy == oracle.mfe);
	const SeededOracleInteraction * trace = findOptimalTrace(energy, oracle,
			interaction);
	REQUIRE(trace != nullptr);
	REQUIRE(interaction.seed != nullptr);
	REQUIRE(interaction.seed->size() == trace->seeds.size());
	for (const SeedOccurrence & expected : trace->seeds) {
		const auto reported = std::find_if(interaction.seed->begin(),
				interaction.seed->end(), [&](const Interaction::Seed & seed) {
					return matchesSeed(seed, energy, trace->chain, expected);
				});
		REQUIRE(reported != interaction.seed->end());
	}
}

void compareExactSeededMfePredictors(const InteractionEnergy & energy,
		const SeedConstraint & seedConstraint,
		const OutputConstraint & outputConstraint, const IndexRange & range1,
		const IndexRange & range2, const SeedOracleResult & oracle)
{
	REQUIRE_FALSE(oracle.interactions.empty());

	OutputHandlerInteractionList directOutput(outputConstraint, 1);
	PredictorMfe2dSeed direct(energy, directOutput, nullptr,
			new SeedHandlerMfe(energy, seedConstraint));
	direct.predict(range1, range2);
	REQUIRE_FALSE(directOutput.empty());
	requireValidSeededMfe(energy, oracle, **directOutput.begin());

	OutputHandlerInteractionList extensionOutput(outputConstraint, 1);
	PredictorMfe2dSeedExtension extension(energy, extensionOutput, nullptr,
			new SeedHandlerMfe(energy, seedConstraint));
	extension.predict(range1, range2);
	REQUIRE_FALSE(extensionOutput.empty());
	requireValidSeededMfe(energy, oracle, **extensionOutput.begin());
}

void compareSeedOnlyPredictors(const InteractionEnergy & energy,
		const SeedConstraint & seedConstraint,
		const OutputConstraint & outputConstraint, const IndexRange & range1,
		const IndexRange & range2, const SeedOracleResult & oracle)
{
	const SeedOnlySummary expected = summarizeMfeSeedPerStart(energy, oracle);
	REQUIRE_FALSE(expected.optimalInteractions.empty());

	OutputHandlerInteractionList mfeOutput(outputConstraint, 1);
	PredictorMfeSeedOnly mfe(energy, mfeOutput, nullptr,
			new SeedHandlerMfe(energy, seedConstraint));
	mfe.predict(range1, range2);
	REQUIRE_FALSE(mfeOutput.empty());
	const Interaction & mfeInteraction = **mfeOutput.begin();
	REQUIRE(mfeInteraction.energy == expected.mfe);
	REQUIRE(std::any_of(expected.optimalInteractions.begin(),
			expected.optimalInteractions.end(),
			[&](const SeedOnlyOracleInteraction & candidate) {
				return seedOracleBasePairs(energy, candidate.chain)
						== mfeInteraction.basePairs;
			}));
	REQUIRE(mfeInteraction.seed != nullptr);
	REQUIRE_FALSE(mfeInteraction.seed->empty());

	OutputHandlerInteractionList ensembleOutput(outputConstraint, 1);
	PredictorMfeEnsSeedOnly ensemble(energy, ensembleOutput, nullptr,
			new SeedHandlerMfe(energy, seedConstraint));
	ensemble.predict(range1, range2);
	REQUIRE(ensemble.getZall()
			== Approx(expected.partition).epsilon(1e-12));
}

SeedConstraint makeSeedConstraint(const size_t basePairs,
		const size_t maxUnpairedOverall, const size_t maxUnpaired1,
		const size_t maxUnpaired2, const bool noLP)
{
	return SeedConstraint(basePairs, maxUnpairedOverall, maxUnpaired1,
			maxUnpaired2, E_INF, Accessibility::ED_UPPER_BOUND, E_INF,
			IndexRangeList(), IndexRangeList(), "", false, false, noLP);
}

OutputConstraint makeOutputConstraint(const bool noLP) {
	return OutputConstraint(1, OutputConstraint::OVERLAP_BOTH,
			E_INF, E_INF, false, noLP, false, true, true);
}

void requireUniqueAnchorEnsemblePartition(const InteractionEnergy & energy,
		const SeedConstraint & seedConstraint,
		const OutputConstraint & outputConstraint, const IndexRange & range1,
		const IndexRange & range2, const SeedOracleResult & oracle)
{
	// Seed-extension Z accounting is only asserted in this proved-trivial
	// domain.  General multi-anchor X-mode Z is not an exact-oracle contract.
	REQUIRE(oracle.interactions.size() == 1);
	REQUIRE(oracle.interactions.front().seeds.size() == 1);
	OutputHandlerInteractionList output(outputConstraint, 1);
	PredictorMfeEns2dSeedExtension predictor(energy, output, nullptr,
			new SeedHandlerMfe(energy, seedConstraint));
	predictor.predict(range1, range2);
	REQUIRE(predictor.getZall()
			== Approx(oracle.unionPartition).epsilon(1e-12));
}

} // namespace

TEST_CASE("tiny exhaustive seed oracle for exact predictor families",
		"[PredictorSeedOracle]")
{
	#include "testEasyLoggingSetup.icc"

	SECTION("a unique bulged LP seed is required and traced") {
		RnaSequence target("target", "GAGG");
		RnaSequence query("query", "CCC");
		AccessibilityDisabled targetAcc(target, 0, nullptr);
		AccessibilityDisabled queryAcc(query, 0, nullptr);
		ReverseAccessibility reverseQueryAcc(queryAcc);
		InteractionEnergyBasePair energy(targetAcc, reverseQueryAcc, 2, 2);
		const SeedConstraint seedConstraint = makeSeedConstraint(3, 1, 1, 0,
				false);
		const OutputConstraint outputConstraint = makeOutputConstraint(false);
		const IndexRange targetRange(0, 3);
		const IndexRange queryRange(0, 2);
		const SeedOracleResult oracle = enumerateSeededInteractions(energy,
				seedConstraint, outputConstraint, targetRange, queryRange);

		REQUIRE(oracle.interactions.size() == 1);
		REQUIRE(oracle.interactions.front().seeds.size() == 1);
		const SeedOccurrence & seed = oracle.interactions.front().seeds.front();
		REQUIRE(oracle.interactions.front().chain[seed.lastPair].i1
				-oracle.interactions.front().chain[seed.firstPair].i1+1
				-seedConstraint.getBasePairs() == 1);
		compareExactSeededMfePredictors(energy, seedConstraint,
				outputConstraint, targetRange, queryRange, oracle);
		compareSeedOnlyPredictors(energy, seedConstraint, outputConstraint,
				targetRange, queryRange, oracle);
		requireUniqueAnchorEnsemblePartition(energy, seedConstraint,
				outputConstraint, targetRange, queryRange, oracle);
	}

	SECTION("a noLP bulged seed is two stacked blocks") {
		RnaSequence target("target", "GGAGG");
		RnaSequence query("query", "CCCC");
		AccessibilityDisabled targetAcc(target, 0, nullptr);
		AccessibilityDisabled queryAcc(query, 0, nullptr);
		ReverseAccessibility reverseQueryAcc(queryAcc);
		InteractionEnergyBasePair energy(targetAcc, reverseQueryAcc, 2, 2);
		const SeedConstraint seedConstraint = makeSeedConstraint(4, 1, 1, 0,
				true);
		const OutputConstraint outputConstraint = makeOutputConstraint(true);
		const IndexRange targetRange(0, 4);
		const IndexRange queryRange(0, 3);
		const SeedOracleResult oracle = enumerateSeededInteractions(energy,
				seedConstraint, outputConstraint, targetRange, queryRange);

		REQUIRE(oracle.interactions.size() == 1);
		REQUIRE(hasNoLonelyPairs(oracle.interactions.front().chain, 0, 3));
		REQUIRE_FALSE(isStacked(oracle.interactions.front().chain[1],
				oracle.interactions.front().chain[2]));
		compareExactSeededMfePredictors(energy, seedConstraint,
				outputConstraint, targetRange, queryRange, oracle);
		compareSeedOnlyPredictors(energy, seedConstraint, outputConstraint,
				targetRange, queryRange, oracle);
		requireUniqueAnchorEnsemblePartition(energy, seedConstraint,
				outputConstraint, targetRange, queryRange, oracle);
	}

	SECTION("overlapping equal-energy seeds preserve every witness") {
		RnaSequence target("target", "GGGGG");
		RnaSequence query("query", "CCCCC");
		AccessibilityDisabled targetAcc(target, 0, nullptr);
		AccessibilityDisabled queryAcc(query, 0, nullptr);
		ReverseAccessibility reverseQueryAcc(queryAcc);
		InteractionEnergyBasePair energy(targetAcc, reverseQueryAcc, 2, 2);
		const SeedConstraint seedConstraint = makeSeedConstraint(3, 0, 0, 0,
				false);
		const OutputConstraint outputConstraint = makeOutputConstraint(false);
		const IndexRange full(0, 4);
		const SeedOracleResult oracle = enumerateSeededInteractions(energy,
				seedConstraint, outputConstraint, full, full);

		const auto optimum = std::find_if(oracle.interactions.begin(),
				oracle.interactions.end(), [&](const SeededOracleInteraction & i) {
					return i.energy == oracle.mfe && i.chain.size() == 5;
				});
		REQUIRE(optimum != oracle.interactions.end());
		REQUIRE(optimum->seeds.size() == 3);
		REQUIRE(std::all_of(optimum->seeds.begin(), optimum->seeds.end(),
				[&](const SeedOccurrence & seed) {
					return seed.totalEnergy == optimum->seeds.front().totalEnergy;
				}));

		Z_type naivePerSeedPartition = 0.0;
		for (const SeededOracleInteraction & interaction : oracle.interactions) {
			naivePerSeedPartition += interaction.seeds.size()
					* energy.getBoltzmannWeight(interaction.energy);
		}
		REQUIRE(naivePerSeedPartition > oracle.unionPartition);

		compareExactSeededMfePredictors(energy, seedConstraint,
				outputConstraint, full, full, oracle);
		compareSeedOnlyPredictors(energy, seedConstraint, outputConstraint,
				full, full, oracle);
	}

	SECTION("interactions without a valid seed are not reported") {
		RnaSequence target("target", "GGGG");
		RnaSequence query("query", "CCCC");
		AccessibilityDisabled targetAcc(target, 0, nullptr);
		AccessibilityDisabled queryAcc(query, 0, nullptr);
		ReverseAccessibility reverseQueryAcc(queryAcc);
		InteractionEnergyBasePair energy(targetAcc, reverseQueryAcc, 2, 2);
		const SeedConstraint seedConstraint = makeSeedConstraint(5, 0, 0, 0,
				false);
		const OutputConstraint outputConstraint = makeOutputConstraint(false);
		const IndexRange full(0, 3);
		const SeedOracleResult oracle = enumerateSeededInteractions(energy,
				seedConstraint, outputConstraint, full, full);

		REQUIRE(oracle.interactions.empty());
		REQUIRE(oracle.seedOnlyInteractions.empty());
		REQUIRE(oracle.unionPartition == 0.0);

		OutputHandlerInteractionList directOutput(outputConstraint, 1);
		PredictorMfe2dSeed direct(energy, directOutput, nullptr,
				new SeedHandlerMfe(energy, seedConstraint));
		direct.predict(full, full);
		REQUIRE(directOutput.empty());

		OutputHandlerInteractionList extensionOutput(outputConstraint, 1);
		PredictorMfe2dSeedExtension extension(energy, extensionOutput, nullptr,
				new SeedHandlerMfe(energy, seedConstraint));
		extension.predict(full, full);
		REQUIRE(extensionOutput.empty());

		OutputHandlerInteractionList seedOnlyOutput(outputConstraint, 1);
		PredictorMfeEnsSeedOnly seedOnly(energy, seedOnlyOutput, nullptr,
				new SeedHandlerMfe(energy, seedConstraint));
		seedOnly.predict(full, full);
		REQUIRE(seedOnlyOutput.empty());
		REQUIRE(seedOnly.getZall() == 0.0);
	}

	SECTION("exact seeded predictors can be reused on offset subranges") {
		RnaSequence target("target", "GGGGGG");
		RnaSequence query("query", "CCCCCC");
		AccessibilityDisabled targetAcc(target, 0, nullptr);
		AccessibilityDisabled queryAcc(query, 0, nullptr);
		ReverseAccessibility reverseQueryAcc(queryAcc);
		InteractionEnergyBasePair energy(targetAcc, reverseQueryAcc, 2, 2);
		const SeedConstraint seedConstraint = makeSeedConstraint(3, 0, 0, 0,
				false);
		const OutputConstraint outputConstraint = makeOutputConstraint(false);
		const IndexRange full(0, 5);
		const IndexRange sub1(1, 4);
		const IndexRange sub2(2, 5);
		const SeedOracleResult oracle = enumerateSeededInteractions(energy,
				seedConstraint, outputConstraint, sub1, sub2);
		REQUIRE_FALSE(oracle.interactions.empty());

		SeedRecordingOutput directOutput(outputConstraint);
		PredictorMfe2dSeed direct(energy, directOutput, nullptr,
				new SeedHandlerMfe(energy, seedConstraint));
		direct.predict(full, full);
		direct.predict(sub1, sub2);
		REQUIRE(directOutput.calls.size() == 2);
		requireValidSeededMfe(energy, oracle, directOutput.last());

		SeedRecordingOutput extensionOutput(outputConstraint);
		PredictorMfe2dSeedExtension extension(energy, extensionOutput, nullptr,
				new SeedHandlerMfe(energy, seedConstraint));
		extension.predict(full, full);
		extension.predict(sub1, sub2);
		REQUIRE(extensionOutput.calls.size() == 2);
		requireValidSeededMfe(energy, oracle, extensionOutput.last());

		const SeedOnlySummary seedOnlyExpected = summarizeMfeSeedPerStart(energy,
				oracle);
		SeedRecordingOutput seedOnlyOutput(outputConstraint);
		PredictorMfeEnsSeedOnly seedOnly(energy, seedOnlyOutput, nullptr,
				new SeedHandlerMfe(energy, seedConstraint));
		seedOnly.predict(full, full);
		seedOnly.predict(sub1, sub2);
		REQUIRE(seedOnlyOutput.calls.size() == 2);
		REQUIRE(seedOnlyOutput.last().energy == seedOnlyExpected.mfe);
		REQUIRE(std::any_of(seedOnlyExpected.optimalInteractions.begin(),
				seedOnlyExpected.optimalInteractions.end(),
				[&](const SeedOnlyOracleInteraction & candidate) {
					return seedOracleBasePairs(energy, candidate.chain)
							== seedOnlyOutput.last().basePairs;
				}));
		REQUIRE(seedOnly.getZall()
				== Approx(seedOnlyExpected.partition).epsilon(1e-12));
	}
}

#include "IntaRNA/AccessibilityBasePair.h"
#include "IntaRNA/AccessibilityVrna.h"
#include "IntaRNA/InteractionEnergyVrna.h"
#include "IntaRNA/SeedHandlerNoBulge.h"
#include "IntaRNA/OutputHandlerHub.h"

namespace {
class ExactSeedProbe : public PredictorMfeEns2dSeedExtension {
public:
	using PredictorMfeEns2dSeedExtension::PredictorMfeEns2dSeedExtension;
	std::map<std::array<size_t,4>,Z_type> boundaries;
	bool throwUpdate=false;
	size_t calls=0;
	bool restored() const { return exactContribution==nullptr && !updateZisComplete; }
	size_t retained() const { return Z_partition.size(); }
protected:
	void updateZ(size_t i,size_t j,size_t k,size_t l,Z_type z,bool hybrid) override {
		++calls;
		if (throwUpdate) throw std::runtime_error("injected update failure");
		if (exactContribution) boundaries[{i,k,j,l}]=z*exactContribution->second;
		PredictorMfeEns2dSeedExtension::updateZ(i,j,k,l,z,hybrid);
	}
};
void checkNativeStackPartition(const InteractionEnergy & energy,bool noLP,bool noGU) {
	const SeedConstraint sc=makeSeedConstraint(2,0,0,0,false);
	OutputConstraint oc(3,OutputConstraint::OVERLAP_BOTH,E_INF,E_INF,false,noLP,noGU,true,false);
	IndexRange r1(0,energy.size1()-1), r2(0,energy.size2()-1);
	auto oracle=enumerateSeededInteractions(energy,sc,oc,r1,r2);
	OutputHandlerInteractionList output(oc,3);
	ExactSeedProbe p(energy,output,nullptr,new SeedHandlerMfe(energy,sc));
	p.predict();
	REQUIRE(p.getZall()==Approx(oracle.unionPartition).epsilon(2e-12));
	REQUIRE(output.getZ()==p.getZall());
	REQUIRE(output.usesExactPartition());
	REQUIRE(p.retained()==0);
	std::map<std::array<size_t,4>,Z_type> expected;
	for(const auto & i:oracle.interactions) expected[{i.chain.front().i1,i.chain.front().i2,i.chain.back().i1,i.chain.back().i2}]+=energy.getBoltzmannWeight(i.energy);
	for(const auto & [b,z]:p.boundaries) REQUIRE(z==Approx(expected[b]).epsilon(2e-12));
	for(const auto & [b,z]:expected) REQUIRE(p.boundaries[b]==Approx(z).epsilon(2e-12));
}
}

TEST_CASE("native exact stacked partitions agree beyond the unique-anchor domain", "[PredictorSeedOracle][SeededPartitionFunction]") {
	#include "testEasyLoggingSetup.icc"
	RnaSequence t("t","GUGGCGC"),q("q","GCGCCAC");
	AccessibilityBasePair at(t,6,nullptr), aq(q,5,nullptr);
	ReverseAccessibility ar(aq);
	InteractionEnergyBasePair simple(at,ar,2,1);
	VrnaHandler vrna;
	AccessibilityVrna nativeTarget(t,6,nullptr,vrna),nativeQuery(q,5,nullptr,vrna);
	ReverseAccessibility nativeReverse(nativeQuery);
	InteractionEnergyVrna nearest(nativeTarget,nativeReverse,vrna,2,1);
	for(bool noLP:{false,true}) for(bool noGU:{false,true}) {
		checkNativeStackPartition(simple,noLP,noGU);
		checkNativeStackPartition(nearest,noLP,noGU);
	}
}
TEST_CASE("exact seeded updates retain tiny weights and restore scopes on exceptions", "[PredictorSeedOracle][SeededPartitionFunction]") {
	#include "testEasyLoggingSetup.icc"
	RnaSequence t("t","GGGG"),q("q","CCCC");
	AccessibilityDisabled at(t,0,nullptr),aq(q,0,nullptr);ReverseAccessibility ar(aq);
	InteractionEnergyBasePair energy(at,ar,2,2,false,1,-100,3,10000);
	const auto sc=makeSeedConstraint(2,0,0,0,false);
	const auto oc=makeOutputConstraint(false);
	OutputHandlerInteractionList output(oc,1);
	ExactSeedProbe p(energy,output,nullptr,new SeedHandlerMfe(energy,sc));
	p.throwUpdate=true;REQUIRE_THROWS(p.predict());REQUIRE(p.restored());
	p.throwUpdate=false;p.predict();REQUIRE(p.restored());
	REQUIRE(p.calls>0);
	REQUIRE(p.getZall()>0);
	REQUIRE(p.getZall()<1e-19);
	REQUIRE(output.hasNonzeroZ());
	auto oracle=enumerateSeededInteractions(energy,sc,oc,IndexRange(0,3),IndexRange(0,3));
	REQUIRE(double(p.getZall()/oracle.unionPartition)==Approx(1).epsilon(2e-12));
	OutputHandlerHub hub(oc,false);hub.addOutputHandler(&output);hub.setExactPartition(true);
	REQUIRE(output.usesExactPartition());
	p.predict(IndexRange(1,3),IndexRange(1,3));
	auto sub=enumerateSeededInteractions(energy,sc,oc,IndexRange(1,3),IndexRange(1,3));
	REQUIRE(double(p.getZall()/sub.unionPartition)==Approx(1).epsilon(2e-12));
}

namespace {
class PairTestTracker : public PredictionTracker {
public:
	void updateOptimumCalled(size_t,size_t,size_t,size_t,E_type) override {}
};
}
TEST_CASE("raw native pair masses share the exact forward denominator", "[BasePairProbabilities][PredictorSeedOracle]") {
	#include "testEasyLoggingSetup.icc"
	RnaSequence t("target","GGAGGG"),q("query","CCUCCC");
	AccessibilityBasePair at(t,5,nullptr),aq(q,4,nullptr);ReverseAccessibility ar(aq);
	InteractionEnergyBasePair energy(at,ar,2,1);
	const auto sc=makeSeedConstraint(2,0,0,0,false);
	for(bool noLP:{false,true}) for(bool tracker:{false,true}) {
		auto oc=makeOutputConstraint(noLP);
		BasePairProbabilities result(t.size(),q.size());
		OutputHandlerInteractionList output(oc,1);
		ExactSeedProbe p(energy,output,tracker?new PairTestTracker():nullptr,new SeedHandlerMfe(energy,sc),&result);
		const IndexRange r1(1,5),r2(1,4);
		p.predict(r1,r2);result.finalize();
		auto oracle=enumerateSeededInteractions(energy,sc,oc,r1,r2);
		REQUIRE(result.getZ()==p.getZall());
		REQUIRE(result.getZ()==Approx(oracle.unionPartition).epsilon(2e-12));
		REQUIRE((p.retained()>0)==tracker);
		Matrix<Z_type> expected(t.size(),q.size(),0);
		for(const auto & chain:oracle.interactions) for(const auto & pair:chain.chain) {
			const auto bp=energy.getBasePair(pair.i1,pair.i2);
			expected(bp.first,bp.second)+=energy.getBoltzmannWeight(chain.energy);
		}
		for(size_t i=0;i<t.size();++i) for(size_t j=0;j<q.size();++j)
			REQUIRE(result.rawMasses()(i,j)==Approx(expected(i,j)).epsilon(2e-12));
	}
}
TEST_CASE("raw pair owner excludes region crossings and invalidates partial results", "[BasePairProbabilities][PredictorSeedOracle]") {
	#include "testEasyLoggingSetup.icc"
	RnaSequence t("t","GGGGGG"),q("q","CCCCCC");
	AccessibilityDisabled at(t,0,nullptr),aq(q,0,nullptr);ReverseAccessibility ar(aq);
	InteractionEnergyBasePair energy(at,ar,2,2);
	const auto sc=makeSeedConstraint(2,0,0,0,false);const auto oc=makeOutputConstraint(false);
	BasePairProbabilities owner(6,6);OutputHandlerInteractionList output(oc,1);
	ExactSeedProbe p(energy,output,nullptr,new SeedHandlerMfe(energy,sc),&owner);
	p.predict({0,2},{0,2});const Z_type first=p.getZall();
	p.predict({3,5},{3,5});REQUIRE(p.getZall()==first);
	owner.finalize();REQUIRE(owner.getZ()==2*first);
	REQUIRE(owner.probabilities()(0,0)==0);
	BasePairProbabilities failed(6,6);OutputHandlerInteractionList output2(oc,1);
	ExactSeedProbe p2(energy,output2,nullptr,new SeedHandlerMfe(energy,sc),&failed);
	p2.predict({0,2},{0,2});p2.throwUpdate=true;
	REQUIRE_THROWS(p2.predict({3,5},{3,5}));
	REQUIRE(failed.status()==BasePairProbabilities::Status::failed);
	REQUIRE_THROWS(failed.finalize());REQUIRE(p2.restored());
}

#include "SeededChainOracle.h"
#include "IntaRNA/SeedHandlerExplicit.h"
TEST_CASE("native singleton and mixed explicit seeds count overlapping witnesses once", "[BasePairProbabilities][PredictorSeedOracle]") {
	#include "testEasyLoggingSetup.icc"
	using namespace seeded_oracle;
	RnaSequence t("t","GGGG"),q("q","CCCC");
	AccessibilityDisabled at(t,0,nullptr),aq(q,0,nullptr);ReverseAccessibility ar(aq);
	InteractionEnergyBasePair energy(at,ar,2,2);
	SeedConstraint sc(1,0,0,0,-10000,0,-10000,IndexRangeList(),IndexRangeList(),"1|&4|,2||&2||",true,true,true);
	for(bool noLP:{false,true}) {
		BasePairProbabilities result(4,4);auto oc=makeOutputConstraint(noLP);OutputHandlerInteractionList output(oc,1);
		ExactSeedProbe p(energy,output,nullptr,new SeedHandlerExplicit(energy,sc),&result);
		p.predict();result.finalize();
		Model model{4,4,4,4,noLP,{{{0,0}},{{1,1},{2,2}}},[](Pair){return true;},
			[](Pair){return std::exp(1.L);},[](Pair,Pair){return std::exp(1.L);},[](Pair,Pair){return 1.L;}};
		auto ref=enumerate(model);
		REQUIRE(double(result.getZ())==Approx(double(ref.z)).epsilon(2e-12));
		for(size_t i=0;i<4;++i) for(size_t j=0;j<4;++j)
			REQUIRE(double(result.rawMasses()(i,3-j))==Approx(double(ref.mass[{i,j}])).epsilon(2e-12));
	}
}
