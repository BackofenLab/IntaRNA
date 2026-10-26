#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/OutputHandler.h"
#include "IntaRNA/OutputHandlerInteractionList.h"
#include "IntaRNA/PredictorMfe2dSeed.h"
#include "IntaRNA/PredictorMfe2dSeedExtension.h"
#include "IntaRNA/PredictorMfe2dHeuristicSeedExtension.h"
#include "IntaRNA/PredictorMfeEns2dSeedExtension.h"
#include "IntaRNA/PredictorMfeEns2dHeuristicSeedExtension.h"
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

const SeededOracleInteraction * findSeededTrace(
		const InteractionEnergy & energy, const SeedOracleResult & oracle,
		const Interaction & interaction)
{
	for (const SeededOracleInteraction & candidate : oracle.interactions) {
		if (candidate.energy == interaction.energy
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

void requireValidSeededInteraction(const InteractionEnergy & energy,
		const SeedOracleResult & oracle, const Interaction & interaction)
{
	REQUIRE_FALSE(interaction.isEmpty());
	REQUIRE(interaction.isValid());
	const SeededOracleInteraction * trace = findSeededTrace(energy, oracle,
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

void requireValidSeededMfe(const InteractionEnergy & energy,
		const SeedOracleResult & oracle, const Interaction & interaction)
{
	REQUIRE(interaction.energy == oracle.mfe);
	requireValidSeededInteraction(energy, oracle, interaction);
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
	// This helper also checks gapped seeds, whose legacy partition path is
	// asserted only in a unique-anchor domain. Stack-only multi-anchor cases
	// have independent exact partition and marginal checks below.
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

template<class PredictorType>
void checkMfeExtensionWorkspaceReuse(const InteractionEnergy & energy,
		const SeedConstraint & seedConstraint,
		const OutputConstraint & outputConstraint, const bool exact)
{
	// Change both dimensions and offsets, including a seedless prediction,
	// before returning to the original shape. Each call must overwrite every
	// active cell even when its backing allocation came from a larger matrix.
	const std::vector<std::pair<IndexRange, IndexRange>> ranges{
		{IndexRange(0, 6), IndexRange(0, 6)},
		{IndexRange(1, 5), IndexRange(2, 5)},
		{IndexRange(2, 2), IndexRange(1, 4)},
		{IndexRange(0, 6), IndexRange(1, 6)},
		{IndexRange(0, 6), IndexRange(0, 6)}
	};
	SeedRecordingOutput reusedOutput(outputConstraint);
	PredictorType reused(energy, reusedOutput, nullptr,
			new SeedHandlerMfe(energy, seedConstraint));
	for (const auto & [range1, range2] : ranges) {
		CAPTURE(range1.from, range1.to, range2.from, range2.to);
		const auto oracle = enumerateSeededInteractions(energy, seedConstraint,
				outputConstraint, range1, range2);
		const size_t previousCalls = reusedOutput.calls.size();
		reused.predict(range1, range2);
		REQUIRE(reusedOutput.calls.size() == previousCalls+1);

		SeedRecordingOutput freshOutput(outputConstraint);
		PredictorType fresh(energy, freshOutput, nullptr,
				new SeedHandlerMfe(energy, seedConstraint));
		fresh.predict(range1, range2);
		REQUIRE(freshOutput.calls.size() == 1);
		const Interaction & actual = reusedOutput.last();
		const Interaction & expected = freshOutput.last();
		REQUIRE(actual.isEmpty() == expected.isEmpty());
		if (oracle.interactions.empty()) {
			REQUIRE(actual.isEmpty());
			continue;
		}

		// The heuristic's retained right extension is not assumed globally
		// optimal. Its reported energy, full traceback and every seed witness
		// must nevertheless describe one independently enumerated interaction.
		if (exact) requireValidSeededMfe(energy, oracle, actual);
		else requireValidSeededInteraction(energy, oracle, actual);
		REQUIRE(actual.energy == expected.energy);
		REQUIRE(actual.basePairs == expected.basePairs);
		REQUIRE(expected.seed != nullptr);
		REQUIRE(*actual.seed == *expected.seed);
	}
	// Equal-energy choices must also survive the intervening range changes.
	REQUIRE(reusedOutput.calls.front()->basePairs == reusedOutput.last().basePairs);
	REQUIRE(reusedOutput.calls.front()->energy == reusedOutput.last().energy);
}

} // namespace

TEST_CASE("MFE seed extensions preserve bounded recurrences and reused workspaces",
		"[PredictorSeedOracle][MfeSeedExtensionWorkspace]")
{
	#include "testEasyLoggingSetup.icc"

	RnaSequence target("target", "GUGGCGG"), query("query", "CCGCCUC");
	AccessibilityBasePair targetAcc(target, 6, nullptr), queryAcc(query, 5, nullptr);
	ReverseAccessibility reverseQuery(queryAcc);
	VrnaHandler vrna;
	AccessibilityVrna nativeTarget(target, 6, nullptr, vrna);
	AccessibilityVrna nativeQuery(query, 5, nullptr, vrna);
	ReverseAccessibility nativeReverse(nativeQuery);
	const SeedConstraint seedConstraint = makeSeedConstraint(2, 0, 0, 0, false);
	for (const auto & [loop1, loop2] : std::vector<std::pair<size_t, size_t>>{
			{0, 0}, {0, 2}, {2, 0}, {2, 1}})
	{
		InteractionEnergyBasePair simple(targetAcc, reverseQuery, loop1, loop2);
		InteractionEnergyVrna nearest(nativeTarget, nativeReverse, vrna, loop1, loop2);
		for (bool native : {false, true}) for (bool noLP : {false, true})
		for (bool noGUend : {false, true}) {
			CAPTURE(loop1, loop2, native, noLP, noGUend);
			const InteractionEnergy & energy = native
					? static_cast<const InteractionEnergy &>(nearest) : simple;
			const OutputConstraint outputConstraint(1,
					OutputConstraint::OVERLAP_BOTH, E_INF, E_INF, false,
					noLP, noGUend, false, true);
			checkMfeExtensionWorkspaceReuse<PredictorMfe2dSeedExtension>(
					energy, seedConstraint, outputConstraint, true);
			checkMfeExtensionWorkspaceReuse<PredictorMfe2dHeuristicSeedExtension>(
					energy, seedConstraint, outputConstraint, false);
		}
	}
}

TEST_CASE("dense MFE seed extension ties retain optimal chains and all seeds",
		"[PredictorSeedOracle][MfeSeedExtensionWorkspace]")
{
	#include "testEasyLoggingSetup.icc"

	RnaSequence target("target", "GGGGGGG"), query("query", "CCCCCCC");
	AccessibilityDisabled targetAcc(target, 6, nullptr), queryAcc(query, 5, nullptr);
	ReverseAccessibility reverseQuery(queryAcc);
	InteractionEnergyBasePair energy(targetAcc, reverseQuery, 2, 1);
	const SeedConstraint seedConstraint = makeSeedConstraint(2, 0, 0, 0, false);
	for (bool noLP : {false, true}) {
		CAPTURE(noLP);
		const auto outputConstraint = makeOutputConstraint(noLP);
		const auto oracle = enumerateSeededInteractions(energy, seedConstraint,
				outputConstraint, IndexRange(0, 6), IndexRange(0, 6));
		// Five stacked pairs attain the maximum query interaction length.
		// Both predictors must find this optimum despite many equal-energy seeds.
		REQUIRE(oracle.mfe == Ekcal_2_E(-5));
		SeedRecordingOutput output(outputConstraint);
		PredictorMfe2dHeuristicSeedExtension predictor(energy, output, nullptr,
				new SeedHandlerMfe(energy, seedConstraint));
		predictor.predict();
		REQUIRE(output.calls.size() == 1);
		requireValidSeededMfe(energy, oracle, output.last());
		REQUIRE(output.last().seed->size() == 4);
		checkMfeExtensionWorkspaceReuse<PredictorMfe2dSeedExtension>(
				energy, seedConstraint, outputConstraint, true);
		checkMfeExtensionWorkspaceReuse<PredictorMfe2dHeuristicSeedExtension>(
				energy, seedConstraint, outputConstraint, true);
	}
}

TEST_CASE("MFE seed extensions include the maximal loop on both sides of a seed",
		"[PredictorSeedOracle][MfeSeedExtensionWorkspace]")
{
	#include "testEasyLoggingSetup.icc"

	for (bool gapInTarget : {false, true}) for (bool noLP : {false, true})
	for (size_t gapLimit : {0, 1}) {
		CAPTURE(gapInTarget, noLP, gapLimit);
		RnaSequence target("target", gapInTarget ? "GGAGGAGG" : "GGGGGG");
		RnaSequence query("query", gapInTarget ? "CCCCCC" : "CCACCACC");
		AccessibilityDisabled targetAcc(target, 0, nullptr), queryAcc(query, 0, nullptr);
		ReverseAccessibility reverseQuery(queryAcc);
		InteractionEnergyBasePair energy(targetAcc, reverseQuery,
				gapInTarget ? gapLimit : 0, gapInTarget ? 0 : gapLimit);
		IndexRangeList seedRange1, seedRange2;
		seedRange1.push_back(gapInTarget ? IndexRange(3, 4) : IndexRange(2, 3));
		seedRange2.push_back(gapInTarget ? IndexRange(2, 3) : IndexRange(3, 4));
		const SeedConstraint seedConstraint(2, 0, 0, 0, E_INF,
				Accessibility::ED_UPPER_BOUND, E_INF, seedRange1, seedRange2,
				"", false, false, noLP);
		const auto outputConstraint = makeOutputConstraint(noLP);
		const IndexRange range1(0, target.size()-1), range2(0, query.size()-1);
		const auto oracle = enumerateSeededInteractions(energy, seedConstraint,
				outputConstraint, range1, range2);
		// Only the central block can be a seed. Extending to the flanking
		// blocks therefore requires a maximal-size loop on each side.
		REQUIRE(oracle.mfe == Ekcal_2_E((gapLimit == 1 ? -6 : -2)));
		SeedRecordingOutput exactOutput(outputConstraint);
		PredictorMfe2dSeedExtension exact(energy, exactOutput, nullptr,
				new SeedHandlerMfe(energy, seedConstraint));
		exact.predict(range1, range2);
		REQUIRE(exactOutput.calls.size() == 1);
		requireValidSeededMfe(energy, oracle, exactOutput.last());
		REQUIRE(exactOutput.last().seed->size() == 1);
		SeedRecordingOutput output(outputConstraint);
		PredictorMfe2dHeuristicSeedExtension predictor(energy, output, nullptr,
				new SeedHandlerMfe(energy, seedConstraint));
		predictor.predict(range1, range2);
		REQUIRE(output.calls.size() == 1);
		requireValidSeededMfe(energy, oracle, output.last());
		REQUIRE(output.last().seed->size() == 1);
	}
}

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
TEST_CASE("native first-hit marginals respect filtered seeds wobble ends and active ranges", "[BasePairProbabilities][PredictorSeedOracle][SeededPartitionFunction]") {
	#include "testEasyLoggingSetup.icc"
	RnaSequence t("target","GUGGGUG"),q("query","CGCCGUC");
	AccessibilityBasePair at(t,6,nullptr),aq(q,5,nullptr);ReverseAccessibility ar(aq);
	InteractionEnergyBasePair energy(at,ar,2,1);
	IndexRangeList seedRange1,seedRange2;
	seedRange1.push_back({2,5});seedRange2.push_back({2,5});
	for(bool noLP:{false,true}) for(bool noGUend:{false,true})
	for(bool noSeedGU:{false,true}) for(size_t offset:{0,1}) {
		CAPTURE(noLP,noGUend,noSeedGU,offset);
		const SeedConstraint sc(2,0,0,0,E_INF,Accessibility::ED_UPPER_BOUND,E_INF,
				seedRange1,seedRange2,"",noSeedGU,false,false);
		const OutputConstraint oc(1,OutputConstraint::OVERLAP_BOTH,E_INF,E_INF,
				false,noLP,noGUend,true,false);
		const IndexRange r1(offset,6),r2(offset,6);
		const auto oracle=enumerateSeededInteractions(energy,sc,oc,r1,r2);
		REQUIRE(oracle.unionPartition>0);
		BasePairProbabilities result(t.size(),q.size());
		OutputHandlerInteractionList output(oc,1);
		ExactSeedProbe p(energy,output,nullptr,new SeedHandlerMfe(energy,sc),&result);
		p.predict(r1,r2);result.finalize();
		REQUIRE(result.getZ()==Approx(oracle.unionPartition).epsilon(2e-12));
		Matrix<Z_type> masses(t.size(),q.size(),0);
		std::map<std::array<size_t,4>,Z_type> boundaries;
		for(const auto & interaction:oracle.interactions) {
			const Z_type weight=energy.getBoltzmannWeight(interaction.energy);
			const auto & first=interaction.chain.front();const auto & last=interaction.chain.back();
			boundaries[{first.i1-offset,first.i2-offset,last.i1-offset,last.i2-offset}]+=weight;
			for(const auto & pair:interaction.chain) {
				const auto bp=energy.getBasePair(pair.i1,pair.i2);
				masses(bp.first,bp.second)+=weight;
			}
		}
		for(size_t i=0;i<t.size();++i) for(size_t j=0;j<q.size();++j)
			REQUIRE(result.rawMasses()(i,j)==Approx(masses(i,j)).epsilon(2e-12));
		for(const auto & [boundary,weight]:boundaries)
			REQUIRE(p.boundaries[boundary]==Approx(weight).epsilon(2e-12));
		for(const auto & [boundary,weight]:p.boundaries)
			REQUIRE(weight==Approx(boundaries[boundary]).epsilon(2e-12));
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
TEST_CASE("exact explicit seed regions do not skip seeds when the iterator cursor collides", "[BasePairProbabilities][PredictorSeedOracle]") {
	#include "testEasyLoggingSetup.icc"
	using namespace seeded_oracle;
	RnaSequence t("target","GGGGGG"),q("query","CCCCCC");
	AccessibilityDisabled at(t,0,nullptr),aq(q,0,nullptr);ReverseAccessibility ar(aq);
	InteractionEnergyBasePair energy(at,ar,2,2);
	const SeedConstraint sc(5,0,0,0,E_INF,Accessibility::ED_UPPER_BOUND,E_INF,
			IndexRangeList(),IndexRangeList(),"1||&5||,4||&2||",false,false,false);
	const Model model{3,3,3,3,false,{{{0,0},{1,1}}},[](Pair){return true;},
		[](Pair){return std::exp(1.L);},[](Pair,Pair){return std::exp(1.L);},
		[](Pair,Pair){return 1.L;}};
	const auto expected=enumerate(model);
	const auto oc=makeOutputConstraint(false);
	BasePairProbabilities result(6,6);
	OutputHandlerInteractionList output(oc,1);
	ExactSeedProbe predictor(energy,output,nullptr,new SeedHandlerExplicit(energy,sc),&result);
	for(size_t offset:{0,3}) {
		CAPTURE(offset);
		predictor.predict({offset,offset+2},{offset,offset+2});
		REQUIRE(double(predictor.getZall())==Approx(double(expected.z)).epsilon(2e-12));
	}
	result.finalize();
	REQUIRE(double(result.getZ())==Approx(double(2*expected.z)).epsilon(2e-12));
	Matrix<Z_type> masses(6,6,0);
	for(size_t offset:{0,3}) for(const auto & [pair,weight]:expected.mass)
		masses(offset+pair.first,5-offset-pair.second)=Z_type(weight);
	for(size_t i=0;i<6;++i) for(size_t j=0;j<6;++j)
		REQUIRE(result.rawMasses()(i,j)==Approx(masses(i,j)).epsilon(2e-12));
}

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

TEST_CASE("native exact nested explicit seeds preserve cropped boundary ensembles", "[BasePairProbabilities][PredictorSeedOracle][SeededPartitionFunction]") {
	#include "testEasyLoggingSetup.icc"
	using namespace seeded_oracle;
	RnaSequence target("target","GGCGGCGG"),query("query","CCGCCGCCC");
	AccessibilityBasePair targetAcc(target,6,nullptr),queryAcc(query,5,nullptr);
	ReverseAccessibility reverseQuery(queryAcc);
	InteractionEnergyBasePair energy(targetAcc,reverseQuery,1,2,false,1,-100,3,37);
	const SeedConstraint sc(1,0,0,0,E_INF,Accessibility::ED_UPPER_BOUND,E_INF,
			IndexRangeList(),IndexRangeList(),
			"2|||||&3|||||,3|&6|,4||&4||,6||&2||,1||&8||",false,false,false);
	// These are independent declared occurrences in increasing energy coordinates.
	// The five-pair seed starts before, but completes after, the inner singleton.
	const std::vector<Chain> declared{
		{{1,2},{2,3},{3,4},{4,5},{5,6}},{{2,3}},
		{{3,4},{4,5}},{{5,6},{6,7}},{{0,0},{1,1}}};
	for(bool cropped:{false,true}) for(bool noLP:{false,true}) {
		CAPTURE(cropped,noLP);
		const IndexRange r1=cropped?IndexRange(1,6):IndexRange(0,7);
		const IndexRange r2=cropped?IndexRange(2,7):IndexRange(0,8);
		const size_t n=r1.to-r1.from+1,m=r2.to-r2.from+1;
		auto absolute=[&](Pair p){return Pair{p.first+r1.from,p.second+r2.from};};
		auto weight=[&](E_type e){return E_isINF(e)?0.L:static_cast<long double>(energy.getBoltzmannWeight(e));};
		Model model{n,m,6,5,noLP,{},
			[&](Pair p){p=absolute(p);return energy.areComplementary(p.first,p.second);},
			[&](Pair){return weight(energy.getE_init());},
			[&](Pair p,Pair q){
				if(q.first-p.first-1>1 || q.second-p.second-1>2) return 0.L;
				p=absolute(p);q=absolute(q);
				return weight(energy.getE_interLeft(p.first,q.first,p.second,q.second));
			},
			[&](Pair p,Pair q){
				p=absolute(p);q=absolute(q);
				return weight(energy.getE(p.first,q.first,p.second,q.second,E_type(0)));
			}};
		for(const auto & seed:declared) {
			if(seed.front().first<r1.from || seed.front().second<r2.from
					|| seed.back().first>r1.to || seed.back().second>r2.to) continue;
			Chain local;
			for(const auto p:seed) local.push_back({p.first-r1.from,p.second-r2.from});
			model.seeds.push_back(local);
		}
		const auto expected=enumerate(model);
		REQUIRE(expected.z>0);
		const auto oc=makeOutputConstraint(noLP);
		BasePairProbabilities result(target.size(),query.size());
		OutputHandlerInteractionList output(oc,1);
		ExactSeedProbe predictor(energy,output,nullptr,new SeedHandlerExplicit(energy,sc),&result);
		predictor.predict(r1,r2);result.finalize();
		REQUIRE(double(result.getZ())==Approx(double(expected.z)).epsilon(2e-12));
		REQUIRE(result.getZ()==predictor.getZall());
		REQUIRE(output.getZ()==predictor.getZall());
		for(const auto & [boundary,z]:expected.boundaries)
			REQUIRE(double(predictor.boundaries[boundary])==Approx(double(z)).epsilon(2e-12));
		for(const auto & [boundary,z]:predictor.boundaries) {
			const auto found=expected.boundaries.find(boundary);
			REQUIRE(found!=expected.boundaries.end());
			REQUIRE(double(z)==Approx(double(found->second)).epsilon(2e-12));
		}
		Matrix<Z_type> masses(target.size(),query.size(),0);
		for(const auto & [pair,z]:expected.mass) {
			const Pair p=absolute(pair);
			masses(p.first,query.size()-1-p.second)=Z_type(z);
		}
		for(size_t i=0;i<target.size();++i) for(size_t j=0;j<query.size();++j)
			REQUIRE(result.rawMasses()(i,j)==Approx(masses(i,j)).epsilon(2e-12));

		// Ensemble output encodes the optimal boundary, whose energy is the
		// sum of the rounded hybrid free energy and its native boundary energy.
		E_type best=E_INF;
		std::map<Boundary,E_type> energies;
		for(const auto & [boundary,z]:expected.boundaries) {
			if(z==0) continue;
			const Pair p=absolute({boundary[0],boundary[1]});
			const Pair q=absolute({boundary[2],boundary[3]});
			const E_type boundaryE=energy.getE(p.first,q.first,p.second,q.second,E_type(0));
			const E_type ensemble=energy.getE(Z_type(z/weight(boundaryE)))+boundaryE;
			energies[boundary]=ensemble;best=std::min(best,ensemble);
		}
		REQUIRE_FALSE(output.empty());
		const Interaction & reported=**output.begin();
		REQUIRE(reported.isValid());
		REQUIRE(reported.energy==best);
		const auto & first=reported.basePairs.front();
		const auto & last=reported.basePairs.back();
		const Boundary reportedBoundary{first.first-r1.from,query.size()-1-first.second-r2.from,
				last.first-r1.from,query.size()-1-last.second-r2.from};
		REQUIRE(energies.at(reportedBoundary)==best);
	}
}

namespace {
class AnchoredHeuristicProbe : public PredictorMfeEns2dHeuristicSeedExtension {
public:
	using PredictorMfeEns2dHeuristicSeedExtension::PredictorMfeEns2dHeuristicSeedExtension;
	using Pair=Interaction::BasePair;
	std::map<Pair,Pair> choices;
	std::map<Pair,size_t> fillCalls;
	size_t rankingCalls=0;
	std::map<seeded_oracle::Boundary,std::pair<Z_type,Z_type>> coefficients;
	std::map<seeded_oracle::Boundary,size_t> completeCalls;
protected:
	void fillHybridZ_right(size_t sj1,size_t sj2,size_t si1,size_t si2) override;
	void updateOptRightZ(size_t i1,size_t j1,size_t i2,size_t j2,E_type hybrid) override;
	void updateZ(size_t i1,size_t j1,size_t i2,size_t j2,Z_type z,bool hybrid) override;
};

void AnchoredHeuristicProbe::fillHybridZ_right(size_t sj1,size_t sj2,size_t si1,size_t si2)
{
	++fillCalls[{si1,si2}];
	PredictorMfeEns2dHeuristicSeedExtension::fillHybridZ_right(sj1,sj2,si1,si2);
	// A custom native hook fixes the retained endpoint after ordinary ranking.
	// The shared kernel must use that choice without selecting it again.
	const auto selected=choices.at({si1,si2});
	j1opt=selected.first;j2opt=selected.second;E_right_opt=0;
}

void AnchoredHeuristicProbe::updateOptRightZ(size_t i1,size_t j1,size_t i2,size_t j2,E_type hybrid)
{
	++rankingCalls;
	PredictorMfeEns2dHeuristicSeedExtension::updateOptRightZ(i1,j1,i2,j2,hybrid);
}

void AnchoredHeuristicProbe::updateZ(size_t i1,size_t j1,size_t i2,size_t j2,Z_type z,bool hybrid)
{
	// Right-selection hooks must not publish partial partitions. The common
	// anchor contraction reports each complete boundary with its final factor.
	REQUIRE(exactContribution!=nullptr);
	REQUIRE(updateZisComplete);
	const seeded_oracle::Boundary boundary{i1,i2,j1,j2};
	coefficients[boundary]=*exactContribution;
	++completeCalls[boundary];
	PredictorMfeEns2dHeuristicSeedExtension::updateZ(i1,j1,i2,j2,z,hybrid);
}
}

TEST_CASE("native anchored heuristics preserve right-choice hooks and complete boundary coefficients", "[BasePairProbabilities][PredictorSeedOracle][SeededPartitionFunction]") {
	#include "testEasyLoggingSetup.icc"
	using namespace seeded_oracle;
	RnaSequence target("target","GGGGGGG"),query("query","CCCCCCC");
	AccessibilityBasePair targetAcc(target,6,nullptr),queryAcc(query,5,nullptr);
	ReverseAccessibility reverseQuery(queryAcc);
	InteractionEnergyBasePair energy(targetAcc,reverseQuery,2,1);
	const SeedConstraint sc(2,0,0,0,E_INF,Accessibility::ED_UPPER_BOUND,E_INF,
			IndexRangeList(),IndexRangeList(),"1|||||&3|||||,2||&5||,4||&3||",false,false,false);
	const std::vector<Chain> declared{
		{{0,0},{1,1},{2,2},{3,3},{4,4}},{{1,1},{2,2}},{{3,3},{4,4}}};
	const std::map<Pair,Pair> choices{{{0,0},{6,6}},{{1,1},{5,5}},{{3,3},{6,6}}};
	auto weight=[&](E_type e){return E_isINF(e)?0.L:static_cast<long double>(energy.getBoltzmannWeight(e));};
	for(bool noLP:{false,true}) {
		CAPTURE(noLP);
		Model model{7,7,6,5,noLP,declared,[](Pair){return true;},
			[&](Pair){return weight(energy.getE_init());},
			[&](Pair p,Pair q){
				if(q.first-p.first-1>2 || q.second-p.second-1>1) return 0.L;
				return weight(energy.getE_interLeft(p.first,q.first,p.second,q.second));
			},
			[&](Pair p,Pair q){return weight(energy.getE(p.first,q.first,p.second,q.second,E_type(0)));},
			[&](Pair p,Pair){return choices.at(p);}};
		const auto expected=enumerate(model);
		REQUIRE(expected.z>0);
		auto hybridModel=model;hybridModel.boundary=[](Pair,Pair){return 1.L;};
		const auto hybrids=enumerate(hybridModel);
		const auto oc=makeOutputConstraint(noLP);
		std::map<Boundary,std::pair<Z_type,Z_type>> previousCoefficients;
		for(bool pairs:{false,true}) {
			CAPTURE(pairs);
			BasePairProbabilities result(target.size(),query.size());
			OutputHandlerInteractionList output(oc,1);
			AnchoredHeuristicProbe predictor(energy,output,nullptr,new SeedHandlerExplicit(energy,sc),pairs?&result:nullptr);
			predictor.choices=choices;predictor.predict();
			REQUIRE(predictor.rankingCalls>0);
			REQUIRE(predictor.fillCalls.size()==declared.size());
			for(const auto & seed:declared) REQUIRE(predictor.fillCalls.at(seed.front())==1);
			REQUIRE(double(predictor.getZall())==Approx(double(expected.z)).epsilon(2e-12));
			REQUIRE(output.getZ()==predictor.getZall());
			REQUIRE(predictor.coefficients.size()==hybrids.boundaries.size());
			for(const auto & [boundary,h]:hybrids.boundaries) {
				REQUIRE(predictor.completeCalls.at(boundary)==1);
				const auto [actualH,b]=predictor.coefficients.at(boundary);
				REQUIRE(double(actualH)==Approx(double(h)).epsilon(2e-12));
				REQUIRE(double(b)==Approx(double(model.boundary({boundary[0],boundary[1]},{boundary[2],boundary[3]}))).epsilon(2e-12));
			}
			if(!pairs) { previousCoefficients=predictor.coefficients;continue; }
			REQUIRE(predictor.coefficients==previousCoefficients);
			result.finalize();REQUIRE(result.isApproximate());
			REQUIRE(result.getZ()==predictor.getZall());
			for(size_t i=0;i<target.size();++i) for(size_t j=0;j<query.size();++j) {
				const auto found=expected.mass.find({i,j});
				REQUIRE(double(result.rawMasses()(i,query.size()-1-j))==Approx(double(found==expected.mass.end()?0:found->second)).epsilon(2e-12));
			}
		}
	}
}

namespace {
using PerSeedSite = std::pair<Interaction::BasePair,Interaction::BasePair>;
using PerSeedResults = std::map<PerSeedSite,E_type>;

PerSeedResults perSeedExpected(const InteractionEnergy & energy, const SeedOracleResult & oracle,
		const size_t quota, const bool ensemble, const bool seedOnly = false)
{
	// Enumerate chains independently, then group MFE conformations by every
	// contained seed, or ensemble weights by the first (left-most) seed only.
	using Pair = Interaction::BasePair;
	std::map<Pair,std::map<PerSeedSite,E_type>> minima;
	std::map<Pair,std::map<PerSeedSite,Z_type>> partitions;
	for (const auto & interaction : oracle.interactions) {
		if (seedOnly && interaction.chain.size() != 2) continue;
		const auto & first = interaction.chain.front(), & last = interaction.chain.back();
		const PerSeedSite site(energy.getBasePair(first.i1,first.i2),energy.getBasePair(last.i1,last.i2));
		for (size_t k=0;k<interaction.seeds.size();++k) {
			if (ensemble && k) break;
			const auto & start = interaction.chain[interaction.seeds[k].firstPair];
			const Pair seed = energy.getBasePair(start.i1,start.i2);
			auto [it,inserted] = minima[seed].try_emplace(site,interaction.energy);
			if (!inserted) it->second = std::min(it->second,interaction.energy);
			partitions[seed][site] += energy.getBoltzmannWeight(interaction.energy);
		}
	}
	PerSeedResults result;
	for (const auto & [seed,sites] : minima) {
		std::vector<Interaction> ranked;
		for (const auto & [site,mfe] : sites) {
			Interaction candidate(energy.getAccessibility1().getSequence(),energy.getAccessibility2().getAccessibilityOrigin().getSequence());
			candidate.basePairs = {site.first,site.second};
			candidate.energy = mfe;
			if (ensemble) {
				const size_t i=energy.getIndex1(site.first),j=energy.getIndex1(site.second);
				const size_t k=energy.getIndex2(site.first),l=energy.getIndex2(site.second);
				const E_type boundary = energy.getE(i,j,k,l,E_type(0));
				candidate.energy = boundary + energy.getE(partitions.at(seed).at(site)/energy.getBoltzmannWeight(boundary));
			}
			ranked.push_back(candidate);
		}
		std::sort(ranked.begin(),ranked.end());
		for (size_t k=0;k<std::min(quota,ranked.size());++k) {
			const auto & candidate=ranked[k];
			const PerSeedSite site(candidate.basePairs.front(),candidate.basePairs.back());
			auto [it,inserted]=result.try_emplace(site,candidate.energy);
			if (!inserted) it->second=std::min(it->second,candidate.energy);
		}
	}
	return result;
}

PerSeedResults perSeedActual(const OutputHandlerInteractionList & output)
{
	PerSeedResults result;
	for (const Interaction * candidate : output) {
		REQUIRE(candidate->isValid());
		const auto [it,inserted]=result.emplace(PerSeedSite(candidate->basePairs.front(),candidate->basePairs.back()),candidate->energy);
		REQUIRE(inserted);
	}
	return result;
}
} // namespace

TEST_CASE("per-seed exact MFE and anchor ensembles match independent chain enumeration", "[PredictorSeedOracle][OutputPerSeed]") {
	#include "testEasyLoggingSetup.icc"
	RnaSequence target("target","AGCGGCG"),query("query","CCGCCGC");
	AccessibilityDisabled tAcc(target,6,nullptr),qAcc(query,6,nullptr);
	ReverseAccessibility reverse(qAcc);
	InteractionEnergyBasePair energy(tAcc,reverse,2,2,false,1,-100,3,37);
	const auto sc=makeSeedConstraint(2,0,0,0,false);
	for (bool noLP : {false,true}) for (bool cropped : {false,true}) for (size_t quota : {1,2,1000}) {
		CAPTURE(noLP,cropped,quota);
		const IndexRange r1=cropped?IndexRange(1,5):IndexRange(0,6);
		const IndexRange r2=cropped?IndexRange(1,6):IndexRange(0,6);
		const OutputConstraint c(1000,OutputConstraint::OVERLAP_BOTH,E_INF,E_INF,false,noLP,false,true,true,Accessibility::ED_UPPER_BOUND,quota);
		const auto oracle=enumerateSeededInteractions(energy,sc,c,r1,r2);
		REQUIRE_FALSE(oracle.interactions.empty());
		for (char model : {'S','X','P'}) {
			CAPTURE(model);
			OutputHandlerInteractionList output(c,1000);
			std::unique_ptr<Predictor> predictor;
			if (model=='S') predictor=std::make_unique<PredictorMfe2dSeed>(energy,output,nullptr,new SeedHandlerMfe(energy,sc));
			if (model=='X') predictor=std::make_unique<PredictorMfe2dSeedExtension>(energy,output,nullptr,new SeedHandlerMfe(energy,sc));
			if (model=='P') predictor=std::make_unique<PredictorMfeEns2dSeedExtension>(energy,output,nullptr,new SeedHandlerMfe(energy,sc));
			predictor->predict(r1,r2);
			const auto actual=perSeedActual(output), expected=perSeedExpected(energy,oracle,quota,model=='P');
			std::ostringstream diagnostic;
			for (const auto & [site,e] : actual) diagnostic << " actual " << site.first.first << ',' << site.first.second << '-' << site.second.first << ',' << site.second.second << ':' << e;
			for (const auto & [site,e] : expected) diagnostic << " expected " << site.first.first << ',' << site.first.second << '-' << site.second.first << ',' << site.second.second << ':' << e;
			CAPTURE(diagnostic.str());
			REQUIRE(actual==expected);
			if (model=='P') REQUIRE(double(predictor->getZall())==Approx(double(oracle.unionPartition)).epsilon(1e-12));
			if (model!='P') for (const Interaction * interaction : output) {
				// Each retained structure must be a physical chain of its stated
				// energy, even when another seed has a better trace at that site.
				bool found=false;
				for (const auto & chain : oracle.interactions) {
					Interaction::PairingVec pairs;
					for (const auto & p : chain.chain) pairs.push_back(energy.getBasePair(p.i1,p.i2));
					found=found || (pairs==interaction->basePairs && chain.energy==interaction->energy);
				}
				REQUIRE(found);
			}
		}
	}
}

TEST_CASE("per-seed seed-only ranking retains every seed site", "[PredictorSeedOracle][OutputPerSeed]") {
	#include "testEasyLoggingSetup.icc"
	RnaSequence target("target","GGGGG"),query("query","CCCCC");
	AccessibilityDisabled tAcc(target,5,nullptr),qAcc(query,5,nullptr);ReverseAccessibility reverse(qAcc);
	InteractionEnergyBasePair energy(tAcc,reverse,1,1);
	const auto sc=makeSeedConstraint(2,0,0,0,false);
	for (bool ensemble : {false,true}) {
		const OutputConstraint c(1000,OutputConstraint::OVERLAP_BOTH,E_INF,E_INF,false,false,false,true,true,Accessibility::ED_UPPER_BOUND,1);
		const auto oracle=enumerateSeededInteractions(energy,sc,c,IndexRange(0,4),IndexRange(0,4));
		OutputHandlerInteractionList output(c,1000);
		std::unique_ptr<Predictor> predictor;
		if (ensemble) predictor=std::make_unique<PredictorMfeEnsSeedOnly>(energy,output,nullptr,new SeedHandlerMfe(energy,sc));
		else predictor=std::make_unique<PredictorMfeSeedOnly>(energy,output,nullptr,new SeedHandlerMfe(energy,sc));
		predictor->predict();
		REQUIRE(perSeedActual(output)==perSeedExpected(energy,oracle,1,false,true));
	}
}

namespace {
class PerSeedHeuristicChoices : public PredictorMfeEns2dHeuristicSeedExtension {
public:
	using PredictorMfeEns2dHeuristicSeedExtension::PredictorMfeEns2dHeuristicSeedExtension;
	std::map<Interaction::BasePair,Interaction::BasePair> choices() const;
};
inline std::map<Interaction::BasePair,Interaction::BasePair> PerSeedHeuristicChoices::choices() const
{ return {stackSeedRight.begin(),stackSeedRight.end()}; }

PerSeedResults expectedAnchorSites(const InteractionEnergy & energy, const seeded_oracle::Result & oracle, size_t quota)
{
	PerSeedResults expected;
	for (const auto & [seed,sites] : oracle.anchors) {
		std::vector<Interaction> ranked;
		for (const auto & [site,z] : sites) {
			Interaction candidate(energy.getAccessibility1().getSequence(),energy.getAccessibility2().getAccessibilityOrigin().getSequence());
			candidate.basePairs={energy.getBasePair(site[0],site[1]),energy.getBasePair(site[2],site[3])};
			candidate.energy=energy.getE(Z_type(z)); // disabled accessibility and zero boundary terms
			ranked.push_back(candidate);
		}
		std::sort(ranked.begin(),ranked.end());
		for (size_t i=0;i<std::min(quota,ranked.size());++i) {
			const auto & c=ranked[i];
			auto [it,inserted]=expected.try_emplace(PerSeedSite(c.basePairs.front(),c.basePairs.back()),c.energy);
			if (!inserted) it->second=std::min(it->second,c.energy);
		}
	}
	return expected;
}
}

TEST_CASE("per-seed ensemble selection uses anchor partitions in heuristic and bulged modes", "[PredictorSeedOracle][OutputPerSeed]") {
	#include "testEasyLoggingSetup.icc"
	using namespace seeded_oracle;
	RnaSequence target("target","GGGGGG"),query("query","CCCCCC");
	AccessibilityDisabled tAcc(target,6,nullptr),qAcc(query,6,nullptr);ReverseAccessibility reverse(qAcc);
	InteractionEnergyBasePair energy(tAcc,reverse,2,2);
	auto weight=[&](E_type e) { return E_isINF(e)?0.L:static_cast<long double>(energy.getBoltzmannWeight(e)); };
	for (bool gapped : {false,true}) for (bool heuristic : {false,true}) for (size_t quota : {1,2,1000}) {
		CAPTURE(gapped,heuristic,quota);
		const SeedConstraint sc(gapped?3:2,0,0,0,E_INF,Accessibility::ED_UPPER_BOUND,E_INF,
				IndexRangeList(),IndexRangeList(),gapped?"2|||&2||.|":"",false,false,false);
		Model model{6,6,6,6,false,{},[](Pair){return true;},
				[&](Pair){return weight(energy.getE_init());},
				[&](Pair p,Pair q){return weight(energy.getE_interLeft(p.first,q.first,p.second,q.second));},
				[](Pair,Pair){return 1.L;}};
		if (gapped) model.seeds={{{1,1},{2,3},{3,4}}};
		else for (size_t i=0;i<5;++i) for (size_t j=0;j<5;++j) model.seeds.push_back({{i,j},{i+1,j+1}});
		OutputConstraint c(1000,OutputConstraint::OVERLAP_BOTH,E_INF,E_INF,false,false,false,true,true,Accessibility::ED_UPPER_BOUND,quota);
		OutputHandlerInteractionList output(c,1000);
		std::unique_ptr<PredictorMfeEns2dSeedExtension> predictor;
		SeedHandler * seeds=gapped?static_cast<SeedHandler *>(new SeedHandlerExplicit(energy,sc)):new SeedHandlerMfe(energy,sc);
		if (heuristic) predictor=std::make_unique<PerSeedHeuristicChoices>(energy,output,nullptr,seeds);
		else predictor=std::make_unique<PredictorMfeEns2dSeedExtension>(energy,output,nullptr,seeds);
		predictor->predict();
		if (heuristic) {
			// Freeze the native right-endpoint choice, then independently enumerate
			// physical chains admitted by that heuristic's exact retained domain.
			if (gapped) model.selectRight=[](Pair,Pair){return Pair(4,5);};
			else {
				const auto choices=static_cast<PerSeedHeuristicChoices *>(predictor.get())->choices();
				model.selectRight=[choices](Pair p,Pair){return choices.at(p);};
			}
		}
		const auto expected=enumerate(model);
		REQUIRE(perSeedActual(output)==expectedAnchorSites(energy,expected,quota));
		REQUIRE(double(predictor->getZall())==Approx(double(expected.z)).epsilon(1e-12));
	}
}
