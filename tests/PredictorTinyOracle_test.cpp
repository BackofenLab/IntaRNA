#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/AccessibilityConstraint.h"
#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/OutputHandlerInteractionList.h"
#include "IntaRNA/PredictorMfe2d.h"
#include "IntaRNA/PredictorMfeEns2d.h"
#include "IntaRNA/ReverseAccessibility.h"
#include "IntaRNA/RnaSequence.h"

#include <algorithm>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace IntaRNA;

namespace {

/**
 * A deliberately simple accessibility model for testing site-width and ED
 * filters.  It is independent of the interaction predictors: an accessible
 * interval of width w has ED = w * edPerNt.
 */
class WidthAccessibility final : public Accessibility {
public:
	WidthAccessibility(const RnaSequence & sequence, const size_t maxLength,
			const E_type edPerNt)
	 : Accessibility(sequence, maxLength, NULL)
	 , edPerNt(edPerNt)
	{}

	E_type getED(const size_t from, const size_t to) const override {
		checkIndices(from, to);
		if (to-from+1 > getMaxLength()
				|| !getAccConstraint().isAccessible(from)
				|| !getAccConstraint().isAccessible(to))
		{
			return ED_UPPER_BOUND;
		}
		return static_cast<E_type>((to-from+1) * edPerNt);
	}

private:
	const E_type edPerNt;
};

struct OraclePair {
	size_t i1;
	size_t i2;
};

using OracleChain = std::vector<OraclePair>;

struct OracleInteraction {
	OracleChain chain;
	E_type energy;
};

struct OracleResult {
	std::vector<OracleInteraction> interactions;
	E_type mfe = E_INF;
	Z_type partition = 0.0;
};

size_t lastIndex(const IndexRange & range, const size_t sequenceSize) {
	return range.to == RnaSequence::lastPos ? sequenceSize-1 : range.to;
}

bool hasNoLonelyPairs(const OracleChain & chain) {
	for (size_t i = 0; i < chain.size(); ++i) {
		const bool stackedLeft = i > 0
				&& chain[i-1].i1+1 == chain[i].i1
				&& chain[i-1].i2+1 == chain[i].i2;
		const bool stackedRight = i+1 < chain.size()
				&& chain[i].i1+1 == chain[i+1].i1
				&& chain[i].i2+1 == chain[i+1].i2;
		if (!stackedLeft && !stackedRight) {
			return false;
		}
	}
	return !chain.empty();
}

void recordChain(const InteractionEnergy & energy,
		const OutputConstraint & constraint, const OracleChain & chain,
		const E_type hybridEnergy, OracleResult & result)
{
	if (constraint.noLP && !hasNoLonelyPairs(chain)) {
		return;
	}

	const OraclePair & left = chain.front();
	const OraclePair & right = chain.back();
	if (constraint.noGUend
			&& (energy.isGU(left.i1, left.i2)
					|| energy.isGU(right.i1, right.i2)))
	{
		return;
	}

	const E_type ed1 = energy.getED1(left.i1, right.i1);
	const E_type ed2 = energy.getED2(left.i2, right.i2);
	if (ed1 > constraint.maxED || ed2 > constraint.maxED) {
		return;
	}

	const E_type totalEnergy = energy.getE(left.i1, right.i1,
			left.i2, right.i2, hybridEnergy);
	if (E_isINF(totalEnergy)) {
		return;
	}

	result.interactions.push_back(OracleInteraction{chain, totalEnergy});
	result.mfe = std::min(result.mfe, totalEnergy);
	result.partition += energy.getBoltzmannWeight(totalEnergy);
}

void enumerateExtensions(const InteractionEnergy & energy,
		const OutputConstraint & constraint, const size_t to1, const size_t to2,
		OracleChain & chain, const E_type hybridEnergy, OracleResult & result)
{
	recordChain(energy, constraint, chain, hybridEnergy, result);

	const OraclePair previous = chain.back();
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
			chain.push_back(OraclePair{next1, next2});
			enumerateExtensions(energy, constraint, to1, to2, chain,
					hybridEnergy+loopEnergy, result);
			chain.pop_back();
		}
	}
}

/**
 * Exhaustive reference calculation.  Unlike PredictorMfe2d and
 * PredictorMfeEns2d, this has no dynamic-programming state or recurrence: it
 * visits every strictly increasing pair chain and evaluates it once.
 */
OracleResult enumerateInteractions(const InteractionEnergy & energy,
		const OutputConstraint & constraint, const IndexRange & range1,
		const IndexRange & range2)
{
	OracleResult result;
	const size_t to1 = lastIndex(range1, energy.size1());
	const size_t to2 = lastIndex(range2, energy.size2());
	for (size_t i1 = range1.from; i1 <= to1; ++i1) {
		for (size_t i2 = range2.from; i2 <= to2; ++i2) {
			if (!energy.areComplementary(i1, i2)) {
				continue;
			}
			OracleChain chain(1, OraclePair{i1, i2});
			enumerateExtensions(energy, constraint, to1, to2, chain,
					energy.getE_init(), result);
		}
	}
	return result;
}

Interaction::PairingVec toBasePairs(const InteractionEnergy & energy,
		const OracleChain & chain)
{
	Interaction::PairingVec basePairs;
	basePairs.reserve(chain.size());
	for (const OraclePair & pair : chain) {
		basePairs.push_back(energy.getBasePair(pair.i1, pair.i2));
	}
	return basePairs;
}

bool isOptimalOracleTrace(const InteractionEnergy & energy,
		const OracleResult & oracle, const Interaction & interaction)
{
	for (const OracleInteraction & candidate : oracle.interactions) {
		if (candidate.energy == oracle.mfe
				&& toBasePairs(energy, candidate.chain) == interaction.basePairs)
		{
			return true;
		}
	}
	return false;
}

void compareExactPredictors(const InteractionEnergy & energy,
		const OutputConstraint & constraint, const IndexRange & range1,
		const IndexRange & range2)
{
	const OracleResult oracle = enumerateInteractions(energy, constraint,
			range1, range2);
	REQUIRE_FALSE(oracle.interactions.empty());

	OutputHandlerInteractionList mfeOutput(constraint, 1);
	PredictorMfe2d mfePredictor(energy, mfeOutput, NULL);
	mfePredictor.predict(range1, range2);
	REQUIRE_FALSE(mfeOutput.empty());
	const Interaction & mfe = **mfeOutput.begin();
	REQUIRE(mfe.energy == oracle.mfe);
	REQUIRE(isOptimalOracleTrace(energy, oracle, mfe));

	OutputHandlerInteractionList ensembleOutput(constraint, 1);
	PredictorMfeEns2d ensemblePredictor(energy, ensembleOutput, NULL);
	ensemblePredictor.predict(range1, range2);
	REQUIRE(ensemblePredictor.getZall()
			== Approx(oracle.partition).epsilon(1e-12));
}

class RecordingOutput final : public OutputHandler {
public:
	explicit RecordingOutput(const OutputConstraint & constraint)
	 : OutputHandler(constraint)
	 , zUpdates(0)
	{}

	void add(const Interaction & interaction) override {
		calls.push_back(std::make_unique<Interaction>(interaction));
		++reportedInteractions;
	}

	void incrementZ(const Z_type subZ) override {
		++zUpdates;
		OutputHandler::incrementZ(subZ);
	}

	const Interaction & last() const {
		return *calls.back();
	}

	std::vector<std::unique_ptr<Interaction> > calls;
	size_t zUpdates;
};

std::string reversed(const std::string & sequence) {
	return std::string(sequence.rbegin(), sequence.rend());
}

} // namespace

TEST_CASE("tiny exhaustive oracle for exact predictors", "[PredictorTinyOracle]") {

	#include "testEasyLoggingSetup.icc"

	SECTION("LP chains with internal loops") {
		RnaSequence target("target", "GGAGG");
		RnaSequence query("query", reversed("CCUCC"));
		AccessibilityDisabled targetAcc(target, 0, NULL);
		AccessibilityDisabled queryAcc(query, 0, NULL);
		ReverseAccessibility reverseQueryAcc(queryAcc);
		InteractionEnergyBasePair energy(targetAcc, reverseQueryAcc, 2, 2);
		OutputConstraint constraint(1, OutputConstraint::OVERLAP_BOTH,
				E_INF, E_INF, false, false, false, true, true);

		compareExactPredictors(energy, constraint, IndexRange(0, 4),
				IndexRange(0, 4));
	}

	SECTION("noLP chains are unions of stacked blocks") {
		RnaSequence target("target", "GGGGG");
		RnaSequence query("query", "CCCCC");
		AccessibilityDisabled targetAcc(target, 0, NULL);
		AccessibilityDisabled queryAcc(query, 0, NULL);
		ReverseAccessibility reverseQueryAcc(queryAcc);
		InteractionEnergyBasePair energy(targetAcc, reverseQueryAcc, 2, 2);
		OutputConstraint constraint(1, OutputConstraint::OVERLAP_BOTH,
				E_INF, E_INF, false, true, false, true, true);

		compareExactPredictors(energy, constraint, IndexRange(0, 4),
				IndexRange(0, 4));
	}

	SECTION("GU ends are included or excluded at site level") {
		RnaSequence target("target", "GGGG");
		RnaSequence query("query", reversed("UCCU"));
		AccessibilityDisabled targetAcc(target, 0, NULL);
		AccessibilityDisabled queryAcc(query, 0, NULL);
		ReverseAccessibility reverseQueryAcc(queryAcc);
		InteractionEnergyBasePair energy(targetAcc, reverseQueryAcc, 2, 2);

		OutputConstraint allowGU(1, OutputConstraint::OVERLAP_BOTH,
				E_INF, E_INF, false, false, false, true, true);
		OutputConstraint noGUend(1, OutputConstraint::OVERLAP_BOTH,
				E_INF, E_INF, false, false, true, true, true);
		const OracleResult allowed = enumerateInteractions(energy, allowGU,
				IndexRange(0, 3), IndexRange(0, 3));
		const OracleResult filtered = enumerateInteractions(energy, noGUend,
				IndexRange(0, 3), IndexRange(0, 3));
		REQUIRE(filtered.partition < allowed.partition);

		compareExactPredictors(energy, allowGU, IndexRange(0, 3),
				IndexRange(0, 3));
		compareExactPredictors(energy, noGUend, IndexRange(0, 3),
				IndexRange(0, 3));
	}

	SECTION("loop limits, maximum interaction length, and maxED") {
		RnaSequence target("target", "GGGGG");
		RnaSequence query("query", "CCCCC");
		WidthAccessibility targetAcc(target, 3, Ekcal_2_E(0.10));
		WidthAccessibility queryAcc(query, 4, Ekcal_2_E(0.15));
		ReverseAccessibility reverseQueryAcc(queryAcc);
		InteractionEnergyBasePair energy(targetAcc, reverseQueryAcc, 0, 1);
		OutputConstraint constraint(1, OutputConstraint::OVERLAP_BOTH,
				E_INF, E_INF, false, false, false, true, true,
				Ekcal_2_E(0.45));

		compareExactPredictors(energy, constraint, IndexRange(0, 4),
				IndexRange(0, 4));
	}

	SECTION("offset subranges retain reversed query coordinates") {
		RnaSequence target("target", "AAGGGAA");
		RnaSequence query("query", reversed("ACCCUAA"));
		AccessibilityConstraint targetConstraint(target, "...b...", 0,
				"", "", "");
		AccessibilityDisabled targetAcc(target, 0, &targetConstraint);
		AccessibilityDisabled queryAcc(query, 0, NULL);
		ReverseAccessibility reverseQueryAcc(queryAcc);
		InteractionEnergyBasePair energy(targetAcc, reverseQueryAcc, 2, 2);
		OutputConstraint constraint(1, OutputConstraint::OVERLAP_BOTH,
				E_INF, E_INF, false, false, false, true, true);

		const IndexRange range1(2, 5);
		const IndexRange range2(1, 4);
		compareExactPredictors(energy, constraint, range1, range2);

		const OracleResult oracle = enumerateInteractions(energy, constraint,
				range1, range2);
		for (const OracleInteraction & interaction : oracle.interactions) {
			REQUIRE(std::none_of(interaction.chain.begin(), interaction.chain.end(),
					[](const OraclePair & pair) { return pair.i1 == 3; }));
			for (const Interaction::BasePair & pair
					: toBasePairs(energy, interaction.chain))
			{
				REQUIRE(pair.first >= range1.from);
				REQUIRE(pair.first <= range1.to);
				REQUIRE(pair.second >= query.size()-1-range2.to);
				REQUIRE(pair.second <= query.size()-1-range2.from);
			}
		}
	}

	SECTION("non-seed exact predictors can be reused with new offsets") {
		RnaSequence target("target", "GGGGGG");
		RnaSequence query("query", "CCCCCC");
		AccessibilityDisabled targetAcc(target, 0, NULL);
		AccessibilityDisabled queryAcc(query, 0, NULL);
		ReverseAccessibility reverseQueryAcc(queryAcc);
		InteractionEnergyBasePair energy(targetAcc, reverseQueryAcc, 2, 2);
		OutputConstraint constraint(1, OutputConstraint::OVERLAP_BOTH,
				E_INF, E_INF, false, false, false, true, true);
		const IndexRange full(0, 5);
		const IndexRange sub1(1, 4);
		const IndexRange sub2(2, 5);

		RecordingOutput mfeOutput(constraint);
		PredictorMfe2d mfePredictor(energy, mfeOutput, NULL);
		mfePredictor.predict(full, full);
		mfePredictor.predict(sub1, sub2);
		const OracleResult oracle = enumerateInteractions(energy, constraint,
				sub1, sub2);
		REQUIRE(mfeOutput.calls.size() == 2);
		REQUIRE(mfeOutput.last().energy == oracle.mfe);
		REQUIRE(isOptimalOracleTrace(energy, oracle, mfeOutput.last()));

		RecordingOutput ensembleOutput(constraint);
		PredictorMfeEns2d ensemblePredictor(energy, ensembleOutput, NULL);
		ensemblePredictor.predict(full, full);
		ensemblePredictor.predict(sub1, sub2);
		REQUIRE(ensemblePredictor.getZall()
				== Approx(oracle.partition).epsilon(1e-12));
	}

	SECTION("partition output is merged only when requested") {
		RnaSequence target("target", "GG");
		RnaSequence query("query", "CC");
		AccessibilityDisabled targetAcc(target, 0, NULL);
		AccessibilityDisabled queryAcc(query, 0, NULL);
		ReverseAccessibility reverseQueryAcc(queryAcc);
		InteractionEnergyBasePair energy(targetAcc, reverseQueryAcc, 0, 0);

		OutputConstraint noPartition(1, OutputConstraint::OVERLAP_BOTH,
				E_INF, E_INF, false, false, false, false, false);
		RecordingOutput noPartitionOutput(noPartition);
		PredictorMfe2d noPartitionPredictor(energy, noPartitionOutput, NULL);
		noPartitionPredictor.predict();
		REQUIRE(noPartitionOutput.zUpdates == 0);
		REQUIRE(noPartitionOutput.getZ() == 0.0);
		REQUIRE(noPartitionPredictor.getZall() == 0.0);

		RecordingOutput noPartitionEnsembleOutput(noPartition);
		PredictorMfeEns2d noPartitionEnsemblePredictor(energy,
				noPartitionEnsembleOutput, NULL);
		noPartitionEnsemblePredictor.predict();
		REQUIRE(noPartitionEnsemblePredictor.getZall() > 0.0);
		REQUIRE(noPartitionEnsembleOutput.zUpdates == 0);
		REQUIRE(noPartitionEnsembleOutput.getZ() == 0.0);

		OutputConstraint withPartition(1, OutputConstraint::OVERLAP_BOTH,
				E_INF, E_INF, false, false, false, true, false);
		RecordingOutput partitionOutput(withPartition);
		PredictorMfe2d partitionPredictor(energy, partitionOutput, NULL);
		partitionPredictor.predict();
		REQUIRE(partitionOutput.zUpdates == 1);
		REQUIRE(partitionPredictor.getZall() > 0.0);
		REQUIRE(partitionOutput.getZ() == partitionPredictor.getZall());
	}
}
