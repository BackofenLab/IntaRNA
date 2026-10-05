#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/OutputHandlerInteractionList.h"
#include "IntaRNA/PredictorMfe2dHeuristic.h"
#include "IntaRNA/PredictorMfe2dHeuristicSeed.h"
#include "IntaRNA/PredictorMfeEns2dHeuristic.h"
#include "IntaRNA/PredictorMfe2dHelixBlockHeuristic.h"
#include "IntaRNA/PredictorMfe2dHelixBlockHeuristicSeed.h"
#include "IntaRNA/SeedHandlerMfe.h"
#include "IntaRNA/ReverseAccessibility.h"
#include "IntaRNA/RnaSequence.h"
#include "IntaRNA/SeedConstraint.h"
#include "IntaRNA/SeedHandlerNoBulge.h"

#include <cmath>
#include <memory>
#include <iterator>

using namespace IntaRNA;

namespace {

class InspectableMfeEns2dHeuristic : public PredictorMfeEns2dHeuristic {
public:
	InspectableMfeEns2dHeuristic(const InteractionEnergy & energy,
			OutputHandler & output)
	 : PredictorMfeEns2dHeuristic(energy, output, NULL)
	{}

	Z_type getBoundaryZ(const size_t i1, const size_t j1,
			const size_t i2, const size_t j2) const
	{
		const auto entry = Z_partition.find(Interaction::Boundary(i1, j1, i2, j2));
		return entry == Z_partition.end() ? Z_type(0) : entry->second;
	}
};

void requireThreePairStack(const OutputHandlerInteractionList & output) {
	REQUIRE_FALSE(output.empty());
	const Interaction & interaction = **output.begin();
	REQUIRE(interaction.energy == Ekcal_2_E(-3.0));
	REQUIRE(interaction.basePairs.front() == Interaction::BasePair(0, 3));
	REQUIRE(interaction.basePairs.back() == Interaction::BasePair(2, 1));
}

} // namespace

TEST_CASE("heuristic cells reset their incumbent energy", "[PredictorMfeHeuristicCellState]") {

	#include "testEasyLoggingSetup.icc"

	// In reversed query coordinates, the diagonal interaction is GC-GU-GC.
	// The GU pair is internal and therefore valid when terminal GU pairs are
	// filtered. With stacking-only loops, its energy is exactly -3 kcal/mol.
	RnaSequence target("target", "GGG");
	RnaSequence query("query", "CCUC");
	AccessibilityDisabled targetAcc(target, 0, NULL);
	AccessibilityDisabled queryAcc(query, 0, NULL);
	ReverseAccessibility reverseQueryAcc(queryAcc);
	InteractionEnergyBasePair energy(targetAcc, reverseQueryAcc, 0, 0);
	OutputConstraint constraint(1, OutputConstraint::OVERLAP_BOTH,
			E_INF, E_INF, false, false, true, true, true);

	SECTION("unseeded MFE retains an internal GU continuation") {
		OutputHandlerInteractionList output(constraint, 1);
		PredictorMfe2dHeuristic predictor(energy, output, NULL);

		predictor.predict();

		requireThreePairStack(output);
		REQUIRE((**output.begin()).basePairs.size() == 3);
		REQUIRE((**output.begin()).basePairs.at(1)
				== Interaction::BasePair(1, 2));
	}

	SECTION("ensemble heuristic retains the internal GU boundary partition") {
		OutputHandlerInteractionList output(constraint, 1);
		PredictorMfeEns2dHeuristic predictor(energy, output, NULL);

		predictor.predict();

		requireThreePairStack(output);
		// Ensemble traceback intentionally reports boundaries only.
		REQUIRE((**output.begin()).basePairs.size() == 2);
	}

	SECTION("seeded MFE can extend a seed through an internal GU cell") {
		// Only the leading GC-GC seed is admissible.  Its right extension is
		// GC-GU-GC in the unseeded matrix, so the seed cannot bypass the
		// poisoned cell via a later seed start.
		RnaSequence seedTarget("seedTarget", "GGGG");
		RnaSequence seedQuery("seedQuery", "CCUCC");
		AccessibilityDisabled seedTargetAcc(seedTarget, 0, NULL);
		AccessibilityDisabled seedQueryAcc(seedQuery, 0, NULL);
		ReverseAccessibility reverseSeedQueryAcc(seedQueryAcc);
		InteractionEnergyBasePair seedEnergy(
				seedTargetAcc, reverseSeedQueryAcc, 0, 0);
		IndexRangeList targetSeedRange;
		targetSeedRange.push_back(IndexRange(0, 1));
		IndexRangeList querySeedRange;
		querySeedRange.push_back(IndexRange(0, 1));
		SeedConstraint seedConstraint(2, 0, 0, 0,
				E_INF, Accessibility::ED_UPPER_BOUND, E_INF,
				targetSeedRange, querySeedRange, "", false, false, true);
		OutputConstraint seedOutputConstraint(1, OutputConstraint::OVERLAP_BOTH,
				E_INF, E_INF, false, true, true, true, true);
		OutputHandlerInteractionList output(seedOutputConstraint, 1);
		PredictorMfe2dHeuristicSeed predictor(seedEnergy, output, NULL,
				new SeedHandlerNoBulge(seedEnergy, seedConstraint));

		predictor.predict();

		REQUIRE_FALSE(output.empty());
		const Interaction & interaction = **output.begin();
		REQUIRE(interaction.energy == Ekcal_2_E(-4.0));
		REQUIRE(interaction.basePairs.size() == 4);
		REQUIRE(interaction.basePairs.front() == Interaction::BasePair(0, 4));
		REQUIRE((**output.begin()).basePairs.at(1)
				== Interaction::BasePair(1, 3));
		REQUIRE((**output.begin()).basePairs.at(2)
				== Interaction::BasePair(2, 2));
		REQUIRE(interaction.basePairs.back() == Interaction::BasePair(3, 1));
	}
}

TEST_CASE("ensemble noLP heuristic keeps valid non-direct extensions",
		"[PredictorMfeHeuristicCellState][PredictorMfeEns2dHeuristic]") {

	#include "testEasyLoggingSetup.icc"

	SECTION("a GU mandatory stack can be internal to non-GU boundaries") {
		// Reversed query CUCC gives the unique best GC-GU-GC path at
		// internal boundary (0,2,0,2).  The mandatory second pair is GU,
		// but the actual interaction right end is the following GC pair.
		RnaSequence target("target", "GGG");
		RnaSequence query("query", "CCUC");
		AccessibilityDisabled targetAcc(target, 0, NULL);
		AccessibilityDisabled queryAcc(query, 0, NULL);
		ReverseAccessibility reverseQueryAcc(queryAcc);
		InteractionEnergyBasePair energy(targetAcc, reverseQueryAcc, 0, 0);
		OutputConstraint constraint(1, OutputConstraint::OVERLAP_BOTH,
				E_INF, E_INF, false, true, true, true, true);
		OutputHandlerInteractionList output(constraint, 1);
		InspectableMfeEns2dHeuristic predictor(energy, output);

		predictor.predict();

		const Z_type pathZ = std::exp(3.0);
		const Z_type expectedZall = Z_type(2) * std::exp(2.0) + pathZ;
		REQUIRE(predictor.getBoundaryZ(0, 2, 0, 2)
				== Approx(pathZ).epsilon(1e-12));
		REQUIRE(predictor.getZall() == Approx(expectedZall).epsilon(1e-12));
		REQUIRE_FALSE(output.empty());
		const Interaction & interaction = **output.begin();
		REQUIRE(interaction.energy == Ekcal_2_E(-3.0));
		REQUIRE(interaction.basePairs.size() == 2);
		REQUIRE(interaction.basePairs.front() == Interaction::BasePair(0, 3));
		REQUIRE(interaction.basePairs.back() == Interaction::BasePair(2, 1));
	}

	SECTION("an absent direct extension does not suppress a later bulge") {
		// The first GC-GC block cannot continue directly because A-C at
		// internal (2,2) is impossible.  It can still cross the target A
		// bulge to the second GC-GC block via loop offset (2,1).
		RnaSequence target("target", "GGAGG");
		RnaSequence query("query", "CCCC");
		AccessibilityDisabled targetAcc(target, 0, NULL);
		AccessibilityDisabled queryAcc(query, 0, NULL);
		ReverseAccessibility reverseQueryAcc(queryAcc);
		InteractionEnergyBasePair energy(targetAcc, reverseQueryAcc, 1, 1);
		OutputConstraint constraint(1, OutputConstraint::OVERLAP_BOTH,
				E_INF, E_INF, false, true, true, true, true);
		OutputHandlerInteractionList output(constraint, 1);
		InspectableMfeEns2dHeuristic predictor(energy, output);

		predictor.predict();

		const Z_type pathZ = std::exp(4.0);
		const Z_type expectedZall = Z_type(6) * std::exp(2.0) + pathZ;
		REQUIRE(predictor.getBoundaryZ(0, 4, 0, 3)
				== Approx(pathZ).epsilon(1e-12));
		REQUIRE(predictor.getZall() == Approx(expectedZall).epsilon(1e-12));
		REQUIRE_FALSE(output.empty());
		const Interaction & interaction = **output.begin();
		REQUIRE(interaction.energy == Ekcal_2_E(-4.0));
		REQUIRE(interaction.basePairs.size() == 2);
		REQUIRE(interaction.basePairs.front() == Interaction::BasePair(0, 3));
		REQUIRE(interaction.basePairs.back() == Interaction::BasePair(4, 0));
	}
}

namespace {

std::unique_ptr<Predictor> makeOutputFilterPredictor(const InteractionEnergy & energy,
		OutputHandler & output, const bool seeded, const bool helix) {
	static const SeedConstraint seed(2, 2, 2, 2, E_INF, Accessibility::ED_UPPER_BOUND, E_INF,
			IndexRangeList(), IndexRangeList(), "", false, false, true);
	static const HelixConstraint helixConstraint(2, 4, 2, Accessibility::ED_UPPER_BOUND, E_INF, false);
	if (helix) {
		if (seeded) {
			return std::make_unique<PredictorMfe2dHelixBlockHeuristicSeed>(
					energy, output, nullptr, helixConstraint, new SeedHandlerMfe(energy, seed));
		}
		return std::make_unique<PredictorMfe2dHelixBlockHeuristic>(energy, output, nullptr, helixConstraint);
	}
	if (seeded) {
		return std::make_unique<PredictorMfe2dHeuristicSeed>(energy, output, nullptr,
				new SeedHandlerMfe(energy, seed));
	}
	return std::make_unique<PredictorMfe2dHeuristic>(energy, output, nullptr);
}

// A small accessibility penalty confined to one end of the sequence. The
// second disjoint two-pair site is favorable but exceeds an ED limit of zero.
class EndPenaltyAccessibility : public AccessibilityDisabled {
public:
	EndPenaltyAccessibility(const RnaSequence & sequence, const bool atStart)
	 : AccessibilityDisabled(sequence, 0, nullptr), atStart(atStart) {}

	E_type getED(const size_t from, const size_t to) const override {
		const E_type base = AccessibilityDisabled::getED(from, to);
		return base + ((atStart ? from < 2 : to >= 5) ? Ekcal_2_E(0.25) : 0);
	}
private:
	const bool atStart;
};

} // namespace

TEST_CASE("heuristic suboptimals respect terminal GU constraints", "[PredictorMfeHeuristicCellState][Overlap]") {
	#include "testEasyLoggingSetup.icc"

	for (bool seeded : {false, true}) {
		for (bool helix : {false, true}) {
			for (bool trace : {false, true}) {
				for (size_t offset : {size_t(0), size_t(1)}) {
					RnaSequence target("target", offset ? "NUUGAN" : "UUGA");
					RnaSequence query("query", offset ? "NCAUUN" : "CAUU");
					AccessibilityDisabled targetAcc(target, 0, nullptr);
					AccessibilityDisabled queryAcc(query, 0, nullptr);
					ReverseAccessibility reverseQueryAcc(queryAcc);
					InteractionEnergyBasePair energy(targetAcc, reverseQueryAcc);
					for (auto overlap : {OutputConstraint::OVERLAP_NONE, OutputConstraint::OVERLAP_SEQ1,
							OutputConstraint::OVERLAP_SEQ2, OutputConstraint::OVERLAP_BOTH}) {
						CAPTURE(seeded, helix, trace, offset, overlap);
						OutputConstraint constraint(10, overlap, 0, Ekcal_2_E(100),
								false, false, true, false, trace);
						OutputHandlerInteractionList output(constraint, 10);
						auto predictor = makeOutputFilterPredictor(energy, output, seeded, helix);
						predictor->predict(IndexRange(offset, offset+3), IndexRange(offset, offset+3));
						REQUIRE_FALSE(output.empty());
						for (const Interaction * interaction : output) {
							for (auto bp : {interaction->basePairs.front(), interaction->basePairs.back()}) {
								const char t = target.asString().at(bp.first), q = query.asString().at(bp.second);
								REQUIRE_FALSE((t == 'G' && q == 'U'));
								REQUIRE_FALSE((t == 'U' && q == 'G'));
							}
						}
					}
				}
			}
		}
	}
}

TEST_CASE("heuristic suboptimals respect complete-site accessibility limits", "[PredictorMfeHeuristicCellState][Overlap]") {
	#include "testEasyLoggingSetup.icc"

	RnaSequence target("target", "CCCAACC");
	RnaSequence query("query", "GGAAGGG");
	EndPenaltyAccessibility targetAcc(target, false), queryAcc(query, true);
	ReverseAccessibility reverseQueryAcc(queryAcc);
	InteractionEnergyBasePair energy(targetAcc, reverseQueryAcc, 0, 0);
	for (bool seeded : {false, true}) {
		for (bool helix : {false, true}) {
			for (bool trace : {false, true}) {
				for (auto overlap : {OutputConstraint::OVERLAP_NONE,
						OutputConstraint::OVERLAP_SEQ1, OutputConstraint::OVERLAP_SEQ2}) {
					CAPTURE(seeded, helix, trace, overlap);
					OutputConstraint constraint(10, overlap, 0, Ekcal_2_E(100),
							false, false, false, false, trace, 0);
					OutputHandlerInteractionList output(constraint, 10);
					auto predictor = makeOutputFilterPredictor(energy, output, seeded, helix);
					predictor->predict();
					REQUIRE(std::distance(output.begin(), output.end()) == 1);
					for (const Interaction * interaction : output) {
						REQUIRE(targetAcc.getED(interaction->basePairs.front().first,
								interaction->basePairs.back().first) == 0);
						REQUIRE(queryAcc.getED(interaction->basePairs.back().second,
								interaction->basePairs.front().second) == 0);
						REQUIRE(interaction->energy == Ekcal_2_E(-3));
					}
				}
			}
		}
	}
}
