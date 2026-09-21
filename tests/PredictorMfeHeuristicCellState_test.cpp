#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/OutputHandlerInteractionList.h"
#include "IntaRNA/PredictorMfe2dHeuristic.h"
#include "IntaRNA/PredictorMfe2dHeuristicSeed.h"
#include "IntaRNA/PredictorMfeEns2dHeuristic.h"
#include "IntaRNA/ReverseAccessibility.h"
#include "IntaRNA/RnaSequence.h"
#include "IntaRNA/SeedConstraint.h"
#include "IntaRNA/SeedHandlerNoBulge.h"

#include <cmath>

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
