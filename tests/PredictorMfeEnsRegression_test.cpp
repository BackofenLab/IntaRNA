#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/OutputHandlerInteractionList.h"
#include "IntaRNA/PredictorMfeEns2d.h"
#include "IntaRNA/PredictorMfeEns2dHeuristic.h"
#include "IntaRNA/PredictorMfeEns2dHeuristicSeedExtension.h"
#include "IntaRNA/PredictorMfeEns2dSeedExtension.h"
#include "IntaRNA/ReverseAccessibility.h"
#include "IntaRNA/RnaSequence.h"
#include "IntaRNA/SeedHandlerNoBulge.h"

#include <cmath>

using namespace IntaRNA;

template <class PredictorType>
class InspectableEnsemblePredictor : public PredictorType {
public:
	InspectableEnsemblePredictor(
			const InteractionEnergy & energy,
			OutputHandler & output,
			SeedHandler * seedHandler)
	 : PredictorType(energy, output, NULL, seedHandler)
	{}

	size_t getPartitionCount() const {
		return this->Z_partition.size();
	}
};

TEST_CASE("ensemble predictor regressions", "[PredictorMfeEns]") {

	#include "testEasyLoggingSetup.icc"

	SECTION("a pruning heuristic never exceeds the exact noLP partition") {
		RnaSequence target("target", "GGGG");
		RnaSequence query("query", "CCCC");
		AccessibilityDisabled targetAcc(target, 0, NULL);
		AccessibilityDisabled queryAcc(query, 0, NULL);
		ReverseAccessibility reverseQueryAcc(queryAcc);
		InteractionEnergyBasePair energy(targetAcc, reverseQueryAcc);

		OutputConstraint constraint(1, OutputConstraint::OVERLAP_BOTH,
				E_INF, E_INF, false, true, false, true, false);
		OutputHandlerInteractionList exactOut(constraint, 1);
		OutputHandlerInteractionList heuristicOut(constraint, 1);
		PredictorMfeEns2d exact(energy, exactOut, NULL);
		PredictorMfeEns2dHeuristic heuristic(energy, heuristicOut, NULL);

		exact.predict();
		heuristic.predict();

		REQUIRE(exact.getZall() > 0.0);
		REQUIRE(heuristic.getZall() <= exact.getZall() * (1.0 + 1e-12));
	}

	SECTION("terminal GU filtering also applies to the partition") {
		RnaSequence target("target", "G");
		RnaSequence query("query", "U");
		AccessibilityDisabled targetAcc(target, 0, NULL);
		AccessibilityDisabled queryAcc(query, 0, NULL);
		ReverseAccessibility reverseQueryAcc(queryAcc);
		InteractionEnergyBasePair energy(targetAcc, reverseQueryAcc);

		OutputConstraint constraint(1, OutputConstraint::OVERLAP_BOTH,
				E_INF, E_INF, false, false, true, true, false);
		OutputHandlerInteractionList out(constraint, 1);
		PredictorMfeEns2d predictor(energy, out, NULL);

		predictor.predict();

		REQUIRE(predictor.getZall() == 0.0);
	}

	SECTION("noLP seed extension multiplies stacked partition factors") {
		RnaSequence target("target", "GGG");
		RnaSequence query("query", "CCC");
		AccessibilityDisabled targetAcc(target, 0, NULL);
		AccessibilityDisabled queryAcc(query, 0, NULL);
		ReverseAccessibility reverseQueryAcc(queryAcc);
		InteractionEnergyBasePair energy(targetAcc, reverseQueryAcc);
		IndexRangeList targetSeedRange;
		targetSeedRange.push_back(IndexRange(1, 2));
		IndexRangeList querySeedRange;
		querySeedRange.push_back(IndexRange(1, 2));
		SeedConstraint seedConstraint(2, 0, 0, 0,
				E_INF, Accessibility::ED_UPPER_BOUND, E_INF,
				targetSeedRange, querySeedRange, "", false, false, true);
		OutputConstraint constraint(1, OutputConstraint::OVERLAP_BOTH,
				E_INF, E_INF, false, true, false, true, false);
		OutputHandlerInteractionList out(constraint, 1);
		PredictorMfeEns2dSeedExtension predictor(
				energy, out, NULL, new SeedHandlerNoBulge(energy, seedConstraint));

		predictor.predict();

		// The only seed starts at internal (1,1), contributing exp(2).
		// Its sole noLP extension stacks (0,0) to the left, contributing exp(3).
		const Z_type expectedZ = std::exp(2.0) + std::exp(3.0);
		REQUIRE(predictor.getZall() == Approx(expectedZ).epsilon(1e-12));
	}

	SECTION("seed-extension predictor reuse clears an empty range partition") {
		RnaSequence target("target", "GGG");
		RnaSequence query("query", "CCC");
		AccessibilityDisabled targetAcc(target, 0, NULL);
		AccessibilityDisabled queryAcc(query, 0, NULL);
		ReverseAccessibility reverseQueryAcc(queryAcc);
		InteractionEnergyBasePair energy(targetAcc, reverseQueryAcc);
		SeedConstraint seedConstraint(2, 0, 0, 0,
				E_INF, Accessibility::ED_UPPER_BOUND, E_INF,
				IndexRangeList(), IndexRangeList(), "", false, false, true);
		OutputConstraint constraint(1, OutputConstraint::OVERLAP_BOTH,
				E_INF, E_INF, false, true, false, true, false);

		OutputHandlerInteractionList exactOut(constraint, 1);
		InspectableEnsemblePredictor<PredictorMfeEns2dSeedExtension> exact(
				energy, exactOut, new SeedHandlerNoBulge(energy, seedConstraint));
		exact.predict();
		REQUIRE(exact.getPartitionCount() > 0);
		exact.predict(IndexRange(0, 0), IndexRange(0, 0));
		REQUIRE(exact.getPartitionCount() == 0);
		REQUIRE(exact.getZall() == 0.0);

		OutputHandlerInteractionList heuristicOut(constraint, 1);
		InspectableEnsemblePredictor<PredictorMfeEns2dHeuristicSeedExtension> heuristic(
				energy, heuristicOut, new SeedHandlerNoBulge(energy, seedConstraint));
		heuristic.predict();
		REQUIRE(heuristic.getPartitionCount() > 0);
		heuristic.predict(IndexRange(0, 0), IndexRange(0, 0));
		REQUIRE(heuristic.getPartitionCount() == 0);
		REQUIRE(heuristic.getZall() == 0.0);
	}
}
