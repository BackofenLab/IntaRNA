#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/OutputHandlerInteractionList.h"
#include "IntaRNA/PredictorMfeEns2d.h"
#include "IntaRNA/PredictorMfeEns2dHeuristic.h"
#include "IntaRNA/ReverseAccessibility.h"
#include "IntaRNA/RnaSequence.h"

using namespace IntaRNA;

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
}
