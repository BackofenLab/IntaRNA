#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/PredictionTrackerBasePairProb.h"
#include "IntaRNA/RnaSequence.h"
#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/Interaction.h"
#include "IntaRNA/ReverseAccessibility.h"
#include "IntaRNA/PredictorMfeEns2d.h"
#include "IntaRNA/PredictorMfeEns2dSeedExtension.h"
#include "IntaRNA/OutputHandlerInteractionList.h"
#include "IntaRNA/SeedHandlerMfe.h"

#include <stdexcept>

using namespace IntaRNA;

TEST_CASE( "PredictionTrackerBasePairProb", "[PredictionTrackerBasePairProb]" ) {

	// setup easylogging++ stuff if not already done
	#include "testEasyLoggingSetup.icc"

	SECTION("base pair probs - case 1 noseed") {
		RnaSequence r1("r1", "GG");
		RnaSequence r2("r2", "CC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		OutputConstraint outC(1,OutputConstraint::OVERLAP_SEQ2,0,100);
		OutputHandlerInteractionList out(outC, 1);

		IndexRange idx1(0,r1.lastPos);
		IndexRange idx2(0,r2.lastPos);

		PredictionTrackerBasePairProb * tracker = new PredictionTrackerBasePairProb(energy, "");
		PredictorMfeEns2d predictor(energy, out, tracker);

		predictor.predict(idx1,idx2);

		REQUIRE(Z_equal(tracker->getBasePairProb(0, 0, 1, 1, &predictor), energy.getBoltzmannWeight(Ekcal_2_E(-1.0)) / predictor.getZall()));
		REQUIRE(Z_equal(tracker->getBasePairProb(1, 1, 0, 0, &predictor), energy.getBoltzmannWeight(Ekcal_2_E(-1.0)) / predictor.getZall()));
		REQUIRE(Z_equal(tracker->getBasePairProb(0, 0, 0, 0, &predictor), (energy.getBoltzmannWeight(Ekcal_2_E(-1.0)) + energy.getBoltzmannWeight(Ekcal_2_E(-2.0))) / predictor.getZall()));
		REQUIRE(Z_equal(tracker->getBasePairProb(1, 1, 1, 1, &predictor), (energy.getBoltzmannWeight(Ekcal_2_E(-1.0)) + energy.getBoltzmannWeight(Ekcal_2_E(-2.0))) / predictor.getZall()));
	}

	SECTION("base pair probs - case 1 seed") {
		RnaSequence r1("r1", "GG");
		RnaSequence r2("r2", "CC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		OutputConstraint outC(1,OutputConstraint::OVERLAP_SEQ2,0,10000,0,0,0,1,1);
		OutputHandlerInteractionList out(outC, 1);

		IndexRange idx1(0,r1.lastPos);
		IndexRange idx2(0,r2.lastPos);

		PredictionTrackerBasePairProb * tracker = new PredictionTrackerBasePairProb(energy, "");
		SeedConstraint sC(2,0,0,0,0
				, AccessibilityDisabled::ED_UPPER_BOUND
				, 0
				, IndexRangeList("")
				, IndexRangeList("")
				, ""
				, false, false
				);

		PredictorMfeEns2dSeedExtension predictor(energy, out, tracker, new SeedHandlerMfe(energy, sC));
		predictor.predict(idx1,idx2);

		REQUIRE(Z_equal(tracker->getBasePairProb(0, 0, 1, 1, &predictor), 0));
		REQUIRE(Z_equal(tracker->getBasePairProb(1, 1, 0, 0, &predictor), 0));
		REQUIRE(Z_equal(tracker->getBasePairProb(0, 0, 0, 0, &predictor), 1));
		REQUIRE(Z_equal(tracker->getBasePairProb(1, 1, 1, 1, &predictor), 1));
	}

	SECTION("base pair probs - case 2 noseed") {
		RnaSequence r1("r1", "GGG");
		RnaSequence r2("r2", "CCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		OutputConstraint outC(1,OutputConstraint::OVERLAP_SEQ2,0,100);
		OutputHandlerInteractionList out(outC, 1);

		IndexRange idx1(0,r1.lastPos);
		IndexRange idx2(0,r2.lastPos);

		PredictionTrackerBasePairProb * tracker = new PredictionTrackerBasePairProb(energy, "");
		PredictorMfeEns2d predictor(energy, out, tracker);

		predictor.predict(idx1,idx2);

		REQUIRE(Z_equal(tracker->getBasePairProb(0, 0, 1, 1, &predictor), (energy.getBoltzmannWeight(Ekcal_2_E(-1.0)) + 2 * energy.getBoltzmannWeight(Ekcal_2_E(-2.0))) / predictor.getZall()));
		REQUIRE(Z_equal(tracker->getBasePairProb(1, 1, 0, 0, &predictor), (energy.getBoltzmannWeight(Ekcal_2_E(-1.0)) + 2 * energy.getBoltzmannWeight(Ekcal_2_E(-2.0))) / predictor.getZall()));
		REQUIRE(Z_equal(tracker->getBasePairProb(2, 2, 1, 1, &predictor), (energy.getBoltzmannWeight(Ekcal_2_E(-1.0)) + 2 * energy.getBoltzmannWeight(Ekcal_2_E(-2.0))) / predictor.getZall()));
		REQUIRE(Z_equal(tracker->getBasePairProb(1, 1, 2, 2, &predictor), (energy.getBoltzmannWeight(Ekcal_2_E(-1.0)) + 2 * energy.getBoltzmannWeight(Ekcal_2_E(-2.0))) / predictor.getZall()));

		REQUIRE(Z_equal(tracker->getBasePairProb(0, 0, 2, 2, &predictor), energy.getBoltzmannWeight(Ekcal_2_E(-1.0)) / predictor.getZall()));
		REQUIRE(Z_equal(tracker->getBasePairProb(2, 2, 0, 0, &predictor), energy.getBoltzmannWeight(Ekcal_2_E(-1.0)) / predictor.getZall()));

		REQUIRE(Z_equal(tracker->getBasePairProb(0, 0, 0, 0, &predictor), (energy.getBoltzmannWeight(Ekcal_2_E(-1.0)) + 4 * energy.getBoltzmannWeight(Ekcal_2_E(-2.0)) + energy.getBoltzmannWeight(Ekcal_2_E(-3.0))) / predictor.getZall()));
		REQUIRE(Z_equal(tracker->getBasePairProb(2, 2, 2, 2, &predictor), (energy.getBoltzmannWeight(Ekcal_2_E(-1.0)) + 4 * energy.getBoltzmannWeight(Ekcal_2_E(-2.0)) + energy.getBoltzmannWeight(Ekcal_2_E(-3.0))) / predictor.getZall()));
		REQUIRE(Z_equal(tracker->getBasePairProb(1, 1, 1, 1, &predictor), (energy.getBoltzmannWeight(Ekcal_2_E(-1.0)) + 2 * energy.getBoltzmannWeight(Ekcal_2_E(-2.0)) + energy.getBoltzmannWeight(Ekcal_2_E(-3.0))) / predictor.getZall()));
	}

	SECTION("base pair probs - case 2 seed") {
		RnaSequence r1("r1", "GGG");
		RnaSequence r2("r2", "CCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		OutputConstraint outC(1,OutputConstraint::OVERLAP_SEQ2,0,10000,0,0,0,1,1);
		OutputHandlerInteractionList out(outC, 1);

		IndexRange idx1(0,r1.lastPos);
		IndexRange idx2(0,r2.lastPos);

		PredictionTrackerBasePairProb * tracker = new PredictionTrackerBasePairProb(energy, "");
		SeedConstraint sC(3,0,0,0,0
				, AccessibilityDisabled::ED_UPPER_BOUND
				, 0
				, IndexRangeList("")
				, IndexRangeList("")
				, ""
				, false, false
				);

		PredictorMfeEns2dSeedExtension predictor(energy, out, tracker, new SeedHandlerMfe(energy, sC));
		predictor.predict(idx1,idx2);

		REQUIRE(Z_equal(tracker->getBasePairProb(0, 0, 1, 1, &predictor), 0));
		REQUIRE(Z_equal(tracker->getBasePairProb(1, 1, 0, 0, &predictor), 0));
		REQUIRE(Z_equal(tracker->getBasePairProb(2, 2, 1, 1, &predictor), 0));
		REQUIRE(Z_equal(tracker->getBasePairProb(1, 1, 2, 2, &predictor), 0));

		REQUIRE(Z_equal(tracker->getBasePairProb(0, 0, 2, 2, &predictor), 0));
		REQUIRE(Z_equal(tracker->getBasePairProb(2, 2, 0, 0, &predictor), 0));

		REQUIRE(Z_equal(tracker->getBasePairProb(0, 0, 0, 0, &predictor), 1));
		REQUIRE(Z_equal(tracker->getBasePairProb(2, 2, 2, 2, &predictor), 1));
		REQUIRE(Z_equal(tracker->getBasePairProb(1, 1, 1, 1, &predictor), 1));
	}

	SECTION("base pair probs - case 3") {
		RnaSequence r1("r1", "GGCGC");
		RnaSequence r2("r2", "GGCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		OutputConstraint outC(1,OutputConstraint::OVERLAP_SEQ2,0,100);
		OutputHandlerInteractionList out(outC, 1);

		IndexRange idx1(0,r1.lastPos);
		IndexRange idx2(0,r2.lastPos);

		PredictionTrackerBasePairProb * tracker = new PredictionTrackerBasePairProb(energy, "");
		PredictorMfeEns2d predictor(energy, out, tracker);

    predictor.predict(idx1,idx2);

		REQUIRE(Z_equal(tracker->getBasePairProb(1, 1, 1, 1, &predictor), (energy.getBoltzmannWeight(Ekcal_2_E(-1.0)) + 5 * energy.getBoltzmannWeight(Ekcal_2_E(-2.0)) + 5 * energy.getBoltzmannWeight(Ekcal_2_E(-3.0)) + energy.getBoltzmannWeight(Ekcal_2_E(-4.0))) / predictor.getZall()));
	}

}
