#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/RnaSequence.h"
#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/Interaction.h"
#include "IntaRNA/ReverseAccessibility.h"
#include "IntaRNA/PredictorMfeEns2d.h"
#include "IntaRNA/OutputHandlerInteractionList.h"

using namespace IntaRNA;

#include <sstream>

TEST_CASE( "PredictorMfeEns2d", "[PredictorMfeEns2d]" ) {

	// setup easylogging++ stuff if not already done
	#include "testEasyLoggingSetup.icc"

	SECTION("Zall case 1: check value") {
		RnaSequence r1("r1", "GG");
		RnaSequence r2("r2", "CC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

    OutputConstraint outC(1,OutputConstraint::OVERLAP_SEQ2,0,100);
		OutputHandlerInteractionList out(outC, 1);

		PredictorMfeEns2d predictor(energy, out, NULL);

		IndexRange idx1(0,r1.lastPos);
		IndexRange idx2(0,r2.lastPos);

		predictor.predict(idx1,idx2);

		Z_type boltzmannSum = 4 * energy.getBoltzmannWeight(Ekcal_2_E(-1.0))
		                    + 1 * energy.getBoltzmannWeight(Ekcal_2_E(-2.0));

		REQUIRE(Z_equal(predictor.getZall(), boltzmannSum));
	}

	SECTION("Zall case 2: check value") {
		RnaSequence r1("r1", "GGG");
		RnaSequence r2("r2", "CCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

    OutputConstraint outC(1,OutputConstraint::OVERLAP_SEQ2,0,100);
		OutputHandlerInteractionList out(outC, 1);

		PredictorMfeEns2d predictor(energy, out, NULL);

		IndexRange idx1(0,r1.lastPos);
		IndexRange idx2(0,r2.lastPos);

		predictor.predict(idx1,idx2);

		Z_type boltzmannSum = 9 * energy.getBoltzmannWeight(Ekcal_2_E(-1.0))
		                    + 9 * energy.getBoltzmannWeight(Ekcal_2_E(-2.0))
								        + 1 * energy.getBoltzmannWeight(Ekcal_2_E(-3.0));

		REQUIRE(Z_equal(predictor.getZall(), boltzmannSum));
	}

	SECTION("Zall case 3: check value") {
		RnaSequence r1("r1", "GGGG");
		RnaSequence r2("r2", "CCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

    OutputConstraint outC(1,OutputConstraint::OVERLAP_SEQ2,0,100);
		OutputHandlerInteractionList out(outC, 1);

		PredictorMfeEns2d predictor(energy, out, NULL);

		IndexRange idx1(0,r1.lastPos);
		IndexRange idx2(0,r2.lastPos);

		predictor.predict(idx1,idx2);

		Z_type boltzmannSum = 16 * energy.getBoltzmannWeight(Ekcal_2_E(-1.0))
		                    + 36 * energy.getBoltzmannWeight(Ekcal_2_E(-2.0))
								        + 16 * energy.getBoltzmannWeight(Ekcal_2_E(-3.0))
												+ 1 * energy.getBoltzmannWeight(Ekcal_2_E(-4.0));

		REQUIRE(Z_equal(predictor.getZall(), boltzmannSum));
	}

	SECTION("Zall case 4: check value") {
		RnaSequence r1("r1", "GGC");
		RnaSequence r2("r2", "GCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

    OutputConstraint outC(1,OutputConstraint::OVERLAP_SEQ2,0,100);
		OutputHandlerInteractionList out(outC, 1);

		PredictorMfeEns2d predictor(energy, out, NULL);

		IndexRange idx1(0,r1.lastPos);
		IndexRange idx2(0,r2.lastPos);

		predictor.predict(idx1,idx2);

		Z_type boltzmannSum = 5 * energy.getBoltzmannWeight(Ekcal_2_E(-1.0))
		                    + 5 * energy.getBoltzmannWeight(Ekcal_2_E(-2.0))
								        + 1 * energy.getBoltzmannWeight(Ekcal_2_E(-3.0));

		REQUIRE(Z_equal(predictor.getZall(), boltzmannSum));
	}

}
