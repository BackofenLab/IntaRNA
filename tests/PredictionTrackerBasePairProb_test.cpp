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
#include "IntaRNA/SeedHandlerNoBulge.h"

#include <stdexcept>

using namespace IntaRNA;

//! Helper function to generate Boltzmann weights for different numbers of base pairs
inline std::vector<Z_type> getBPWeights(const InteractionEnergyBasePair & energy, const size_t maxBasepairs)
{
	std::vector<Z_type> bpWeights;
	for (size_t i = 0; i <= maxBasepairs; i++) {
		bpWeights.push_back(energy.getBoltzmannWeight(i * energy.getE_basePair()));
	}
	return bpWeights;
}

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

		std::vector<Z_type> wBP = getBPWeights(energy, 2);

		REQUIRE(Z_equal(tracker->getBasePairProb(0, 1, &predictor), wBP[1] / predictor.getZall()));
		REQUIRE(Z_equal(tracker->getBasePairProb(1, 0, &predictor), wBP[1] / predictor.getZall()));
		REQUIRE(Z_equal(tracker->getBasePairProb(0, 0, &predictor), (wBP[1] + wBP[2]) / predictor.getZall()));
		REQUIRE(Z_equal(tracker->getBasePairProb(1, 1, &predictor), (wBP[1] + wBP[2]) / predictor.getZall()));
	}

	SECTION("base pair probs - case 1 seed") {
		RnaSequence r1("r1", "GG");
		RnaSequence r2("r2", "CC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		OutputConstraint outC(1,OutputConstraint::OVERLAP_SEQ2,0,10000,0,0,0,true,false);
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

		PredictorMfeEns2dSeedExtension predictor(energy, out, tracker, new SeedHandlerNoBulge(energy, sC));
		predictor.predict(idx1,idx2);

		REQUIRE(Z_equal(tracker->getBasePairProb(0, 1, &predictor), 0));
		REQUIRE(Z_equal(tracker->getBasePairProb(1, 0, &predictor), 0));
		REQUIRE(Z_equal(tracker->getBasePairProb(0, 0, &predictor), 1));
		REQUIRE(Z_equal(tracker->getBasePairProb(1, 1, &predictor), 1));
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

		std::vector<Z_type> wBP = getBPWeights(energy, 3);

		REQUIRE(Z_equal(tracker->getBasePairProb(0, 1, &predictor), (wBP[1] + 2 * wBP[2]) / predictor.getZall()));
		REQUIRE(Z_equal(tracker->getBasePairProb(1, 0, &predictor), (wBP[1] + 2 * wBP[2]) / predictor.getZall()));
		REQUIRE(Z_equal(tracker->getBasePairProb(2, 1, &predictor), (wBP[1] + 2 * wBP[2]) / predictor.getZall()));
		REQUIRE(Z_equal(tracker->getBasePairProb(1, 2, &predictor), (wBP[1] + 2 * wBP[2]) / predictor.getZall()));

		REQUIRE(Z_equal(tracker->getBasePairProb(0, 2, &predictor), wBP[1] / predictor.getZall()));
		REQUIRE(Z_equal(tracker->getBasePairProb(2, 0, &predictor), wBP[1] / predictor.getZall()));

		REQUIRE(Z_equal(tracker->getBasePairProb(0, 0, &predictor), (wBP[1] + 4 * wBP[2] + wBP[3]) / predictor.getZall()));
		REQUIRE(Z_equal(tracker->getBasePairProb(2, 2, &predictor), (wBP[1] + 4 * wBP[2] + wBP[3]) / predictor.getZall()));
		REQUIRE(Z_equal(tracker->getBasePairProb(1, 1, &predictor), (wBP[1] + 2 * wBP[2] + wBP[3]) / predictor.getZall()));
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

		PredictorMfeEns2dSeedExtension predictor(energy, out, tracker, new SeedHandlerNoBulge(energy, sC));
		predictor.predict(idx1,idx2);

		REQUIRE(Z_equal(tracker->getBasePairProb(0, 1, &predictor), 0));
		REQUIRE(Z_equal(tracker->getBasePairProb(1, 0, &predictor), 0));
		REQUIRE(Z_equal(tracker->getBasePairProb(2, 1, &predictor), 0));
		REQUIRE(Z_equal(tracker->getBasePairProb(1, 2, &predictor), 0));

		REQUIRE(Z_equal(tracker->getBasePairProb(0, 2, &predictor), 0));
		REQUIRE(Z_equal(tracker->getBasePairProb(2, 0, &predictor), 0));

		REQUIRE(Z_equal(tracker->getBasePairProb(0, 0, &predictor), 1));
		REQUIRE(Z_equal(tracker->getBasePairProb(2, 2, &predictor), 1));
		REQUIRE(Z_equal(tracker->getBasePairProb(1, 1, &predictor), 1));
	}

	SECTION("base pair probs - case 3 noseed") {
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

		std::vector<Z_type> wBP = getBPWeights(energy, 4);

		REQUIRE(Z_equal(tracker->getBasePairProb(1, 1, &predictor), (wBP[1] + 5 * wBP[2] + 5 * wBP[3] + wBP[4]) / predictor.getZall()));
	}

	SECTION("base pair probs - case 4 seed") {
		LOG(DEBUG)<<"\n\n test case 4\n";
		RnaSequence r1("r1", "GGGG");
		RnaSequence r2("r2", "CCCC");
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


		PredictorMfeEns2dSeedExtension predictor(energy, out, tracker, new SeedHandlerNoBulge(energy, sC));
		predictor.predict(idx1,idx2);

		std::vector<Z_type> wBP = getBPWeights(energy, 4);
		LOG(DEBUG) <<" w 1 "<<wBP[1] <<" 2 "<<wBP[2]<<" 3 "<<wBP[3] <<" 4 "<<wBP[4];

		REQUIRE(Z_equal(tracker->getBasePairProb(0, 0, &predictor), (wBP[3] + wBP[4]) / predictor.getZall()));
		LOG(DEBUG) <<"Z(1,2)"<<tracker->getBasePairProb(1, 2, &predictor)*predictor.getZall() <<"  Zall "<<predictor.getZall();
		REQUIRE(Z_equal(tracker->getBasePairProb(1, 1, &predictor), (2 * wBP[3] + wBP[4]) / predictor.getZall()));
		REQUIRE(Z_equal(tracker->getBasePairProb(2, 2, &predictor), (2 * wBP[3] + wBP[4]) / predictor.getZall()));
		REQUIRE(Z_equal(tracker->getBasePairProb(3, 3, &predictor), (wBP[3] + wBP[4]) / predictor.getZall()));

		REQUIRE(Z_equal(tracker->getBasePairProb(0, 1, &predictor), wBP[3] / predictor.getZall()));
		REQUIRE(Z_equal(tracker->getBasePairProb(1, 2, &predictor), wBP[3] / predictor.getZall()));
		REQUIRE(Z_equal(tracker->getBasePairProb(2, 3, &predictor), wBP[3] / predictor.getZall()));

		REQUIRE(Z_equal(tracker->getBasePairProb(3, 2, &predictor), wBP[3] / predictor.getZall()));
		REQUIRE(Z_equal(tracker->getBasePairProb(2, 1, &predictor), wBP[3] / predictor.getZall()));
		REQUIRE(Z_equal(tracker->getBasePairProb(1, 0, &predictor), wBP[3] / predictor.getZall()));
	}

	SECTION("base pair probs - case 5 seed") {
		RnaSequence r1("r1", "GGGG");
		RnaSequence r2("r2", "CCCC");
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

		PredictorMfeEns2dSeedExtension predictor(energy, out, tracker, new SeedHandlerNoBulge(energy, sC));
		predictor.predict(idx1,idx2);

		std::vector<Z_type> wBP = getBPWeights(energy, 4);

		REQUIRE(Z_equal(tracker->getBasePairProb(0, 0, &predictor), (1 * wBP[2] + 7 * wBP[3] + wBP[4]) / predictor.getZall()));
		REQUIRE(Z_equal(tracker->getBasePairProb(1, 1, &predictor), (2 * wBP[2] + 5 * wBP[3] + wBP[4]) / predictor.getZall()));
		REQUIRE(Z_equal(tracker->getBasePairProb(2, 2, &predictor), (2 * wBP[2] + 5 * wBP[3] + wBP[4]) / predictor.getZall()));
		REQUIRE(Z_equal(tracker->getBasePairProb(3, 3, &predictor), (1 * wBP[2] + 7 * wBP[3] + wBP[4]) / predictor.getZall()));
	}
//
//	SECTION("base pair probs - case 6 seed") {
//		RnaSequence r1("r1", "GGGGGGG");
//		RnaSequence r2("r2", "CCCCCCC");
//		AccessibilityDisabled acc1(r1, 0, NULL);
//		AccessibilityDisabled acc2(r2, 0, NULL);
//		ReverseAccessibility racc(acc2);
//		InteractionEnergyBasePair energy(acc1, racc);
//
//		OutputConstraint outC(1,OutputConstraint::OVERLAP_SEQ2,0,10000,0,0,0,1,1);
//		OutputHandlerInteractionList out(outC, 1);
//
//		IndexRange idx1(0,r1.lastPos);
//		IndexRange idx2(0,r2.lastPos);
//
//		PredictionTrackerBasePairProb * tracker = new PredictionTrackerBasePairProb(energy, "");
//		SeedConstraint sC(5,0,0,0,0
//				, AccessibilityDisabled::ED_UPPER_BOUND
//				, 0
//				, IndexRangeList("")
//				, IndexRangeList("")
//				, ""
//				, false, false
//				);
//
//		PredictorMfeEns2dSeedExtension predictor(energy, out, tracker, new SeedHandlerNoBulge(energy, sC));
//		predictor.predict(idx1,idx2);
//
//		std::vector<Z_type> wBP = getBPWeights(energy, 7);
//
//		REQUIRE(Z_equal(tracker->getBasePairProb(3, 3, &predictor), (3 * wBP[5] + 6 * wBP[6] + wBP[7]) / predictor.getZall()));
//	}
//
	// SECTION("base pair probs - case 7 seed") {
	// 	RnaSequence r1("r1", "GGGGGG");
	// 	RnaSequence r2("r2", "CCCCCC");
	// 	AccessibilityDisabled acc1(r1, 0, NULL);
	// 	AccessibilityDisabled acc2(r2, 0, NULL);
	// 	ReverseAccessibility racc(acc2);
	// 	InteractionEnergyBasePair energy(acc1, racc);
	//
	// 	OutputConstraint outC(1,OutputConstraint::OVERLAP_SEQ2,0,10000,0,0,0,1,1);
	// 	OutputHandlerInteractionList out(outC, 1);
	//
	// 	IndexRange idx1(0,r1.lastPos);
	// 	IndexRange idx2(0,r2.lastPos);
	//
	// 	PredictionTrackerBasePairProb * tracker = new PredictionTrackerBasePairProb(energy, "");
	// 	SeedConstraint sC(4,0,0,0,0
	// 			, AccessibilityDisabled::ED_UPPER_BOUND
	// 			, 0
	// 			, IndexRangeList("")
	// 			, IndexRangeList("")
	// 			, ""
	// 			, false, false
	// 			);
	//
	// 	PredictorMfeEns2dSeedExtension predictor(energy, out, tracker, new SeedHandlerNoBulge(energy, sC));
	// 	predictor.predict(idx1,idx2);
	//
	// 	std::vector<Z_type> wBP = getBPWeights(energy, 6);
	//
	// 	REQUIRE(Z_equal(tracker->getBasePairProb(2, 2, &predictor), (3 * wBP[4] + 6 * wBP[5] + wBP[6]) / predictor.getZall()));
	// 	REQUIRE(Z_equal(tracker->getBasePairProb(3, 3, &predictor), (3 * wBP[4] + 6 * wBP[5] + wBP[6]) / predictor.getZall()));
	// }

}
