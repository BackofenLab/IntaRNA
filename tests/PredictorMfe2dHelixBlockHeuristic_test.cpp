
#include "catch.hpp"

#undef NDEBUG
#define protected public

#include "IntaRNA/RnaSequence.h"
#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/HelixConstraint.h"
#include "IntaRNA/HelixHandlerNoBulgeMax.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/Interaction.h"
#include "IntaRNA/ReverseAccessibility.h"
#include "IntaRNA/PredictorMfe2dHelixBlockHeuristic.h"
#include "IntaRNA/OutputHandlerInteractionList.h"

using namespace IntaRNA;

TEST_CASE( "PredictorMfe2dHelixBlockHeuristic", "[PredictorMfe2dHelixBlockHeuristic]") {

	SECTION("Predictor: Case 1", "[PredictorMfe2dHelixBlockHeuristic]") {

		RnaSequence r1("r1", "GGGAAGG");
		RnaSequence r2("r2", "CCAACCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hc(2, 4,  0, 999, 0, false);

		OutputConstraint outC(1,OutputConstraint::OVERLAP_SEQ2,0,100);
		OutputHandlerInteractionList out(outC, 1);

		PredictorMfe2dHelixBlockHeuristic pLSH(energy, out, NULL, hc);

		IndexRange idx1(0,r1.lastPos);
		IndexRange idx2(0,r2.lastPos);

		pLSH.predict(idx1,idx2);

		REQUIRE_FALSE(out.empty());
		REQUIRE(out.reported() == 1);

		const Interaction * interaction((*out.begin()));
		REQUIRE(interaction->basePairs.begin()->first == 0);
		REQUIRE(interaction->basePairs.begin()->second == 6);

		REQUIRE(interaction->basePairs.rbegin()->first == 6);
		REQUIRE(interaction->basePairs.rbegin()->second == 0);

		REQUIRE(interaction->dotBracket(*interaction) == "(((..((&))..)))");
	}

	SECTION("Predictor: Case 2", "[PredictorMfe2dHelixBlockHeuristic]") {

		RnaSequence r1("r1", "AAGGGAAGGA");
		RnaSequence r2("r2", "ACCAACCCAA");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hc(2, 4,  0, 999, 0, false);

		OutputConstraint outC(1,OutputConstraint::OVERLAP_SEQ2,0,100);
		OutputHandlerInteractionList out(outC,1);

		PredictorMfe2dHelixBlockHeuristic pLSH(energy, out, NULL, hc);

		IndexRange idx1(0,r1.lastPos);
		IndexRange idx2(0,r2.lastPos);

		pLSH.predict(idx1,idx2);

		REQUIRE_FALSE(out.empty());
		REQUIRE(out.reported() == 1);

		const Interaction * interaction((*out.begin()));
		REQUIRE(interaction->basePairs.begin()->first == 2);
		REQUIRE(interaction->basePairs.begin()->second == 7);

		REQUIRE(interaction->basePairs.rbegin()->first == 8);
		REQUIRE(interaction->basePairs.rbegin()->second == 1);

		REQUIRE(interaction->dotBracket(*interaction) == "(((..((&))..)))");
	}

	SECTION("Predictor: Case 3 - no favorable output", "[PredictorMfe2dHelixBlockHeuristic]") {

		RnaSequence r1("r1", "AAAAAAAA");
		RnaSequence r2("r2", "AAAAAAAA");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hc(2, 4,  0, 999, 0, false);

		OutputConstraint outC(1,OutputConstraint::OVERLAP_SEQ2
				,Ekcal_2_E(0) // maxE
				,Ekcal_2_E(1.0) // deltaE
				);
		OutputHandlerInteractionList out(outC,1);

		PredictorMfe2dHelixBlockHeuristic pLSH(energy, out, NULL, hc);

		IndexRange idx1(0,r1.lastPos);
		IndexRange idx2(0,r2.lastPos);

		pLSH.predict(idx1,idx2);

		REQUIRE(out.empty());

	}

	SECTION("Predictor: Case 4", "[PredictorMfe2dHelixBlockHeuristic]") {

		RnaSequence r1("r1", "AAAAAAAAGGG");
		RnaSequence r2("r2", "AAAAAACCCAA");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hc(2, 4,  0, 999, 0, false);

		OutputConstraint outC(1,OutputConstraint::OVERLAP_SEQ2
				,Ekcal_2_E(0) // maxE
				,Ekcal_2_E(1.0) // deltaE
				);
		OutputHandlerInteractionList out(outC,1);

		PredictorMfe2dHelixBlockHeuristic pLSH(energy, out, NULL, hc);

		IndexRange idx1(0,r1.lastPos);
		IndexRange idx2(0,r2.lastPos);

		pLSH.predict(idx1,idx2);

		REQUIRE_FALSE(out.empty());
		REQUIRE(out.reported() == 1);

		const Interaction * interaction((*out.begin()));
		REQUIRE(interaction->basePairs.begin()->first == 8);
		REQUIRE(interaction->basePairs.begin()->second == 8);

		REQUIRE(interaction->basePairs.rbegin()->first == 10);
		REQUIRE(interaction->basePairs.rbegin()->second == 6);

		REQUIRE(interaction->dotBracket(*interaction) == "(((&)))");
	}

	SECTION("Predictor: Case 5 - unpaired bases allowed", "[PredictorMfe2dHelixBlockHeuristic]") {

		RnaSequence r1("r1", "GGGAGGGG");
		RnaSequence r2("r2", "CCCCACCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hc(2, 4, 2, 999, 0, false);

		OutputConstraint outC(1
				,OutputConstraint::OVERLAP_SEQ2
				,Ekcal_2_E(0) // maxE
				,Ekcal_2_E(1.0) // deltaE
				);
		OutputHandlerInteractionList out(outC,1);

		PredictorMfe2dHelixBlockHeuristic pLSH(energy, out, NULL, hc);

		IndexRange idx1(0,r1.lastPos);
		IndexRange idx2(0,r2.lastPos);

		pLSH.predict(idx1,idx2);

		REQUIRE_FALSE(out.empty());
		REQUIRE(out.reported() == 1);

		const Interaction * interaction((*out.begin()));

		REQUIRE(interaction->dotBracket(*interaction) == "(((.(..(&).)).))");

		REQUIRE(interaction->basePairs.begin()->first == 0);
		REQUIRE(interaction->basePairs.begin()->second == 6);

		REQUIRE(interaction->basePairs.rbegin()->first == 7);
		REQUIRE(interaction->basePairs.rbegin()->second == 0);

	}


	SECTION("Predictor: Case 6 - unpaired bases allowed", "[PredictorMfe2dHelixBlockHeuristic]") {

		RnaSequence r1("r1", "GGGAAGG");
		RnaSequence r2("r2", "CCCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hc(2, 5, 2, 999, 0, false);

		OutputConstraint outC(1,OutputConstraint::OVERLAP_SEQ2
				,Ekcal_2_E(0) // maxE
				,Ekcal_2_E(1.0) // deltaE
				);
		OutputHandlerInteractionList out(outC,1);

		PredictorMfe2dHelixBlockHeuristic pLSH(energy, out, NULL, hc);

		IndexRange idx1(0,r1.lastPos);
		IndexRange idx2(0,r2.lastPos);

		pLSH.predict(idx1,idx2);

		REQUIRE_FALSE(out.empty());
		REQUIRE(out.reported() == 1);

		const Interaction * interaction((*out.begin()));
		REQUIRE(interaction->basePairs.begin()->first == 0);
		REQUIRE(interaction->basePairs.begin()->second == 4);

		REQUIRE(interaction->basePairs.rbegin()->first == 6);
		REQUIRE(interaction->basePairs.rbegin()->second == 0);

		REQUIRE(interaction->dotBracket(*interaction) == "(((..((&)))))");
	}

	SECTION("Predictor: Case 7 - only stacking allowed", "[PredictorMfe2dHelixBlockHeuristic]") {

		RnaSequence r1("r1", "GGGGGG");
		RnaSequence r2("r2", "CCCCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hc(2, 4,  0, 999, 0, false);

		OutputConstraint outC(1,OutputConstraint::OVERLAP_SEQ2
				,Ekcal_2_E(0) // maxE
				,Ekcal_2_E(1.0) // deltaE
				);
		OutputHandlerInteractionList out(outC,1);

		PredictorMfe2dHelixBlockHeuristic pLSH(energy, out, NULL, hc);

		IndexRange idx1(0,r1.lastPos);
		IndexRange idx2(0,r2.lastPos);

		pLSH.predict(idx1,idx2);

		REQUIRE_FALSE(out.empty());
		REQUIRE(out.reported() == 1);

		const Interaction * interaction((*out.begin()));
		REQUIRE(interaction->dotBracket(*interaction) == "((((.(&)))))");
		REQUIRE(interaction->basePairs.begin()->first == 0);
		REQUIRE(interaction->basePairs.begin()->second == 4);

		REQUIRE(interaction->basePairs.rbegin()->first == 5);
		REQUIRE(interaction->basePairs.rbegin()->second == 0);

	}

}
