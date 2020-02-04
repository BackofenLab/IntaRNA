
#include "catch.hpp"

#undef NDEBUG
#define protected public

#include "IntaRNA/RnaSequence.h"
#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/HelixConstraint.h"
#include "IntaRNA/HelixHandlerNoBulgeMax.h"
#include "IntaRNA/HelixHandlerIdxOffset.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/SeedHandlerMfe.h"


using namespace IntaRNA;

TEST_CASE( "HelixHandlerIdxOffset for NoBulgeMaxSeed", "[HelixHandlerIdxOffset]") {

	SECTION("HelixSeed with Offset: Case 1 - offset 1", "[HelixHandlerIdxOffset]") {

		RnaSequence r1("r1", "GGGGGG");
		RnaSequence r2("r2", "CCCCCG");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4,  0, 999, 0, false);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ / seedNoGU
		SeedConstraint sC(3, 0, 0, 0, 0, AccessibilityDisabled::ED_UPPER_BOUND, 0, IndexRangeList(""), IndexRangeList(""),
						  "", false, false, false);

		HelixHandler *hhS = new HelixHandlerNoBulgeMax(energy, hC);
		SeedHandler *sH = new SeedHandlerMfe(energy, sC);

		hhS->setSeedHandler(*sH);

		SeedHandlerIdxOffset sHIO(sH);
		HelixHandlerIdxOffset hhIO(hhS);

		// Set the offsets
		sHIO.setOffset1(1);
		sHIO.setOffset2(1);
		hhIO.setOffset1(1);
		hhIO.setOffset2(1);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////   FILLHELIXSEED  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		hhIO.fillHelix(0, energy.size1() - sHIO.getOffset1() - 1, 0, energy.size2() - sHIO.getOffset2() - 1);
		REQUIRE(sHIO.fillSeed(0, energy.size1() - sHIO.getOffset1() - 1, 0, energy.size2() - sHIO.getOffset2() - 1) ==
				9);
		REQUIRE(hhIO.fillHelixSeed(0, energy.size1() - hhIO.getOffset1() - 1, 0,
								   energy.size2() - hhIO.getOffset2() - 1) == 9);

		// (0,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(0, 0)) == -3);
		REQUIRE(hhIO.getHelixSeedLength1(0, 0) == 4);
		REQUIRE(hhIO.getHelixSeedLength2(0, 0) == 4);

		// (0,1)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(0, 1)) == -3);
		REQUIRE(hhIO.getHelixSeedLength1(0, 1) == 4);
		REQUIRE(hhIO.getHelixSeedLength2(0, 1) == 4);

		// (0,2)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(0, 2)) == -2);
		REQUIRE(hhIO.getHelixSeedLength1(0, 2) == 3);
		REQUIRE(hhIO.getHelixSeedLength2(0, 2) == 3);

		// (0,3)
		REQUIRE(hhIO.getHelixSeedE(0, 3) == E_INF);
		REQUIRE(hhIO.getHelixSeedLength1(0, 3) == 0);
		REQUIRE(hhIO.getHelixSeedLength2(0, 3) == 0);

		// (1,1)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(1, 1)) == -3);
		REQUIRE(hhIO.getHelixSeedLength1(1, 1) == 4);
		REQUIRE(hhIO.getHelixSeedLength2(1, 1) == 4);

		// (2,2)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(2, 2)) == -2);
		REQUIRE(hhIO.getHelixSeedLength1(2, 2) == 3);
		REQUIRE(hhIO.getHelixSeedLength2(2, 2) == 3);

		// (4,4)
		REQUIRE(hhIO.getHelixSeedE(4, 4) == E_INF);
		REQUIRE(hhIO.getHelixSeedLength1(4, 4) == 0);
		REQUIRE(hhIO.getHelixSeedLength2(4, 4) == 0);


		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		Interaction interaction(r1, r2);
		hhIO.traceBackHelixSeed(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 2);

		// First / last base pair of helixSeed
		REQUIRE(interaction.basePairs.begin()->first == 2);
		REQUIRE(interaction.basePairs.begin()->second == 3);

		REQUIRE(interaction.basePairs.rbegin()->first == 3);
		REQUIRE(interaction.basePairs.rbegin()->second == 2);

		// Case (4,4)
		interaction.clear();
		hhIO.traceBackHelixSeed(interaction, 4, 4);

		REQUIRE(interaction.basePairs.size() == 0);
	}

	SECTION("HelixSeed with Offset: Case 2 - 'A' disrupting complementarity", "[HelixHandlerIdxOffset]") {

		RnaSequence r1("r1", "GGGAGG");
		RnaSequence r2("r2", "CCACCG");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4,  0, 999, 0, false);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ / seedNoGU
		SeedConstraint sC(3, 0, 0, 0, 0, AccessibilityDisabled::ED_UPPER_BOUND, 0, IndexRangeList(""), IndexRangeList(""),
						  "", false, false, false);

		HelixHandler *hhS = new HelixHandlerNoBulgeMax(energy, hC);
		SeedHandler *sH = new SeedHandlerMfe(energy, sC);

		hhS->setSeedHandler(*sH);

		SeedHandlerIdxOffset sHIO(sH);
		HelixHandlerIdxOffset hhIO(hhS);

		// Set the offsets
		sHIO.setOffset1(1);
		sHIO.setOffset2(1);
		hhIO.setOffset1(1);
		hhIO.setOffset2(1);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////   FILLHELIXSEED  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		hhIO.fillHelix(0, energy.size1() - sHIO.getOffset1() - 1, 0, energy.size2() - sHIO.getOffset2() - 1);
		REQUIRE(sHIO.fillSeed(0, energy.size1() - sHIO.getOffset1() - 1, 0, energy.size2() - sHIO.getOffset2() - 1) ==
				0);
		REQUIRE(hhIO.fillHelixSeed(0, energy.size1() - hhIO.getOffset1() - 1, 0,
								   energy.size2() - hhIO.getOffset2() - 1) == 0);

		// (0,0)
		REQUIRE(hhIO.getHelixSeedE(0, 0) == E_INF);
		REQUIRE(hhIO.getHelixSeedLength1(0, 0) == 0);
		REQUIRE(hhIO.getHelixSeedLength2(0, 0) == 0);

		// (0,2)
		REQUIRE(hhIO.getHelixSeedE(0, 2) == E_INF);
		REQUIRE(hhIO.getHelixSeedLength1(0, 2) == 0);
		REQUIRE(hhIO.getHelixSeedLength2(0, 2) == 0);

		// (1,3)
		REQUIRE(hhIO.getHelixSeedE(1, 3) == E_INF);
		REQUIRE(hhIO.getHelixSeedLength1(1, 3) == 0);
		REQUIRE(hhIO.getHelixSeedLength2(1, 3) == 0);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		Interaction interaction(r1, r2);
		hhIO.traceBackHelixSeed(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 0);
	}

	SECTION("HelixSeed with Offset: Case 3 - only seed possible", "[HelixHandlerIdxOffset]") {

		RnaSequence r1("r1", "AAGGGA");
		RnaSequence r2("r2", "ACCCAA");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4,  0, 999, 0, false);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ / seedNoGU
		SeedConstraint sC(3, 0, 0, 0, 0, AccessibilityDisabled::ED_UPPER_BOUND, 0, IndexRangeList(""), IndexRangeList(""),
						  "", false, false, false);

		HelixHandler *hhS = new HelixHandlerNoBulgeMax(energy, hC);
		SeedHandler *sH = new SeedHandlerMfe(energy, sC);

		hhS->setSeedHandler(*sH);

		SeedHandlerIdxOffset sHIO(sH);
		HelixHandlerIdxOffset hhIO(hhS);

		// Set the offsets
		sHIO.setOffset1(1);
		sHIO.setOffset2(1);
		hhIO.setOffset1(1);
		hhIO.setOffset2(1);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////   FILLHELIXSEED  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		hhIO.fillHelix(0, energy.size1() - sHIO.getOffset1() - 1, 0, energy.size2() - sHIO.getOffset2() - 1);
		REQUIRE(sHIO.fillSeed(0, energy.size1() - sHIO.getOffset1() - 1, 0,
							  energy.size2() - sHIO.getOffset2() - 1) == 1);
		REQUIRE(hhIO.fillHelixSeed(0, energy.size1() - hhIO.getOffset1() - 1, 0,
								   energy.size2() - hhIO.getOffset2() - 1) == 1);

		// (0,0)
		REQUIRE(hhIO.getHelixSeedE(0, 0) == E_INF);
		REQUIRE(hhIO.getHelixSeedLength1(0, 0) == 0);
		REQUIRE(hhIO.getHelixSeedLength2(0, 0) == 0);

		// (0,2)
		REQUIRE(hhIO.getHelixSeedE(1, 0) == E_INF);
		REQUIRE(hhIO.getHelixSeedLength1(1, 0) == 0);
		REQUIRE(hhIO.getHelixSeedLength2(1, 0) == 0);

		// (1,3)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(1, 1)) == -2);
		REQUIRE(hhIO.getHelixSeedLength1(1, 1) == 3);
		REQUIRE(hhIO.getHelixSeedLength2(1, 1) == 3);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		Interaction interaction(r1, r2);
		hhIO.traceBackHelixSeed(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 0);

		// Case (1,1)
		interaction.clear();
		hhIO.traceBackHelixSeed(interaction, 1, 1);
		REQUIRE(interaction.basePairs.size() == 1);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 3);
		REQUIRE(interaction.basePairs.begin()->second == 2);

		REQUIRE(interaction.basePairs.rbegin()->first == 3);
		REQUIRE(interaction.basePairs.rbegin()->second == 2);
	}

	SECTION("HelixSeed: Case 4 - unpaired allowed in seed", "[HelixHandlerNoBulgeMax]") {

		RnaSequence r1("r1", "AGGAGG");
		RnaSequence r2("r2", "CCACCA");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4,  0, 999, 0, false);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ / seedNoGU
		SeedConstraint sC(3, 2, 1, 1, 0, AccessibilityDisabled::ED_UPPER_BOUND, 0, IndexRangeList(""), IndexRangeList(""),
						  "", false, false, false);

		HelixHandler *hhS = new HelixHandlerNoBulgeMax(energy, hC);
		SeedHandler *sH = new SeedHandlerMfe(energy, sC);

		hhS->setSeedHandler(*sH);

		SeedHandlerIdxOffset sHIO(sH);
		HelixHandlerIdxOffset hhIO(hhS);

		// Set the offsets
		sHIO.setOffset1(1);
		sHIO.setOffset2(1);
		hhIO.setOffset1(1);
		hhIO.setOffset2(1);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////   FILLHELIXSEED  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		hhIO.fillHelix(0, energy.size1() - sHIO.getOffset1() - 1, 0, energy.size2() - sHIO.getOffset2() - 1);
		REQUIRE(sHIO.fillSeed(0, energy.size1() - sHIO.getOffset1() - 1, 0,
							  energy.size2() - sHIO.getOffset2() - 1) == 4);
		REQUIRE(hhIO.fillHelixSeed(0, energy.size1() - hhIO.getOffset1() - 1, 0,
								   energy.size2() - hhIO.getOffset2() - 1) == 4);

		// (0,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(0, 0)) == -3);
		REQUIRE(hhIO.getHelixSeedLength1(0, 0) == 5);
		REQUIRE(hhIO.getHelixSeedLength2(0, 0) == 5);

		// (0,2)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(1, 0)) == -2);
		REQUIRE(hhIO.getHelixSeedLength1(1, 0) == 4);
		REQUIRE(hhIO.getHelixSeedLength2(1, 0) == 4);

		// (1,3)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(1, 1)) == -2);
		REQUIRE(hhIO.getHelixSeedLength1(1, 1) == 4);
		REQUIRE(hhIO.getHelixSeedLength2(1, 1) == 4);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		Interaction interaction(r1, r2);
		hhIO.traceBackHelixSeed(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 2);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 2);
		REQUIRE(interaction.basePairs.begin()->second == 3);

		REQUIRE(interaction.basePairs.rbegin()->first == 4);
		REQUIRE(interaction.basePairs.rbegin()->second == 1);

		// Case (1,1)
		interaction.clear();
		hhIO.traceBackHelixSeed(interaction, 1, 1);
		REQUIRE(interaction.basePairs.size() == 1);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 4);
		REQUIRE(interaction.basePairs.begin()->second == 1);

		REQUIRE(interaction.basePairs.rbegin()->first == 4);
		REQUIRE(interaction.basePairs.rbegin()->second == 1);

		// Case (0,1)
		interaction.clear();
		hhIO.traceBackHelixSeed(interaction, 0, 1);
		REQUIRE(interaction.basePairs.size() == 1);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 2);
		REQUIRE(interaction.basePairs.begin()->second == 1);

		REQUIRE(interaction.basePairs.rbegin()->first == 2);
		REQUIRE(interaction.basePairs.rbegin()->second == 1);

	}

	SECTION("HelixSeed with Offset: Case 5 - uneven offset", "[HelixHandlerIdxOffset]") {

		RnaSequence r1("r1", "AGGGGG");
		RnaSequence r2("r2", "CCCCCAA");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4,  0, 999, 0, false);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ / seedNoGU
		SeedConstraint sC(3, 0, 0, 0, 0, AccessibilityDisabled::ED_UPPER_BOUND, 0, IndexRangeList(""), IndexRangeList(""),
						  "", false, false, false);

		HelixHandler *hhS = new HelixHandlerNoBulgeMax(energy, hC);
		SeedHandler *sH = new SeedHandlerMfe(energy, sC);

		hhS->setSeedHandler(*sH);

		SeedHandlerIdxOffset sHIO(sH);
		HelixHandlerIdxOffset hhIO(hhS);

		// Set the offsets
		sHIO.setOffset1(1);
		sHIO.setOffset2(2);
		hhIO.setOffset1(1);
		hhIO.setOffset2(2);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////   FILLHELIXSEED  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		hhIO.fillHelix(0, energy.size1() - sHIO.getOffset1() - 1, 0, energy.size2() - sHIO.getOffset2() - 1);
		REQUIRE(sHIO.fillSeed(0, energy.size1() - sHIO.getOffset1() - 1, 0, energy.size2() - sHIO.getOffset2() - 1) ==
				9);
		REQUIRE(hhIO.fillHelixSeed(0, energy.size1() - hhIO.getOffset1() - 1, 0,
								   energy.size2() - hhIO.getOffset2() - 1) == 9);

		// (0,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(0, 0)) == -3);
		REQUIRE(hhIO.getHelixSeedLength1(0, 0) == 4);
		REQUIRE(hhIO.getHelixSeedLength2(0, 0) == 4);

		// (0,1)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(0, 1)) == -3);
		REQUIRE(hhIO.getHelixSeedLength1(0, 1) == 4);
		REQUIRE(hhIO.getHelixSeedLength2(0, 1) == 4);

		// (0,2)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(0, 2)) == -2);
		REQUIRE(hhIO.getHelixSeedLength1(0, 2) == 3);
		REQUIRE(hhIO.getHelixSeedLength2(0, 2) == 3);

		// (0,3)
		REQUIRE(hhIO.getHelixSeedE(0, 3) == E_INF);
		REQUIRE(hhIO.getHelixSeedLength1(0, 3) == 0);
		REQUIRE(hhIO.getHelixSeedLength2(0, 3) == 0);

		// (1,1)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(1, 1)) == -3);
		REQUIRE(hhIO.getHelixSeedLength1(1, 1) == 4);
		REQUIRE(hhIO.getHelixSeedLength2(1, 1) == 4);

		// (2,2)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(2, 2)) == -2);
		REQUIRE(hhIO.getHelixSeedLength1(2, 2) == 3);
		REQUIRE(hhIO.getHelixSeedLength2(2, 2) == 3);

		// (4,4)
		REQUIRE(hhIO.getHelixSeedE(4, 4) == E_INF);
		REQUIRE(hhIO.getHelixSeedLength1(4, 4) == 0);
		REQUIRE(hhIO.getHelixSeedLength2(4, 4) == 0);


		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		Interaction interaction(r1, r2);
		hhIO.traceBackHelixSeed(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 2);

		// First / last base pair of helixSeed
		REQUIRE(interaction.basePairs.begin()->first == 2);
		REQUIRE(interaction.basePairs.begin()->second == 3);

		REQUIRE(interaction.basePairs.rbegin()->first == 3);
		REQUIRE(interaction.basePairs.rbegin()->second == 2);

		// Case (4,4)
		interaction.clear();
		hhIO.traceBackHelixSeed(interaction, 4, 4);

		REQUIRE(interaction.basePairs.size() == 0);
	}
}
