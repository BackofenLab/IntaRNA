
#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/RnaSequence.h"
#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/HelixConstraint.h"
#include "IntaRNA/HelixHandlerUnpaired.h"
#include "IntaRNA/HelixHandlerIdxOffset.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/SeedHandlerMfe.h"

using namespace IntaRNA;

TEST_CASE( "HelixSeed for Unpaired with offset", "[HelixHandlerUnpaired]" ) {

	SECTION("getter", "[HelixHandlerIdxOffset]") {

		RnaSequence r1("r1", "GGGGG");
		RnaSequence r2("r2", "CCCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 10, 2, 999, 0, false);
		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ
		SeedConstraint sC(3, 0, 0, 0, 0, AccessibilityDisabled::ED_UPPER_BOUND, 0, IndexRangeList(""), IndexRangeList(""),
						  "", false, false, false );

		HelixHandler *hhU = new HelixHandlerUnpaired(energy, hC);
		SeedHandler *sH = new SeedHandlerMfe(energy, sC);

		hhU->setSeedHandler(*sH);

		SeedHandlerIdxOffset sHIO (sH);
		HelixHandlerIdxOffset hhIO(hhU);


		// Initial offset of 0
		REQUIRE(hhIO.getOffset1() == 0);
		REQUIRE(hhIO.getOffset2() == 0);

		// Set offset
		hhIO.setOffset1(1);
		hhIO.setOffset2(4);
		REQUIRE(hhIO.getOffset1() == 1);
		REQUIRE(hhIO.getOffset2() == 4);

		// get Constraints
		REQUIRE(hhIO.getConstraint().getMinBasePairs() == 2);
		REQUIRE(hhIO.getConstraint().getMaxBasePairs() == 10);

	}

	SECTION("HelixSeed (Unpaired) with Offset: Case 1 - offset 1", "[HelixHandlerIdxOffset]") {

		RnaSequence r1("r1", "GGGGGG");
		RnaSequence r2("r2", "CCCCCG");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 2, 999, 0, false);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ
		SeedConstraint sC(3, 0, 0, 0, 0, AccessibilityDisabled::ED_UPPER_BOUND, 0, IndexRangeList(""), IndexRangeList(""),
						  "", false, false, false );

		HelixHandler *hhU = new HelixHandlerUnpaired(energy, hC);
		SeedHandler *sH = new SeedHandlerMfe(energy, sC);

		hhU->setSeedHandler(*sH);

		SeedHandlerIdxOffset sHIO(sH);
		HelixHandlerIdxOffset hhIO(hhU);

		// Set the offsets
		sHIO.setOffset1(1);
		sHIO.setOffset2(1);
		hhIO.setOffset1(1);
		hhIO.setOffset2(1);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////   FILLHELIXSEED  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		hhIO.fillHelix(0, energy.size1() - sHIO.getOffset1() - 1, 0, energy.size2() - sHIO.getOffset2() - 1);
		sHIO.fillSeed(0, energy.size1() - sHIO.getOffset1() - 1, 0, energy.size2() - sHIO.getOffset2() - 1);
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

	SECTION("HelixSeed (Unpaired) with Offset: Case 2 - 'A' disrupting complementarity", "[HelixHandlerIdxOffset]") {

		RnaSequence r1("r1", "GGGAGG");
		RnaSequence r2("r2", "CCACCG");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 2, 999, 0, false);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ
		SeedConstraint sC(3, 0, 0, 0, 0, AccessibilityDisabled::ED_UPPER_BOUND, 0, IndexRangeList(""), IndexRangeList(""),
						  "", false, false, false );

		HelixHandler *hhU = new HelixHandlerUnpaired(energy, hC);
		SeedHandler *sH = new SeedHandlerMfe(energy, sC);

		hhU->setSeedHandler(*sH);

		SeedHandlerIdxOffset sHIO(sH);
		HelixHandlerIdxOffset hhIO(hhU);

		// Set the offsets
		sHIO.setOffset1(1);
		sHIO.setOffset2(1);
		hhIO.setOffset1(1);
		hhIO.setOffset2(1);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////   FILLHELIXSEED  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		hhIO.fillHelix(0, energy.size1() - sHIO.getOffset1() - 1, 0, energy.size2() - sHIO.getOffset2() - 1);
		sHIO.fillSeed(0, energy.size1() - sHIO.getOffset1() - 1, 0, energy.size2() - sHIO.getOffset2() - 1);
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

	SECTION("HelixSeed (Unpaired) with Offset: Case 3 - only seed possible", "[HelixHandlerIdxOffset]") {

		RnaSequence r1("r1", "AAGGGA");
		RnaSequence r2("r2", "ACCCAA");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 2, 999, 0, false);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ
		SeedConstraint sC(3, 0, 0, 0, 0, AccessibilityDisabled::ED_UPPER_BOUND, 0, IndexRangeList(""), IndexRangeList(""),
						  "", false, false, false );

		HelixHandler *hhU = new HelixHandlerUnpaired(energy, hC);
		SeedHandler *sH = new SeedHandlerMfe(energy, sC);

		hhU->setSeedHandler(*sH);

		SeedHandlerIdxOffset sHIO(sH);
		HelixHandlerIdxOffset hhIO(hhU);

		// Set the offsets
		sHIO.setOffset1(1);
		sHIO.setOffset2(1);
		hhIO.setOffset1(1);
		hhIO.setOffset2(1);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////   FILLHELIXSEED  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		hhIO.fillHelix(0, energy.size1() - sHIO.getOffset1() - 1, 0, energy.size2() - sHIO.getOffset2() - 1);
		sHIO.fillSeed(0, energy.size1() - sHIO.getOffset1() - 1, 0, energy.size2() - sHIO.getOffset2() - 1);
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

	SECTION("HelixSeed (Unpaired) with Offset : Case 4 - unpaired allowed in seed", "[HelixHandlerUnpaired]") {

		RnaSequence r1("r1", "AGGAGG");
		RnaSequence r2("r2", "CCACCA");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 2, 999, 0, false);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ
		SeedConstraint sC(3, 2, 1, 1, 0, AccessibilityDisabled::ED_UPPER_BOUND, 0, IndexRangeList(""), IndexRangeList(""),
						  "", false, false, false );

		HelixHandler *hhU = new HelixHandlerUnpaired(energy, hC);
		SeedHandler *sH = new SeedHandlerMfe(energy, sC);

		hhU->setSeedHandler(*sH);

		SeedHandlerIdxOffset sHIO(sH);
		HelixHandlerIdxOffset hhIO(hhU);

		// Set the offsets
		sHIO.setOffset1(1);
		sHIO.setOffset2(1);
		hhIO.setOffset1(1);
		hhIO.setOffset2(1);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////   FILLHELIXSEED  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		hhIO.fillHelix(0, energy.size1() - sHIO.getOffset1() - 1, 0, energy.size2() - sHIO.getOffset2() - 1);
		sHIO.fillSeed(0, energy.size1() - sHIO.getOffset1() - 1, 0, energy.size2() - sHIO.getOffset2() - 1);
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

	SECTION("HelixSeed (Unpaired) with Offset: Case 5 - uneven offset", "[HelixHandlerIdxOffset]") {

		RnaSequence r1("r1", "AGGGGG");
		RnaSequence r2("r2", "CCCCCAA");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 2, 999, 0, false);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ
		SeedConstraint sC(3, 0, 0, 0, 0, AccessibilityDisabled::ED_UPPER_BOUND, 0, IndexRangeList(""), IndexRangeList(""),
						  "", false, false, false );

		HelixHandler *hhU = new HelixHandlerUnpaired(energy, hC);
		SeedHandler *sH = new SeedHandlerMfe(energy, sC);

		hhU->setSeedHandler(*sH);

		SeedHandlerIdxOffset sHIO(sH);
		HelixHandlerIdxOffset hhIO(hhU);

		// Set the offsets
		sHIO.setOffset1(1);
		sHIO.setOffset2(2);
		hhIO.setOffset1(1);
		hhIO.setOffset2(2);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////   FILLHELIXSEED  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		hhIO.fillHelix(0, energy.size1() - sHIO.getOffset1() - 1, 0, energy.size2() - sHIO.getOffset2() - 1);
		sHIO.fillSeed(0, energy.size1() - sHIO.getOffset1() - 1, 0, energy.size2() - sHIO.getOffset2() - 1);
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

	SECTION("HelixSeed (Unpaired) with Offset: Case 6 - leading unpaired bases 1", "[HelixHandlerIdxOffset]") {

		RnaSequence r1("r1", "AGGAGGG");
		RnaSequence r2("r2", "CCCACCA");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 5, 2, 999, 0, false);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ
		SeedConstraint sC(3, 0, 0, 0, 0, AccessibilityDisabled::ED_UPPER_BOUND, 0, IndexRangeList(""), IndexRangeList(""),
						  "", false, false, false );

		HelixHandler *hhU = new HelixHandlerUnpaired(energy, hC);
		SeedHandler *sH = new SeedHandlerMfe(energy, sC);

		hhU->setSeedHandler(*sH);

		SeedHandlerIdxOffset sHIO(sH);
		HelixHandlerIdxOffset hhIO(hhU);

		// Set the offsets
		sHIO.setOffset1(1);
		sHIO.setOffset2(1);
		hhIO.setOffset1(1);
		hhIO.setOffset2(1);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////   FILLHELIXSEED  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		hhIO.fillHelix(0, energy.size1() - sHIO.getOffset1() - 1, 0, energy.size2() - sHIO.getOffset2() - 1);
		sHIO.fillSeed(0, energy.size1() - sHIO.getOffset1() - 1, 0, energy.size2() - sHIO.getOffset2() - 1);
		REQUIRE(hhIO.fillHelixSeed(0, energy.size1() - hhIO.getOffset1() - 1, 0,
								   energy.size2() - hhIO.getOffset2() - 1) == 3);

		// (0,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(0, 0)) == -4);
		REQUIRE(hhIO.getHelixSeedLength1(0, 0) == 6);
		REQUIRE(hhIO.getHelixSeedLength2(0, 0) == 6);

		// (1,1)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(1, 1)) == -3);
		REQUIRE(hhIO.getHelixSeedLength1(1, 1) == 5);
		REQUIRE(hhIO.getHelixSeedLength2(1, 1) == 5);

		// (3,3)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(3, 3)) == -2);
		REQUIRE(hhIO.getHelixSeedLength1(3, 3) == 3);
		REQUIRE(hhIO.getHelixSeedLength2(3, 3) == 3);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		Interaction interaction(r1, r2);
		hhIO.traceBackHelixSeed(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 3);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 2);
		REQUIRE(interaction.basePairs.begin()->second == 4);

		REQUIRE(interaction.basePairs.rbegin()->first == 5);
		REQUIRE(interaction.basePairs.rbegin()->second == 1);

		// Case (1,1)
		interaction.clear();
		hhIO.traceBackHelixSeed(interaction, 1, 1);
		REQUIRE(interaction.basePairs.size() == 2);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 4);
		REQUIRE(interaction.basePairs.begin()->second == 2);

		REQUIRE(interaction.basePairs.rbegin()->first == 5);
		REQUIRE(interaction.basePairs.rbegin()->second == 1);

		// Case (0,1)
		interaction.clear();
		hhIO.traceBackHelixSeed(interaction, 0, 1);
		REQUIRE(interaction.basePairs.size() == 0);

	}

	SECTION("HelixSeed (Unpaired) with Offset: Case 7 - leading unpaired bases 2", "[HelixHandlerIdxOffset]") {

		RnaSequence r1("r1", "AGAGGGG");
		RnaSequence r2("r2", "CCCACCA");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 5, 2, 999, 0, false);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ
	    SeedConstraint sC(3, 0, 0, 0, 0, AccessibilityDisabled::ED_UPPER_BOUND, 0, IndexRangeList(""), IndexRangeList(""),
						  "", false, false, false );

		HelixHandler *hhU = new HelixHandlerUnpaired(energy, hC);
		SeedHandler *sH = new SeedHandlerMfe(energy, sC);

		hhU->setSeedHandler(*sH);

		SeedHandlerIdxOffset sHIO(sH);
		HelixHandlerIdxOffset hhIO(hhU);

		// Set the offsets
		sHIO.setOffset1(1);
		sHIO.setOffset2(1);
		hhIO.setOffset1(1);
		hhIO.setOffset2(1);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////   FILLHELIXSEED  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		hhIO.fillHelix(0, energy.size1() - sHIO.getOffset1() - 1, 0, energy.size2() - sHIO.getOffset2() - 1);
		sHIO.fillSeed(0, energy.size1() - sHIO.getOffset1() - 1, 0, energy.size2() - sHIO.getOffset2() - 1);
		REQUIRE(hhIO.fillHelixSeed(0, energy.size1() - hhIO.getOffset1() - 1, 0,
								   energy.size2() - hhIO.getOffset2() - 1) == 5);

		// (0,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(0, 0)) == -4);
		REQUIRE(hhIO.getHelixSeedLength1(0, 0) == 6);
		REQUIRE(hhIO.getHelixSeedLength2(0, 0) == 6);

		// (0,1)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(0, 1)) == -3);
		REQUIRE(hhIO.getHelixSeedLength1(0, 1) == 5);
		REQUIRE(hhIO.getHelixSeedLength2(0, 1) == 5);

//		// (2,0)
//		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(2, 0)) == -3);
//		REQUIRE(hhIO.getHelixSeedLength1(2, 0) == 4);
//		REQUIRE(hhIO.getHelixSeedLength2(2, 0) == 6);

		// (2,1)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(2, 1)) == -3);
		REQUIRE(hhIO.getHelixSeedLength1(2, 1) == 4);
		REQUIRE(hhIO.getHelixSeedLength2(2, 1) == 5);

		// (2,3)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(2, 3)) == -2);
		REQUIRE(hhIO.getHelixSeedLength1(2, 3) == 3);
		REQUIRE(hhIO.getHelixSeedLength2(2, 3) == 3);

		// (3,3)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(3, 3)) == -2);
		REQUIRE(hhIO.getHelixSeedLength1(3, 3) == 3);
		REQUIRE(hhIO.getHelixSeedLength2(3, 3) == 3);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		Interaction interaction(r1, r2);
		hhIO.traceBackHelixSeed(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 3);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 3);
		REQUIRE(interaction.basePairs.begin()->second == 4);

		REQUIRE(interaction.basePairs.rbegin()->first == 5);
		REQUIRE(interaction.basePairs.rbegin()->second == 1);

		// Case (1,1)
		interaction.clear();
		hhIO.traceBackHelixSeed(interaction, 1, 1);
		REQUIRE(interaction.basePairs.size() == 0);

		// Case (0,1)
		interaction.clear();
		hhIO.traceBackHelixSeed(interaction, 0, 1);
		REQUIRE(interaction.basePairs.size() == 2);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 3);
		REQUIRE(interaction.basePairs.begin()->second == 2);

		REQUIRE(interaction.basePairs.rbegin()->first == 4);
		REQUIRE(interaction.basePairs.rbegin()->second == 1);

	}

	SECTION("HelixSeed (Unpaired) with Offset: Case 8 - trailing unpaired bases 1", "[HelixHandlerIdxOffset]") {

		RnaSequence r1("r1", "AGGGAGG");
		RnaSequence r2("r2", "CCACCCA");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 5, 2, 999, 0, false);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ
		SeedConstraint sC(3, 0, 0, 0, 0, AccessibilityDisabled::ED_UPPER_BOUND, 0, IndexRangeList(""), IndexRangeList(""),
						  "", false, false, false );

		HelixHandler *hhU = new HelixHandlerUnpaired(energy, hC);
		SeedHandler *sH = new SeedHandlerMfe(energy, sC);

		hhU->setSeedHandler(*sH);

		SeedHandlerIdxOffset sHIO(sH);
		HelixHandlerIdxOffset hhIO(hhU);

		// Set the offsets
		sHIO.setOffset1(1);
		sHIO.setOffset2(1);
		hhIO.setOffset1(1);
		hhIO.setOffset2(1);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////   FILLHELIXSEED  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		hhIO.fillHelix(0, energy.size1() - sHIO.getOffset1() - 1, 0, energy.size2() - sHIO.getOffset2() - 1);
		sHIO.fillSeed(0, energy.size1() - sHIO.getOffset1() - 1, 0, energy.size2() - sHIO.getOffset2() - 1);
		REQUIRE(hhIO.fillHelixSeed(0, energy.size1() - hhIO.getOffset1() - 1,
								   0, energy.size2() - hhIO.getOffset2() - 1) == 1);

		// (0,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(0, 0)) == -4);
		REQUIRE(hhIO.getHelixSeedLength1(0, 0) == 6);
		REQUIRE(hhIO.getHelixSeedLength2(0, 0) == 6);


		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		Interaction interaction(r1, r2);
		hhIO.traceBackHelixSeed(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 3);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 2);
		REQUIRE(interaction.basePairs.begin()->second == 4);

		REQUIRE(interaction.basePairs.rbegin()->first == 5);
		REQUIRE(interaction.basePairs.rbegin()->second == 1);

		// Case (1,1)
		interaction.clear();
		hhIO.traceBackHelixSeed(interaction, 1, 1);
		REQUIRE(interaction.basePairs.size() == 0);
	}

	SECTION("HelixSeed (Unpaired) with Offset: Case 9 - trailing unpaired bases 2", "[HelixHandlerIdxOffset]") {

		RnaSequence r1("r1", "AGGGAGG");
		RnaSequence r2("r2", "CCACCCCA");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 5, 2, 999, 0, false);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ
		SeedConstraint sC(3, 0, 0, 0, 0, AccessibilityDisabled::ED_UPPER_BOUND, 0, IndexRangeList(""), IndexRangeList(""),
						  "", false, false, false );

		HelixHandler *hhU = new HelixHandlerUnpaired(energy, hC);
		SeedHandler *sH = new SeedHandlerMfe(energy, sC);

		hhU->setSeedHandler(*sH);

		SeedHandlerIdxOffset sHIO(sH);
		HelixHandlerIdxOffset hhIO(hhU);

		// Set the offsets
		sHIO.setOffset1(1);
		sHIO.setOffset2(1);
		hhIO.setOffset1(1);
		hhIO.setOffset2(1);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////   FILLHELIXSEED  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		hhIO.fillHelix(0, energy.size1() - sHIO.getOffset1() - 1, 0, energy.size2() - sHIO.getOffset2() - 1);
		sHIO.fillSeed(0, energy.size1() - sHIO.getOffset1() - 1, 0, energy.size2() - sHIO.getOffset2() - 1);
		REQUIRE(hhIO.fillHelixSeed(0, energy.size1() - hhIO.getOffset1() - 1, 0,
								   energy.size2() - hhIO.getOffset2() - 1) == 2);
		// (0,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(0, 0)) == -4);
		REQUIRE(hhIO.getHelixSeedLength1(0, 0) == 6);
		REQUIRE(hhIO.getHelixSeedLength2(0, 0) == 6);

		// (0,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(0, 1)) == -4);
		REQUIRE(hhIO.getHelixSeedLength1(0, 1) == 6);
		REQUIRE(hhIO.getHelixSeedLength2(0, 1) == 6);


		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		Interaction interaction(r1, r2);
		hhIO.traceBackHelixSeed(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 3);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 2);
		REQUIRE(interaction.basePairs.begin()->second == 5);

		REQUIRE(interaction.basePairs.rbegin()->first == 5);
		REQUIRE(interaction.basePairs.rbegin()->second == 3);

		// Case (0,1)
		interaction.clear();
		hhIO.traceBackHelixSeed(interaction, 0, 1);
		REQUIRE(interaction.basePairs.size() == 3);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 2);
		REQUIRE(interaction.basePairs.begin()->second == 4);

		REQUIRE(interaction.basePairs.rbegin()->first == 5);
		REQUIRE(interaction.basePairs.rbegin()->second == 1);

		// Case (1,1)
		interaction.clear();
		hhIO.traceBackHelixSeed(interaction, 1, 1);
		REQUIRE(interaction.basePairs.size() == 0);
	}

	SECTION("HelixSeed (Unpaired) with Offset: Case 10 - Leading + trailing unpaired bases", "[HelixHandlerIdxOffset]") {

		RnaSequence r1("r1", "AGAGGGG");
		RnaSequence r2("r2", "CACCCCA");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 5, 2, 999, 0, false);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ
		SeedConstraint sC(3, 0, 0, 0, 0, AccessibilityDisabled::ED_UPPER_BOUND, 0, IndexRangeList(""), IndexRangeList(""),
						  "", false, false, false );

		HelixHandler *hhU = new HelixHandlerUnpaired(energy, hC);
		SeedHandler *sH = new SeedHandlerMfe(energy, sC);

		hhU->setSeedHandler(*sH);

		SeedHandlerIdxOffset sHIO(sH);
		HelixHandlerIdxOffset hhIO(hhU);

		// Set the offsets
		sHIO.setOffset1(1);
		sHIO.setOffset2(1);
		hhIO.setOffset1(1);
		hhIO.setOffset2(1);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////   FILLHELIXSEED  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		hhIO.fillHelix(0, energy.size1() - sHIO.getOffset1() - 1, 0, energy.size2() - sHIO.getOffset2() - 1);
		sHIO.fillSeed(0, energy.size1() - sHIO.getOffset1() - 1, 0, energy.size2() - sHIO.getOffset2() - 1);
		REQUIRE(hhIO.fillHelixSeed(0, energy.size1() - hhIO.getOffset1() - 1, 0,
								   energy.size2() - hhIO.getOffset2() - 1) == 5);

		// (0,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(0, 0)) == -4);
		REQUIRE(hhIO.getHelixSeedLength1(0, 0) == 6);
		REQUIRE(hhIO.getHelixSeedLength2(0, 0) == 6);

		// (0,1) - Not working
		REQUIRE(hhIO.getHelixSeedE(0, 1) == E_INF);
		REQUIRE(hhIO.getHelixSeedLength1(0, 1) == 0);
		REQUIRE(hhIO.getHelixSeedLength2(0, 1) == 0);

		// (2,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(2, 0)) == -3);
		REQUIRE(hhIO.getHelixSeedLength1(2, 0) == 4);
		REQUIRE(hhIO.getHelixSeedLength2(2, 0) == 4);

		// (2,1)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(2, 1)) == -3);
		REQUIRE(hhIO.getHelixSeedLength1(2, 1) == 4);
		REQUIRE(hhIO.getHelixSeedLength2(2, 1) == 5);

		// (3,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(3, 0)) == -2);
		REQUIRE(hhIO.getHelixSeedLength1(3, 0) == 3);
		REQUIRE(hhIO.getHelixSeedLength2(3, 0) == 3);

		// (3,1)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(3, 1)) == -2);
		REQUIRE(hhIO.getHelixSeedLength1(3, 1) == 3);
		REQUIRE(hhIO.getHelixSeedLength2(3, 1) == 3);


		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		Interaction interaction(r1, r2);
		hhIO.traceBackHelixSeed(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 3);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 3);
		REQUIRE(interaction.basePairs.begin()->second == 4);

		REQUIRE(interaction.basePairs.rbegin()->first == 5);
		REQUIRE(interaction.basePairs.rbegin()->second == 2);

		// Case (3,1)
		interaction.clear();
		hhIO.traceBackHelixSeed(interaction, 3, 1);
		REQUIRE(interaction.basePairs.size() == 1);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 5);
		REQUIRE(interaction.basePairs.begin()->second == 3);

		REQUIRE(interaction.basePairs.rbegin()->first == 5);
		REQUIRE(interaction.basePairs.rbegin()->second == 3);

		// Case (1,1)
		interaction.clear();
		hhIO.traceBackHelixSeed(interaction, 0, 1);
		REQUIRE(interaction.basePairs.size() == 0);
	}

	SECTION("HelixSeed (Unpaired) with Offset: Case 11 - Leading + trailing unpaired bases + seed allows 1 unpaired", "[HelixHandlerIdxOffset]") {

		RnaSequence r1("r1", "AGAGGAGG");
		RnaSequence r2("r2", "CACCCCA");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 5, 2, 999, 0, false);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ
		SeedConstraint sC(3, 1, 1, 0, 0, AccessibilityDisabled::ED_UPPER_BOUND, 0, IndexRangeList(""), IndexRangeList(""),
						  "", false, false, false );

		HelixHandler *hhU = new HelixHandlerUnpaired(energy, hC);
		SeedHandler *sH = new SeedHandlerMfe(energy, sC);

		hhU->setSeedHandler(*sH);

		SeedHandlerIdxOffset sHIO(sH);
		HelixHandlerIdxOffset hhIO(hhU);

		// Set the offsets
		sHIO.setOffset1(1);
		sHIO.setOffset2(1);
		hhIO.setOffset1(1);
		hhIO.setOffset2(1);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////   FILLHELIXSEED  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		hhIO.fillHelix(0, energy.size1() - sHIO.getOffset1() - 1, 0, energy.size2() - sHIO.getOffset2() - 1);
		sHIO.fillSeed(0, energy.size1() - sHIO.getOffset1() - 1, 0, energy.size2() - sHIO.getOffset2() - 1);
		REQUIRE(hhIO.fillHelixSeed(0, energy.size1() - hhIO.getOffset1() - 1, 0,
								   energy.size2() - hhIO.getOffset2() - 1) == 6);

		// (0,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(0, 0)) == -4);
		REQUIRE(hhIO.getHelixSeedLength1(0, 0) == 7);
		REQUIRE(hhIO.getHelixSeedLength2(0, 0) == 6);

		// (0,1)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(0, 1)) == -3);
		REQUIRE(hhIO.getHelixSeedLength1(0, 1) == 6);
		REQUIRE(hhIO.getHelixSeedLength2(0, 1) == 5);

		// (2,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(2, 0)) == -3);
		REQUIRE(hhIO.getHelixSeedLength1(2, 0) == 5);
		REQUIRE(hhIO.getHelixSeedLength2(2, 0) == 4);

		// (2,1)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(2, 1)) == -3);
		REQUIRE(hhIO.getHelixSeedLength1(2, 1) == 5);
		REQUIRE(hhIO.getHelixSeedLength2(2, 1) == 5);

		// (3,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(3, 0)) == -2);
		REQUIRE(hhIO.getHelixSeedLength1(3, 0) == 4);
		REQUIRE(hhIO.getHelixSeedLength2(3, 0) == 3);

		// (3,1)
		REQUIRE(E_2_Ekcal(hhIO.getHelixSeedE(3, 1)) == -2);
		REQUIRE(hhIO.getHelixSeedLength1(3, 1) == 4);
		REQUIRE(hhIO.getHelixSeedLength2(3, 1) == 3);


		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		Interaction interaction(r1, r2);
		hhIO.traceBackHelixSeed(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 3);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 3);
		REQUIRE(interaction.basePairs.begin()->second == 4);

		REQUIRE(interaction.basePairs.rbegin()->first == 6);
		REQUIRE(interaction.basePairs.rbegin()->second == 2);

		// Case (3,1)
		interaction.clear();
		hhIO.traceBackHelixSeed(interaction, 3, 1);
		REQUIRE(interaction.basePairs.size() == 1);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 6);
		REQUIRE(interaction.basePairs.begin()->second == 3);

		REQUIRE(interaction.basePairs.rbegin()->first == 6);
		REQUIRE(interaction.basePairs.rbegin()->second == 3);

		// Case (1,1)
		interaction.clear();
		hhIO.traceBackHelixSeed(interaction, 0, 1);
		REQUIRE(interaction.basePairs.size() == 2);
		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 3);
		REQUIRE(interaction.basePairs.begin()->second == 3);

		REQUIRE(interaction.basePairs.rbegin()->first == 4);
		REQUIRE(interaction.basePairs.rbegin()->second == 2);
	}
}
