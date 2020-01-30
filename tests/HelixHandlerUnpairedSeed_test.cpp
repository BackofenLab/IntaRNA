
#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/RnaSequence.h"
#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/HelixConstraint.h"
#include "IntaRNA/HelixHandlerUnpaired.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/SeedHandlerMfe.h"

using namespace IntaRNA;

TEST_CASE( "HelixHandlerUnpairedSeed", "[HelixHandlerUnpaired]" ) {

	SECTION("HelixSeed: Case 1 - all complementary", "[HelixHandlerUnpaired]") {

		RnaSequence r1("r1", "GGGGG");
		RnaSequence r2("r2", "CCCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 2, 999, 0, false);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ
		SeedConstraint sC(3, 0, 0, 0, 0, AccessibilityDisabled::ED_UPPER_BOUND, 0, IndexRangeList(""), IndexRangeList(""),
						  "", false, false, false );

		SeedHandlerMfe sH(energy, sC);
		HelixHandlerUnpaired hhU(energy, hC);

		sH.fillSeed(0, energy.size1() - 1, 0, energy.size2() - 1);
		hhU.setSeedHandler(sH);

		hhU.fillHelix(0, energy.size1() - 1, 0, energy.size2() - 1);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////   FILLHELIXSEED  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		REQUIRE(hhU.fillHelixSeed(0, energy.size1() - 1, 0, energy.size2() - 1) == 9);

		// (0,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(0, 0)) == -3);
		REQUIRE(hhU.getHelixSeedLength1(0, 0) == 4);
		REQUIRE(hhU.getHelixSeedLength2(0, 0) == 4);

		// (0,1)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(0, 1)) == -3);
		REQUIRE(hhU.getHelixSeedLength1(0, 1) == 4);
		REQUIRE(hhU.getHelixSeedLength2(0, 1) == 4);

		// (0,2)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(0, 2)) == -2);
		REQUIRE(hhU.getHelixSeedLength1(0, 2) == 3);
		REQUIRE(hhU.getHelixSeedLength2(0, 2) == 3);

		// (0,3)
		REQUIRE(hhU.getHelixSeedE(0, 3) == E_INF);
		REQUIRE(hhU.getHelixSeedLength1(0, 3) == 0);
		REQUIRE(hhU.getHelixSeedLength2(0, 3) == 0);

		// (1,1)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(1, 1)) == -3);
		REQUIRE(hhU.getHelixSeedLength1(1, 1) == 4);
		REQUIRE(hhU.getHelixSeedLength2(1, 1) == 4);

		// (2,2)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(2, 2)) == -2);
		REQUIRE(hhU.getHelixSeedLength1(2, 2) == 3);
		REQUIRE(hhU.getHelixSeedLength2(2, 2) == 3);

		// (4,4)
		REQUIRE(hhU.getHelixSeedE(4, 4) == E_INF);
		REQUIRE(hhU.getHelixSeedLength1(4, 4) == 0);
		REQUIRE(hhU.getHelixSeedLength2(4, 4) == 0);


		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		Interaction interaction(r1, r2);
		hhU.traceBackHelixSeed(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 2);

		// First / last base pair of helixSeed
		REQUIRE(interaction.basePairs.begin()->first == 1);
		REQUIRE(interaction.basePairs.begin()->second == 3);

		REQUIRE(interaction.basePairs.rbegin()->first == 2);
		REQUIRE(interaction.basePairs.rbegin()->second == 2);

		// Case (4,4)
		interaction.clear();
		hhU.traceBackHelixSeed(interaction, 4, 4);

		REQUIRE(interaction.basePairs.size() == 0);
	}

	SECTION("HelixSeed: Case 2 - 'A' disrupting complementarity", "[HelixHandlerUnpaired]") {

		RnaSequence r1("r1", "GGAGG");
		RnaSequence r2("r2", "CCACC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 2, 999, 0, false);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ
		SeedConstraint sC(3, 0, 0, 0, 0, AccessibilityDisabled::ED_UPPER_BOUND, 0, IndexRangeList(""), IndexRangeList(""),
						  "", false, false, false );

		SeedHandlerMfe sH(energy, sC);
		HelixHandlerUnpaired hhU(energy, hC);

		sH.fillSeed(0, energy.size1() - 1, 0, energy.size2() - 1);
		hhU.setSeedHandler(sH);
		hhU.fillHelix(0, energy.size1() - 1, 0, energy.size2() - 1);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////   FILLHELIXSEED  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		REQUIRE(hhU.fillHelixSeed(0, energy.size1() - 1, 0, energy.size2() - 1) == 0);

		// (0,0)
		REQUIRE(hhU.getHelixSeedE(0, 0) == E_INF);
		REQUIRE(hhU.getHelixSeedLength1(0, 0) == 0);
		REQUIRE(hhU.getHelixSeedLength2(0, 0) == 0);

		// (0,2)
		REQUIRE(hhU.getHelixSeedE(0, 2) == E_INF);
		REQUIRE(hhU.getHelixSeedLength1(0, 2) == 0);
		REQUIRE(hhU.getHelixSeedLength2(0, 2) == 0);

		// (1,3)
		REQUIRE(hhU.getHelixSeedE(1, 3) == E_INF);
		REQUIRE(hhU.getHelixSeedLength1(1, 3) == 0);
		REQUIRE(hhU.getHelixSeedLength2(1, 3) == 0);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		Interaction interaction(r1, r2);
		hhU.traceBackHelixSeed(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 0);
	}

	SECTION("HelixSeed: Case 3 - only seed possible", "[HelixHandlerUnpaired]") {

		RnaSequence r1("r1", "AGGGA");
		RnaSequence r2("r2", "ACCCA");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 2, 999, 0, false);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ
		SeedConstraint sC(3, 0, 0, 0, 0, AccessibilityDisabled::ED_UPPER_BOUND, 0, IndexRangeList(""), IndexRangeList(""),
						  "", false, false, false );

		SeedHandlerMfe sH(energy, sC);
		HelixHandlerUnpaired hhU(energy, hC);

		sH.fillSeed(0, energy.size1() - 1, 0, energy.size2() - 1);
		hhU.setSeedHandler(sH);

		hhU.fillHelix(0, energy.size1() - 1, 0, energy.size2() - 1);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////   FILLHELIXSEED  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		REQUIRE(hhU.fillHelixSeed(0, energy.size1() - 1, 0, energy.size2() - 1) == 1);

		// (0,0)
		REQUIRE(hhU.getHelixSeedE(0, 0) == E_INF);
		REQUIRE(hhU.getHelixSeedLength1(0, 0) == 0);
		REQUIRE(hhU.getHelixSeedLength2(0, 0) == 0);

		// (0,2)
		REQUIRE(hhU.getHelixSeedE(1, 0) == E_INF);
		REQUIRE(hhU.getHelixSeedLength1(1, 0) == 0);
		REQUIRE(hhU.getHelixSeedLength2(1, 0) == 0);

		// (1,3)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(1, 1)) == -2);
		REQUIRE(hhU.getHelixSeedLength1(1, 1) == 3);
		REQUIRE(hhU.getHelixSeedLength2(1, 1) == 3);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		Interaction interaction(r1, r2);
		hhU.traceBackHelixSeed(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 0);

		// Case (1,1)
		interaction.clear();
		hhU.traceBackHelixSeed(interaction, 1, 1);
		REQUIRE(interaction.basePairs.size() == 1);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 2);
		REQUIRE(interaction.basePairs.begin()->second == 2);

		REQUIRE(interaction.basePairs.rbegin()->first == 2);
		REQUIRE(interaction.basePairs.rbegin()->second == 2);
	}

	SECTION("HelixSeed: Case4 - unpaired allowed in seed", "[HelixHandlerUnpaired]") {

		RnaSequence r1("r1", "GGAGG");
		RnaSequence r2("r2", "CCACC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 2, 999, 0, false);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ
		SeedConstraint sC(3, 2, 1, 1, 0, AccessibilityDisabled::ED_UPPER_BOUND, 0, IndexRangeList(""), IndexRangeList(""),
						  "", false, false, false );

		SeedHandlerMfe sH(energy, sC);
		HelixHandlerUnpaired hhU(energy, hC);

		sH.fillSeed(0, energy.size1() - 1, 0, energy.size2() - 1);
		hhU.setSeedHandler(sH);
		hhU.fillHelix(0, energy.size1() - 1, 0, energy.size2() - 1);


		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////   FILLHELIXSEED  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		REQUIRE(hhU.fillHelixSeed(0, energy.size1() - 1, 0, energy.size2() - 1) == 4);

		// (0,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(0, 0)) == -3);
		REQUIRE(hhU.getHelixSeedLength1(0, 0) == 5);
		REQUIRE(hhU.getHelixSeedLength2(0, 0) == 5);

		// (0,2)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(1, 0)) == -2);
		REQUIRE(hhU.getHelixSeedLength1(1, 0) == 4);
		REQUIRE(hhU.getHelixSeedLength2(1, 0) == 4);

		// (1,3)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(1, 1)) == -2);
		REQUIRE(hhU.getHelixSeedLength1(1, 1) == 4);
		REQUIRE(hhU.getHelixSeedLength2(1, 1) == 4);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		Interaction interaction(r1, r2);
		hhU.traceBackHelixSeed(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 2);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 1);
		REQUIRE(interaction.basePairs.begin()->second == 3);

		REQUIRE(interaction.basePairs.rbegin()->first == 3);
		REQUIRE(interaction.basePairs.rbegin()->second == 1);

		// Case (1,1)
		interaction.clear();
		hhU.traceBackHelixSeed(interaction, 1, 1);
		REQUIRE(interaction.basePairs.size() == 1);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 3);
		REQUIRE(interaction.basePairs.begin()->second == 1);

		REQUIRE(interaction.basePairs.rbegin()->first == 3);
		REQUIRE(interaction.basePairs.rbegin()->second == 1);

		// Case (0,1)
		interaction.clear();
		hhU.traceBackHelixSeed(interaction, 0, 1);
		REQUIRE(interaction.basePairs.size() == 1);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 1);
		REQUIRE(interaction.basePairs.begin()->second == 1);

		REQUIRE(interaction.basePairs.rbegin()->first == 1);
		REQUIRE(interaction.basePairs.rbegin()->second == 1);

	}

	SECTION("HelixSeed: Case 5 - leading unpaired bases 1", "[HelixHandlerUnpaired]") {

		RnaSequence r1("r1", "GGAGGG");
		RnaSequence r2("r2", "CCCACC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 5, 2, 999, 0, false);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ
		SeedConstraint sC(3, 0, 0, 0, 0, AccessibilityDisabled::ED_UPPER_BOUND, 0, IndexRangeList(""), IndexRangeList(""),
						  "", false, false, false );

		SeedHandlerMfe sH(energy, sC);
		HelixHandlerUnpaired hhU(energy, hC);

		sH.fillSeed(0, energy.size1() - 1, 0, energy.size2() - 1);
		hhU.setSeedHandler(sH);
		hhU.fillHelix(0, energy.size1() - 1, 0, energy.size2() - 1);


		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////   FILLHELIXSEED  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		REQUIRE(hhU.fillHelixSeed(0, energy.size1() - 1, 0, energy.size2() - 1) == 3);

		// (0,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(0, 0)) == -4);
		REQUIRE(hhU.getHelixSeedLength1(0, 0) == 6);
		REQUIRE(hhU.getHelixSeedLength2(0, 0) == 6);

		// (1,1)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(1, 1)) == -3);
		REQUIRE(hhU.getHelixSeedLength1(1, 1) == 5);
		REQUIRE(hhU.getHelixSeedLength2(1, 1) == 5);

		// (3,3)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(3, 3)) == -2);
		REQUIRE(hhU.getHelixSeedLength1(3, 3) == 3);
		REQUIRE(hhU.getHelixSeedLength2(3, 3) == 3);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		Interaction interaction(r1, r2);
		hhU.traceBackHelixSeed(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 3);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 1);
		REQUIRE(interaction.basePairs.begin()->second == 4);

		REQUIRE(interaction.basePairs.rbegin()->first == 4);
		REQUIRE(interaction.basePairs.rbegin()->second == 1);

		// Case (1,1)
		interaction.clear();
		hhU.traceBackHelixSeed(interaction, 1, 1);
		REQUIRE(interaction.basePairs.size() == 2);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 3);
		REQUIRE(interaction.basePairs.begin()->second == 2);

		REQUIRE(interaction.basePairs.rbegin()->first == 4);
		REQUIRE(interaction.basePairs.rbegin()->second == 1);

		// Case (0,1)
		interaction.clear();
		hhU.traceBackHelixSeed(interaction, 0, 1);
		REQUIRE(interaction.basePairs.size() == 0);

	}

	SECTION("HelixSeed: Case 6 - leading unpaired bases 2", "[HelixHandlerUnpaired]") {

		RnaSequence r1("r1", "GAGGGG");
		RnaSequence r2("r2", "CCCACC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 5, 2, 999, 0, false);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ
		SeedConstraint sC(3, 0, 0, 0, 0, AccessibilityDisabled::ED_UPPER_BOUND, 0, IndexRangeList(""), IndexRangeList(""),
						  "", false, false, false );

		SeedHandlerMfe sH(energy, sC);
		HelixHandlerUnpaired hhU(energy, hC);

		sH.fillSeed(0, energy.size1() - 1, 0, energy.size2() - 1);
		hhU.setSeedHandler(sH);
		hhU.fillHelix(0, energy.size1() - 1, 0, energy.size2() - 1);


		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////   FILLHELIXSEED  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		REQUIRE(hhU.fillHelixSeed(0, energy.size1() - 1, 0, energy.size2() - 1) == 5);

		// (0,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(0, 0)) == -4);
		REQUIRE(hhU.getHelixSeedLength1(0, 0) == 6);
		REQUIRE(hhU.getHelixSeedLength2(0, 0) == 6);

		// (0,1)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(0, 1)) == -3);
		REQUIRE(hhU.getHelixSeedLength1(0, 1) == 5);
		REQUIRE(hhU.getHelixSeedLength2(0, 1) == 5);

//		// (2,0)
//		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(2, 0)) == -3);
//		REQUIRE(hhU.getHelixSeedLength1(2, 0) == 4);
//		REQUIRE(hhU.getHelixSeedLength2(2, 0) == 6);

		// (2,1)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(2, 1)) == -3);
		REQUIRE(hhU.getHelixSeedLength1(2, 1) == 4);
		REQUIRE(hhU.getHelixSeedLength2(2, 1) == 5);

		// (2,3)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(2, 3)) == -2);
		REQUIRE(hhU.getHelixSeedLength1(2, 3) == 3);
		REQUIRE(hhU.getHelixSeedLength2(2, 3) == 3);

		// (3,3)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(3, 3)) == -2);
		REQUIRE(hhU.getHelixSeedLength1(3, 3) == 3);
		REQUIRE(hhU.getHelixSeedLength2(3, 3) == 3);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		Interaction interaction(r1, r2);
		hhU.traceBackHelixSeed(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 3);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 2);
		REQUIRE(interaction.basePairs.begin()->second == 4);

		REQUIRE(interaction.basePairs.rbegin()->first == 4);
		REQUIRE(interaction.basePairs.rbegin()->second == 1);

		// Case (1,1)
		interaction.clear();
		hhU.traceBackHelixSeed(interaction, 1, 1);
		REQUIRE(interaction.basePairs.size() == 0);

		// Case (0,1)
		interaction.clear();
		hhU.traceBackHelixSeed(interaction, 0, 1);
		REQUIRE(interaction.basePairs.size() == 2);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 2);
		REQUIRE(interaction.basePairs.begin()->second == 2);

		REQUIRE(interaction.basePairs.rbegin()->first == 3);
		REQUIRE(interaction.basePairs.rbegin()->second == 1);

	}

	SECTION("HelixSeed: Case 7 - trailing unpaired bases 1", "[HelixHandlerUnpaired]") {

		RnaSequence r1("r1", "GGGAGG");
		RnaSequence r2("r2", "CCACCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 5, 2, 999, 0, false);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ
		SeedConstraint sC(3, 0, 0, 0, 0, AccessibilityDisabled::ED_UPPER_BOUND, 0, IndexRangeList(""), IndexRangeList(""),
						  "", false, false, false );

		SeedHandlerMfe sH(energy, sC);
		HelixHandlerUnpaired hhU(energy, hC);

		sH.fillSeed(0, energy.size1() - 1, 0, energy.size2() - 1);
		hhU.setSeedHandler(sH);
		hhU.fillHelix(0, energy.size1() - 1, 0, energy.size2() - 1);


		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////   FILLHELIXSEED  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		REQUIRE(hhU.fillHelixSeed(0, energy.size1() - 1, 0, energy.size2() - 1) == 1);

		// (0,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(0, 0)) == -4);
		REQUIRE(hhU.getHelixSeedLength1(0, 0) == 6);
		REQUIRE(hhU.getHelixSeedLength2(0, 0) == 6);


		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		Interaction interaction(r1, r2);
		hhU.traceBackHelixSeed(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 3);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 1);
		REQUIRE(interaction.basePairs.begin()->second == 4);

		REQUIRE(interaction.basePairs.rbegin()->first == 4);
		REQUIRE(interaction.basePairs.rbegin()->second == 1);

		// Case (1,1)
		interaction.clear();
		hhU.traceBackHelixSeed(interaction, 1, 1);
		REQUIRE(interaction.basePairs.size() == 0);
	}

	SECTION("HelixSeed: Case 7 - trailing unpaired bases 2", "[HelixHandlerUnpaired]") {

		RnaSequence r1("r1", "GGGAGG");
		RnaSequence r2("r2", "CCACCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 5, 2, 999, 0, false);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ
		SeedConstraint sC(3, 0, 0, 0, 0, AccessibilityDisabled::ED_UPPER_BOUND, 0, IndexRangeList(""), IndexRangeList(""),
						  "", false, false, false );

		SeedHandlerMfe sH(energy, sC);
		HelixHandlerUnpaired hhU(energy, hC);

		sH.fillSeed(0, energy.size1() - 1, 0, energy.size2() - 1);
		hhU.setSeedHandler(sH);
		hhU.fillHelix(0, energy.size1() - 1, 0, energy.size2() - 1);


		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////   FILLHELIXSEED  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		REQUIRE(hhU.fillHelixSeed(0, energy.size1() - 1, 0, energy.size2() - 1) == 2);

		// (0,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(0, 0)) == -4);
		REQUIRE(hhU.getHelixSeedLength1(0, 0) == 6);
		REQUIRE(hhU.getHelixSeedLength2(0, 0) == 6);

		// (0,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(0, 1)) == -4);
		REQUIRE(hhU.getHelixSeedLength1(0, 1) == 6);
		REQUIRE(hhU.getHelixSeedLength2(0, 1) == 6);


		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		Interaction interaction(r1, r2);
		hhU.traceBackHelixSeed(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 3);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 1);
		REQUIRE(interaction.basePairs.begin()->second == 5);

		REQUIRE(interaction.basePairs.rbegin()->first == 4);
		REQUIRE(interaction.basePairs.rbegin()->second == 3);

		// Case (0,1)
		interaction.clear();
		hhU.traceBackHelixSeed(interaction, 0, 1);
		REQUIRE(interaction.basePairs.size() == 3);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 1);
		REQUIRE(interaction.basePairs.begin()->second == 4);

		REQUIRE(interaction.basePairs.rbegin()->first == 4);
		REQUIRE(interaction.basePairs.rbegin()->second == 1);

		// Case (1,1)
		interaction.clear();
		hhU.traceBackHelixSeed(interaction, 1, 1);
		REQUIRE(interaction.basePairs.size() == 0);
	}

	SECTION("HelixSeed: Case 8 - Leading + trailing unpaired bases", "[HelixHandlerUnpaired]") {

		RnaSequence r1("r1", "GAGGGG");
		RnaSequence r2("r2", "CACCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 5, 2, 999, 0, false);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ
		SeedConstraint sC(3, 0, 0, 0, 0, AccessibilityDisabled::ED_UPPER_BOUND, 0, IndexRangeList(""), IndexRangeList(""),
						  "", false, false, false );

		SeedHandlerMfe sH(energy, sC);
		HelixHandlerUnpaired hhU(energy, hC);

		sH.fillSeed(0, energy.size1() - 1, 0, energy.size2() - 1);
		hhU.setSeedHandler(sH);
		hhU.fillHelix(0, energy.size1() - 1, 0, energy.size2() - 1);


		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////   FILLHELIXSEED  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		REQUIRE(hhU.fillHelixSeed(0, energy.size1() - 1, 0, energy.size2() - 1) == 5);

		// (0,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(0, 0)) == -4);
		REQUIRE(hhU.getHelixSeedLength1(0, 0) == 6);
		REQUIRE(hhU.getHelixSeedLength2(0, 0) == 6);

		// (0,1) - Not working
		REQUIRE(hhU.getHelixSeedE(0, 1) == E_INF);
		REQUIRE(hhU.getHelixSeedLength1(0, 1) == 0);
		REQUIRE(hhU.getHelixSeedLength2(0, 1) == 0);

		// (2,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(2, 0)) == -3);
		REQUIRE(hhU.getHelixSeedLength1(2, 0) == 4);
		REQUIRE(hhU.getHelixSeedLength2(2, 0) == 4);

		// (2,1)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(2, 1)) == -3);
		REQUIRE(hhU.getHelixSeedLength1(2, 1) == 4);
		REQUIRE(hhU.getHelixSeedLength2(2, 1) == 5);

		// (3,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(3, 0)) == -2);
		REQUIRE(hhU.getHelixSeedLength1(3, 0) == 3);
		REQUIRE(hhU.getHelixSeedLength2(3, 0) == 3);

		// (3,1)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(3, 1)) == -2);
		REQUIRE(hhU.getHelixSeedLength1(3, 1) == 3);
		REQUIRE(hhU.getHelixSeedLength2(3, 1) == 3);


		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		Interaction interaction(r1, r2);
		hhU.traceBackHelixSeed(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 3);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 2);
		REQUIRE(interaction.basePairs.begin()->second == 4);

		REQUIRE(interaction.basePairs.rbegin()->first == 4);
		REQUIRE(interaction.basePairs.rbegin()->second == 2);

		// Case (3,1)
		interaction.clear();
		hhU.traceBackHelixSeed(interaction, 3, 1);
		REQUIRE(interaction.basePairs.size() == 1);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 4);
		REQUIRE(interaction.basePairs.begin()->second == 3);

		REQUIRE(interaction.basePairs.rbegin()->first == 4);
		REQUIRE(interaction.basePairs.rbegin()->second == 3);

		// Case (1,1)
		interaction.clear();
		hhU.traceBackHelixSeed(interaction, 0, 1);
		REQUIRE(interaction.basePairs.size() == 0);
	}

	SECTION("HelixSeed: Case 9 - Leading + trailing unpaired bases + seed allows 1 unpaired",
			"[HelixHandlerUnpaired]") {

		RnaSequence r1("r1", "GAGGAGG");
		RnaSequence r2("r2", "CACCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 5, 2, 999, 0, false);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ
		SeedConstraint sC(3, 1, 1, 0, 0, AccessibilityDisabled::ED_UPPER_BOUND, 0, IndexRangeList(""), IndexRangeList(""),
						  "", false, false, false );

		SeedHandlerMfe sH(energy, sC);
		HelixHandlerUnpaired hhU(energy, hC);

		sH.fillSeed(0, energy.size1() - 1, 0, energy.size2() - 1);
		hhU.setSeedHandler(sH);
		hhU.fillHelix(0, energy.size1() - 1, 0, energy.size2() - 1);


		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////   FILLHELIXSEED  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		REQUIRE(hhU.fillHelixSeed(0, energy.size1() - 1, 0, energy.size2() - 1) == 6);

		// (0,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(0, 0)) == -4);
		REQUIRE(hhU.getHelixSeedLength1(0, 0) == 7);
		REQUIRE(hhU.getHelixSeedLength2(0, 0) == 6);

		// (0,1)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(0, 1)) == -3);
		REQUIRE(hhU.getHelixSeedLength1(0, 1) == 6);
		REQUIRE(hhU.getHelixSeedLength2(0, 1) == 5);

		// (2,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(2, 0)) == -3);
		REQUIRE(hhU.getHelixSeedLength1(2, 0) == 5);
		REQUIRE(hhU.getHelixSeedLength2(2, 0) == 4);

		// (2,1)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(2, 1)) == -3);
		REQUIRE(hhU.getHelixSeedLength1(2, 1) == 5);
		REQUIRE(hhU.getHelixSeedLength2(2, 1) == 5);

		// (3,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(3, 0)) == -2);
		REQUIRE(hhU.getHelixSeedLength1(3, 0) == 4);
		REQUIRE(hhU.getHelixSeedLength2(3, 0) == 3);

		// (3,1)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(3, 1)) == -2);
		REQUIRE(hhU.getHelixSeedLength1(3, 1) == 4);
		REQUIRE(hhU.getHelixSeedLength2(3, 1) == 3);


		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		Interaction interaction(r1, r2);
		hhU.traceBackHelixSeed(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 3);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 2);
		REQUIRE(interaction.basePairs.begin()->second == 4);

		REQUIRE(interaction.basePairs.rbegin()->first == 5);
		REQUIRE(interaction.basePairs.rbegin()->second == 2);

		// Case (3,1)
		interaction.clear();
		hhU.traceBackHelixSeed(interaction, 3, 1);
		REQUIRE(interaction.basePairs.size() == 1);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 5);
		REQUIRE(interaction.basePairs.begin()->second == 3);

		REQUIRE(interaction.basePairs.rbegin()->first == 5);
		REQUIRE(interaction.basePairs.rbegin()->second == 3);

		// Case (1,1)
		interaction.clear();
		hhU.traceBackHelixSeed(interaction, 0, 1);
		REQUIRE(interaction.basePairs.size() == 2);
		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 2);
		REQUIRE(interaction.basePairs.begin()->second == 3);

		REQUIRE(interaction.basePairs.rbegin()->first == 3);
		REQUIRE(interaction.basePairs.rbegin()->second == 2);
	}

	SECTION("HelixSeed: Case 10 - seed surrounded by unpaired bases", "[HelixHandlerUnpaired]") {

		RnaSequence r1("r1", "GGAGGGAAGG");
		RnaSequence r2("r2", "CACCCCACAAC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 7, 2, 999, 0, false);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ
		SeedConstraint sC(3, 0, 0, 0, 0, AccessibilityDisabled::ED_UPPER_BOUND, 0, IndexRangeList(""), IndexRangeList(""),
						  "", false, false, false );

		SeedHandlerMfe sH(energy, sC);
		HelixHandlerUnpaired hhU(energy, hC);

		sH.fillSeed(0, energy.size1() - 1, 0, energy.size2() - 1);
		hhU.setSeedHandler(sH);
		hhU.fillHelix(0, energy.size1() - 1, 0, energy.size2() - 1);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////   FILLHELIXSEED  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		REQUIRE(hhU.fillHelixSeed(0, energy.size1() - 1, 0, energy.size2() - 1) == 6);

		// (0,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(0, 0)) == -6);
		REQUIRE(hhU.getHelixSeedLength1(0, 0) == 10);
		REQUIRE(hhU.getHelixSeedLength2(0, 0) == 11);

		// (0,3)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(0, 3)) == -4);
		REQUIRE(hhU.getHelixSeedLength1(0, 3) == 6);
		REQUIRE(hhU.getHelixSeedLength2(0, 3) == 6);

		// (1,3)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(1, 3)) == -5);
		REQUIRE(hhU.getHelixSeedLength1(1, 3) == 9);
		REQUIRE(hhU.getHelixSeedLength2(1, 3) == 8);

		// (1,5)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(1, 5)) == -3);
		REQUIRE(hhU.getHelixSeedLength1(1, 5) == 5);
		REQUIRE(hhU.getHelixSeedLength2(1, 5) == 4);

		// (3,5)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(3, 5)) == -4);
		REQUIRE(hhU.getHelixSeedLength1(3, 5) == 7);
		REQUIRE(hhU.getHelixSeedLength2(3, 5) == 6);

		// (3,6)
		REQUIRE(E_2_Ekcal(hhU.getHelixSeedE(3, 6)) == -2);
		REQUIRE(hhU.getHelixSeedLength1(3, 6) == 3);
		REQUIRE(hhU.getHelixSeedLength2(3, 6) == 3);


		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		Interaction interaction(r1, r2);
		hhU.traceBackHelixSeed(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 5);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 1);
		REQUIRE(interaction.basePairs.begin()->second == 7);

		REQUIRE(interaction.basePairs.rbegin()->first == 8);
		REQUIRE(interaction.basePairs.rbegin()->second == 2);

	}
}
