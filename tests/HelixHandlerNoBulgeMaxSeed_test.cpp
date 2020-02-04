
#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/RnaSequence.h"
#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/HelixConstraint.h"
#include "IntaRNA/HelixHandlerNoBulgeMax.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/SeedHandlerMfe.h"

using namespace IntaRNA;

TEST_CASE( "HelixHandlerNoBulgeMaxSeed", "[HelixHandlerNoBulgeMax]" ) {

	SECTION("HelixSeed: Case 1 - all complementary", "[HelixHandlerNoBulgeMax]") {

		RnaSequence r1("r1", "GGGGG");
		RnaSequence r2("r2", "CCCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2,4, 0, 999, 0, false);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ
		SeedConstraint sC(3,0,0,0,0
				, AccessibilityDisabled::ED_UPPER_BOUND, 0
				, IndexRangeList("")
				, IndexRangeList("")
				, "", false, false, false );

		SeedHandlerMfe sH(energy, sC);
		HelixHandlerNoBulgeMax hhS(energy, hC);

		sH.fillSeed(0, energy.size1()-1, 0,energy.size2()-1);
		hhS.setSeedHandler(sH);

		hhS.fillHelix(0, energy.size1()-1, 0, energy.size2()-1);
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////   FILLHELIXSEED  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		REQUIRE(hhS.fillHelixSeed( 0,energy.size1()-1, 0,energy.size2()-1 ) == 9);

		// (0,0)
		REQUIRE(E_2_Ekcal(hhS.getHelixSeedE(0,0)) == -3);
		REQUIRE(hhS.getHelixSeedLength1(0,0) == 4);
		REQUIRE(hhS.getHelixSeedLength2(0,0) == 4);

		// (0,1)
		REQUIRE(E_2_Ekcal(hhS.getHelixSeedE(0,1)) == -3);
		REQUIRE(hhS.getHelixSeedLength1(0,1) == 4);
		REQUIRE(hhS.getHelixSeedLength2(0,1) == 4);

		// (0,2)
		REQUIRE(E_2_Ekcal(hhS.getHelixSeedE(0,2)) == -2);
		REQUIRE(hhS.getHelixSeedLength1(0,2) == 3);
		REQUIRE(hhS.getHelixSeedLength2(0,2) == 3);

		// (0,3)
		REQUIRE(hhS.getHelixSeedE(0,3) == E_INF);
		REQUIRE(hhS.getHelixSeedLength1(0,3) == 0);
		REQUIRE(hhS.getHelixSeedLength2(0,3) == 0);

		// (1,1)
		REQUIRE(E_2_Ekcal(hhS.getHelixSeedE(1,1)) == -3);
		REQUIRE(hhS.getHelixSeedLength1(1,1) == 4);
		REQUIRE(hhS.getHelixSeedLength2(1,1) == 4);

		// (2,2)
		REQUIRE(E_2_Ekcal(hhS.getHelixSeedE(2,2)) == -2);
		REQUIRE(hhS.getHelixSeedLength1(2,2) == 3);
		REQUIRE(hhS.getHelixSeedLength2(2,2) == 3);

		// (4,4)
		REQUIRE(hhS.getHelixSeedE(4,4) == E_INF);
		REQUIRE(hhS.getHelixSeedLength1(4,4) == 0);
		REQUIRE(hhS.getHelixSeedLength2(4,4) == 0);


		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		Interaction interaction(r1,r2);
		hhS.traceBackHelixSeed(interaction,0,0);

		REQUIRE(interaction.basePairs.size() == 2);

		// First / last base pair of helixSeed
		REQUIRE(interaction.basePairs.begin()->first == 1);
		REQUIRE(interaction.basePairs.begin()->second == 3);

		REQUIRE(interaction.basePairs.rbegin()->first == 2);
		REQUIRE(interaction.basePairs.rbegin()->second == 2);

		// Case (4,4)
		interaction.clear();
		hhS.traceBackHelixSeed(interaction,4,4);

		REQUIRE(interaction.basePairs.size() == 0);
	}

	SECTION("HelixSeed: Case 2 - 'A' disrupting complementarity", "[HelixHandlerNoBulgeMax]") {

		RnaSequence r1("r1", "GGAGG");
		RnaSequence r2("r2", "CCACC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4,  0, 999, 0, false);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ / seedNoGU
		SeedConstraint sC(3, 0, 0, 0, 0, AccessibilityDisabled::ED_UPPER_BOUND, 0, IndexRangeList(""), IndexRangeList(""),
						  "", false, false, false);

		SeedHandlerMfe sH(energy, sC);
		HelixHandlerNoBulgeMax hhS(energy, hC);

		sH.fillSeed(0, energy.size1()-1, 0,energy.size2()-1);
		hhS.setSeedHandler(sH);
		hhS.fillHelix(0, energy.size1()-1, 0, energy.size2()-1);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////   FILLHELIXSEED  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		REQUIRE(hhS.fillHelixSeed(0, energy.size1() - 1, 0, energy.size2() - 1) == 0);

		// (0,0)
		REQUIRE(hhS.getHelixSeedE(0, 0) == E_INF);
		REQUIRE(hhS.getHelixSeedLength1(0, 0) == 0);
		REQUIRE(hhS.getHelixSeedLength2(0, 0) == 0);

		// (0,2)
		REQUIRE(hhS.getHelixSeedE(0, 2) == E_INF);
		REQUIRE(hhS.getHelixSeedLength1(0, 2) == 0);
		REQUIRE(hhS.getHelixSeedLength2(0, 2) == 0);

		// (1,3)
		REQUIRE(hhS.getHelixSeedE(1, 3) == E_INF);
		REQUIRE(hhS.getHelixSeedLength1(1, 3) == 0);
		REQUIRE(hhS.getHelixSeedLength2(1, 3) == 0);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		Interaction interaction(r1, r2);
		hhS.traceBackHelixSeed(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 0);
	}

	SECTION("HelixSeed: Case 3 - only seed possible", "[HelixHandlerNoBulgeMax]") {

		RnaSequence r1("r1", "AGGGA");
		RnaSequence r2("r2", "ACCCA");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4,  0, 999, 0, false);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ / seedNoGU
		SeedConstraint sC(3, 0, 0, 0, 0, AccessibilityDisabled::ED_UPPER_BOUND, 0, IndexRangeList(""), IndexRangeList(""),
						  "", false, false, false);

		SeedHandlerMfe sH(energy, sC);
		HelixHandlerNoBulgeMax hhS(energy, hC);

		sH.fillSeed(0, energy.size1()-1, 0,energy.size2()-1);
		hhS.setSeedHandler(sH);

		hhS.fillHelix(0, energy.size1()-1, 0, energy.size2()-1);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////   FILLHELIXSEED  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		REQUIRE(hhS.fillHelixSeed(0, energy.size1() - 1, 0, energy.size2() - 1) == 1);

		// (0,0)
		REQUIRE(hhS.getHelixSeedE(0, 0) == E_INF);
		REQUIRE(hhS.getHelixSeedLength1(0, 0) == 0);
		REQUIRE(hhS.getHelixSeedLength2(0, 0) == 0);

		// (0,2)
		REQUIRE(hhS.getHelixSeedE(1, 0) == E_INF);
		REQUIRE(hhS.getHelixSeedLength1(1, 0) == 0);
		REQUIRE(hhS.getHelixSeedLength2(1, 0) == 0);

		// (1,3)
		REQUIRE(E_2_Ekcal(hhS.getHelixSeedE(1, 1)) == -2);
		REQUIRE(hhS.getHelixSeedLength1(1, 1) == 3);
		REQUIRE(hhS.getHelixSeedLength2(1, 1) == 3);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		Interaction interaction(r1, r2);
		hhS.traceBackHelixSeed(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 0);

		// Case (1,1)
		interaction.clear();
		hhS.traceBackHelixSeed(interaction, 1, 1);
		REQUIRE(interaction.basePairs.size() == 1);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 2);
		REQUIRE(interaction.basePairs.begin()->second == 2);

		REQUIRE(interaction.basePairs.rbegin()->first == 2);
		REQUIRE(interaction.basePairs.rbegin()->second == 2);
	}

	SECTION("HelixSeed: Case4 - unpaired allowed in seed", "[HelixHandlerNoBulgeMax]") {

		RnaSequence r1("r1", "GGAGG");
		RnaSequence r2("r2", "CCACC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4,  0, 999, 0, false);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ / seedNoGU
		SeedConstraint sC(3, 2, 1, 1, 0, AccessibilityDisabled::ED_UPPER_BOUND, 0, IndexRangeList(""), IndexRangeList(""),
						  "", false, false, false);

		SeedHandlerMfe sH(energy, sC);
		HelixHandlerNoBulgeMax hhS(energy, hC);

		sH.fillSeed(0, energy.size1()-1, 0,energy.size2()-1);
		hhS.setSeedHandler(sH);
		hhS.fillHelix(0, energy.size1()-1, 0, energy.size2()-1);


		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////   FILLHELIXSEED  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		REQUIRE(hhS.fillHelixSeed(0, energy.size1() - 1, 0, energy.size2() - 1) == 4);

		// (0,0)
		REQUIRE(E_2_Ekcal(hhS.getHelixSeedE(0, 0)) == -3);
		REQUIRE(hhS.getHelixSeedLength1(0, 0) == 5);
		REQUIRE(hhS.getHelixSeedLength2(0, 0) == 5);

		// (0,2)
		REQUIRE(E_2_Ekcal(hhS.getHelixSeedE(1, 0)) == -2);
		REQUIRE(hhS.getHelixSeedLength1(1, 0) == 4);
		REQUIRE(hhS.getHelixSeedLength2(1, 0) == 4);

		// (1,3)
		REQUIRE(E_2_Ekcal(hhS.getHelixSeedE(1, 1)) == -2);
		REQUIRE(hhS.getHelixSeedLength1(1, 1) == 4);
		REQUIRE(hhS.getHelixSeedLength2(1, 1) == 4);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		Interaction interaction(r1, r2);
		hhS.traceBackHelixSeed(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 2);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 1);
		REQUIRE(interaction.basePairs.begin()->second == 3);

		REQUIRE(interaction.basePairs.rbegin()->first == 3);
		REQUIRE(interaction.basePairs.rbegin()->second == 1);

		// Case (1,1)
		interaction.clear();
		hhS.traceBackHelixSeed(interaction, 1, 1);
		REQUIRE(interaction.basePairs.size() == 1);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 3);
		REQUIRE(interaction.basePairs.begin()->second == 1);

		REQUIRE(interaction.basePairs.rbegin()->first == 3);
		REQUIRE(interaction.basePairs.rbegin()->second == 1);

		// Case (0,1)
		interaction.clear();
		hhS.traceBackHelixSeed(interaction, 0, 1);
		REQUIRE(interaction.basePairs.size() == 1);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 1);
		REQUIRE(interaction.basePairs.begin()->second == 1);

		REQUIRE(interaction.basePairs.rbegin()->first == 1);
		REQUIRE(interaction.basePairs.rbegin()->second == 1);

	}
}
