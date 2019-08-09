
#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/RnaSequence.h"
#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/HelixConstraint.h"
#include "IntaRNA/HelixHandlerNoBulgeMax.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/Interaction.h"
#include "IntaRNA/ReverseAccessibility.h"


using namespace IntaRNA;

TEST_CASE( "HelixHandlerNoBulgeMax", "[HelixHandlerNoBulgeMax]") {


	SECTION("getter", "[HelixHandlerNoBulgeMax]") {

		RnaSequence r1("r1", "GGGGG");
		RnaSequence r2("r2", "CCCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 10, 0, 999, 0, false);
		HelixHandlerNoBulgeMax hh(energy, hC);

		REQUIRE(&hh.getInteractionEnergy() == &energy);
		REQUIRE(&hh.getConstraint() == &hC);

	}

	SECTION("Helix: Case 1 - Everything is complementary", "[HelixHandlerNoBulgeMax]") {

		// Case 1 Perfect sequence
		const RnaSequence r1("r1", "GGGGG");
		const RnaSequence r2("r2", "CCCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		const HelixConstraint hC(2, 4, 0, 999, 0, false);
		HelixHandlerNoBulgeMax hh(energy, hC);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   FILLHELIX  ////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// When counting all non-inf values for helices -> 29
		// When only counting the optimal helices that are non-inf -> 16
		// (0,0) -> 4 ; (1,0) -> 4 ; (2,0) -> 3 ; (3,0) -> 2
		// (0,1) -> 4 ; (1,1) -> 4 ; (2,1) -> 3 ; (3,1) -> 2
		// (0,2) -> 3 ; (1,2) -> 3 ; (2,2) -> 3 ; (3,2) -> 2
		// (0,3) -> 2 ; (1,3) -> 2 ; (2,3) -> 2 ; (3,3) -> 2
		REQUIRE(hh.fillHelix(0, energy.size1() - 1, 0, energy.size2() - 1) == 16);

		// All optimal combinations

		// (0,0)
		REQUIRE(E_2_Ekcal(hh.getHelixE(0, 0)) == (-3));
		REQUIRE(hh.getHelixLength1(0, 0) == 4);
		REQUIRE(hh.getHelixLength2(0, 0) == hh.getHelixLength1(0, 0));

		// (0,1)
		REQUIRE(E_2_Ekcal(hh.getHelixE(0, 1)) == (-3));
		REQUIRE(hh.getHelixLength1(0, 1) == 4);
		REQUIRE(hh.getHelixLength2(0, 1) == hh.getHelixLength1(0, 1));

		// (0,2)
		REQUIRE(E_2_Ekcal(hh.getHelixE(0, 2)) == (-2));
		REQUIRE(hh.getHelixLength1(0, 2) == 3);
		REQUIRE(hh.getHelixLength2(0, 2) == hh.getHelixLength1(0, 2));

		// (0,3)
		REQUIRE(E_2_Ekcal(hh.getHelixE(0, 3)) == (-1));
		REQUIRE(hh.getHelixLength1(0, 3) == 2);
		REQUIRE(hh.getHelixLength2(0, 3) == hh.getHelixLength1(0, 3));

		// (1,0)
		REQUIRE(E_2_Ekcal(hh.getHelixE(1, 0)) == (-3));
		REQUIRE(hh.getHelixLength1(1, 0) == 4);
		REQUIRE(hh.getHelixLength2(1, 0) == hh.getHelixLength1(1, 0));

		// (1,1)
		REQUIRE(E_2_Ekcal(hh.getHelixE(1, 1)) == (-3));
		REQUIRE(hh.getHelixLength1(1, 1) == 4);
		REQUIRE(hh.getHelixLength2(1, 1) == hh.getHelixLength1(1, 1));

		// (1,2)
		REQUIRE(E_2_Ekcal(hh.getHelixE(1, 2)) == (-2));
		REQUIRE(hh.getHelixLength1(1, 2) == 3);
		REQUIRE(hh.getHelixLength2(1, 2) == hh.getHelixLength1(1, 2));

		// (1,3)
		REQUIRE(E_2_Ekcal(hh.getHelixE(1, 3)) == (-1));
		REQUIRE(hh.getHelixLength1(1, 3) == 2);
		REQUIRE(hh.getHelixLength2(1, 3) == hh.getHelixLength1(1, 3));

		// (2,0)
		REQUIRE(E_2_Ekcal(hh.getHelixE(2, 0)) == (-2));
		REQUIRE(hh.getHelixLength1(2, 0) == 3);
		REQUIRE(hh.getHelixLength2(2, 0) == hh.getHelixLength1(2, 0));

		// (2,1)
		REQUIRE(E_2_Ekcal(hh.getHelixE(2, 1)) == (-2));
		REQUIRE(hh.getHelixLength1(2, 1) == 3);
		REQUIRE(hh.getHelixLength2(2, 1) == hh.getHelixLength1(2, 1));

		// (2,2)
		REQUIRE(E_2_Ekcal(hh.getHelixE(2, 2)) == (-2));
		REQUIRE(hh.getHelixLength1(2, 2) == 3);
		REQUIRE(hh.getHelixLength2(2, 2) == hh.getHelixLength1(2, 2));

		// (2,3)
		REQUIRE(E_2_Ekcal(hh.getHelixE(2, 3)) == (-1));
		REQUIRE(hh.getHelixLength1(2, 3) == 2);
		REQUIRE(hh.getHelixLength2(2, 3) == hh.getHelixLength1(2, 3));

		// (3,0)
		REQUIRE(E_2_Ekcal(hh.getHelixE(3, 0)) == (-1));
		REQUIRE(hh.getHelixLength1(3, 0) == 2);
		REQUIRE(hh.getHelixLength2(3, 0) == hh.getHelixLength1(3, 0));

		// (3,1)
		REQUIRE(E_2_Ekcal(hh.getHelixE(3, 1)) == (-1));
		REQUIRE(hh.getHelixLength1(3, 1) == 2);
		REQUIRE(hh.getHelixLength2(3, 1) == hh.getHelixLength1(3, 1));

		// (3,2)
		REQUIRE(E_2_Ekcal(hh.getHelixE(3, 2)) == (-1));
		REQUIRE(hh.getHelixLength1(3, 2) == 2);
		REQUIRE(hh.getHelixLength2(3, 2) == hh.getHelixLength1(3, 2));

		// (3,3)
		REQUIRE(E_2_Ekcal(hh.getHelixE(3, 3)) == -1);
		REQUIRE(hh.getHelixLength1(3, 3) == 2);
		REQUIRE(hh.getHelixLength2(3, 3) == hh.getHelixLength1(3, 3));


		// Several cases that do not allow a helix

		// (4,4)
		REQUIRE(hh.getHelixE(4, 4) == E_INF);
		REQUIRE(hh.getHelixLength1(4, 4) == 0);
		REQUIRE(hh.getHelixLength2(4, 4) == hh.getHelixLength1(4, 4));

		// (4,3)
		REQUIRE(hh.getHelixE(0, 4) == E_INF);
		REQUIRE(hh.getHelixLength1(0, 4) == 0);
		REQUIRE(hh.getHelixLength2(0, 4) == hh.getHelixLength1(0, 4));


		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		//////////////////////
		Interaction interaction(r1,r2);

		hh.traceBackHelix(interaction, 0, 0);


		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 1);
		REQUIRE(interaction.basePairs.begin()->second == 3);

		REQUIRE(interaction.basePairs.rbegin()->first == 2);
		REQUIRE(interaction.basePairs.rbegin()->second == 2);



		// Case (2,1)
		//////////////////////
		interaction.clear();

		hh.traceBackHelix(interaction, 2, 1);

		REQUIRE(interaction.basePairs.size() == 1);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 3);
		REQUIRE(interaction.basePairs.begin()->second == 2);

		REQUIRE(interaction.basePairs.rbegin()->first == 3);
		REQUIRE(interaction.basePairs.rbegin()->second == 2);
	}


	SECTION("Helix: Case 2 - Sequence 1 contains an A", "[HelixHandlerNoBulgeMax]") {
		// Case 2 - sequence containing an "A" to disrupt perfect stacking
		RnaSequence r1("r1", "GGGAGG");
		RnaSequence r2("r2", "CCCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 0, 999, 0, false);
		HelixHandlerNoBulgeMax hh(energy, hC);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   FILLHELIX  ////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		REQUIRE(hh.fillHelix(0, energy.size1() - 1, 0, energy.size2() - 1) == 12);


		// All optimal combinations
		// (0,0)
		REQUIRE(E_2_Ekcal(hh.getHelixE(0, 0)) == -2);
		REQUIRE(hh.getHelixLength1(0, 0) == 3);
		REQUIRE(hh.getHelixLength2(0, 0) == hh.getHelixLength1(0, 0));

		// (0,1)
		REQUIRE(E_2_Ekcal(hh.getHelixE(0, 1)) == -2);
		REQUIRE(hh.getHelixLength1(0, 1) == 3);
		REQUIRE(hh.getHelixLength2(0, 1) == hh.getHelixLength1(0, 1));

		// (0,2)
		REQUIRE(E_2_Ekcal(hh.getHelixE(0, 2)) == -2);
		REQUIRE(hh.getHelixLength1(0, 2) == 3);
		REQUIRE(hh.getHelixLength2(0, 2) == hh.getHelixLength1(0, 2));

		// (0,3)
		REQUIRE(E_2_Ekcal(hh.getHelixE(0, 3)) == -1);
		REQUIRE(hh.getHelixLength1(0, 3) == 2);
		REQUIRE(hh.getHelixLength2(0, 3) == hh.getHelixLength1(0, 3));

		// (1,0)
		REQUIRE(E_2_Ekcal(hh.getHelixE(1, 0)) == -1);
		REQUIRE(hh.getHelixLength1(1, 0) == 2);
		REQUIRE(hh.getHelixLength2(1, 0) == hh.getHelixLength1(1, 0));

		// (1,1)
		REQUIRE(E_2_Ekcal(hh.getHelixE(1, 1)) == -1);
		REQUIRE(hh.getHelixLength1(1, 1) == 2);
		REQUIRE(hh.getHelixLength2(1, 1) == hh.getHelixLength1(1, 1));

		// (1,2)
		REQUIRE(E_2_Ekcal(hh.getHelixE(1, 2)) == -1);
		REQUIRE(hh.getHelixLength1(1, 2) == 2);
		REQUIRE(hh.getHelixLength2(1, 2) == hh.getHelixLength1(1, 2));

		// (1,3)
		REQUIRE(E_2_Ekcal(hh.getHelixE(1, 3)) == -1);
		REQUIRE(hh.getHelixLength1(1, 3) == 2);
		REQUIRE(hh.getHelixLength2(1, 3) == hh.getHelixLength1(1, 3));

		// (4,0)
		REQUIRE(E_2_Ekcal(hh.getHelixE(4, 0)) == -1);
		REQUIRE(hh.getHelixLength1(4, 0) == 2);
		REQUIRE(hh.getHelixLength2(4, 0) == hh.getHelixLength1(4, 0));

		// (4,1)
		REQUIRE(E_2_Ekcal(hh.getHelixE(4, 1)) == -1);
		REQUIRE(hh.getHelixLength1(4, 1) == 2);
		REQUIRE(hh.getHelixLength2(4, 1) == hh.getHelixLength1(4, 1));

		// (4,2)
		REQUIRE(E_2_Ekcal(hh.getHelixE(4, 2)) == -1);
		REQUIRE(hh.getHelixLength1(4, 2) == 2);
		REQUIRE(hh.getHelixLength2(4, 2) == hh.getHelixLength1(4, 2));

		// (4,3)
		REQUIRE(E_2_Ekcal(hh.getHelixE(4, 3)) == -1);
		REQUIRE(hh.getHelixLength1(4, 3) == 2);
		REQUIRE(hh.getHelixLength2(4, 3) == hh.getHelixLength1(4, 3));


		// Not viable cases
		// (3,0) Not complementary
		REQUIRE(hh.getHelixE(3,0) == E_INF);
		REQUIRE(hh.getHelixLength1(3,0) == 0);
		REQUIRE(hh.getHelixLength2(3,0) == hh.getHelixLength1(3,0));

		// (5,4) no Helix possible (both sides too short)
		REQUIRE(hh.getHelixE(5,4) == E_INF);
		REQUIRE(hh.getHelixLength1(5,4) == 0);
		REQUIRE(hh.getHelixLength2(5,4) == hh.getHelixLength1(5,4));

		// (5,3) no helix possible
		REQUIRE(hh.getHelixE(5,3) == E_INF);
		REQUIRE(hh.getHelixLength1(5,3) == 0);
		REQUIRE(hh.getHelixLength2(5,3) == hh.getHelixLength1(5,3));

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		//////////////////////
		Interaction interaction(r1,r2);

		hh.traceBackHelix(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 1);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 1);
		REQUIRE(interaction.basePairs.begin()->second == 3);

		REQUIRE(interaction.basePairs.rbegin()->first == 1);
		REQUIRE(interaction.basePairs.rbegin()->second == 3);

	}

	SECTION("Helix: Case 3 - A 'wall' of A's disrupts the possible helices", "[HelixHandlerNoBulgeMax]") {
		// Case 2 - sequence containing an "A"-wall to disrupt perfect stacking
		RnaSequence r1("r1", "GGGAAGG");
		RnaSequence r2("r2", "CCAACCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 0, 999, 0, false);
		HelixHandlerNoBulgeMax hh(energy, hC);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   FILLHELIX  ////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		REQUIRE(hh.fillHelix(0, energy.size1() - 1, 0, energy.size2() - 1) == 9);

		REQUIRE_FALSE(energy.areComplementary(5,4));
		// All optimal combinations
		// (0,0)
		REQUIRE(E_2_Ekcal(hh.getHelixE(0, 0)) == -2);
		REQUIRE(hh.getHelixLength1(0, 0) == 3);
		REQUIRE(hh.getHelixLength2(0, 0) == hh.getHelixLength1(0, 0));

		// (0,1)
		REQUIRE(E_2_Ekcal(hh.getHelixE(0, 1)) == -1);
		REQUIRE(hh.getHelixLength1(0, 1) == 2);
		REQUIRE(hh.getHelixLength2(0, 1) == hh.getHelixLength1(0, 1));

		// (0,5)
		REQUIRE(E_2_Ekcal(hh.getHelixE(0, 5)) == -1);
		REQUIRE(hh.getHelixLength1(0, 5) == 2);
		REQUIRE(hh.getHelixLength2(0, 5) == hh.getHelixLength1(0, 5));

		// (1,0)
		REQUIRE(E_2_Ekcal(hh.getHelixE(1, 0)) == -1);
		REQUIRE(hh.getHelixLength1(1, 0) == 2);
		REQUIRE(hh.getHelixLength2(1, 0) == hh.getHelixLength1(1, 0));

		// (1,1)
		REQUIRE(E_2_Ekcal(hh.getHelixE(1, 1)) == -1);
		REQUIRE(hh.getHelixLength1(1, 1) == 2);
		REQUIRE(hh.getHelixLength2(1, 1) == hh.getHelixLength1(1, 1));

		// (1,5)
		REQUIRE(E_2_Ekcal(hh.getHelixE(1, 5)) == -1);
		REQUIRE(hh.getHelixLength1(1, 5) == 2);
		REQUIRE(hh.getHelixLength2(1, 5) == hh.getHelixLength1(1, 5));

		// (5,0)
		REQUIRE(E_2_Ekcal(hh.getHelixE(5, 0)) == -1);
		REQUIRE(hh.getHelixLength1(5, 0) == 2);
		REQUIRE(hh.getHelixLength2(5, 0) == hh.getHelixLength1(5, 0));

		// (5,1)
		REQUIRE(E_2_Ekcal(hh.getHelixE(5, 1)) == -1);
		REQUIRE(hh.getHelixLength1(5, 1) == 2);
		REQUIRE(hh.getHelixLength2(5, 1) == hh.getHelixLength1(5, 1));

		// (5,5)
		REQUIRE(E_2_Ekcal(hh.getHelixE(5, 5)) == -1);
		REQUIRE(hh.getHelixLength1(5, 5) == 2);
		REQUIRE(hh.getHelixLength2(5, 5) == hh.getHelixLength1(5, 5));



		// Not viable cases
		// (3,0) Not complementary
		REQUIRE(hh.getHelixE(3,0) == E_INF);
		REQUIRE(hh.getHelixLength1(3,0) == 0);
		REQUIRE(hh.getHelixLength2(3,0) == hh.getHelixLength1(3,0));

		// (5,4) no Helix possible (both sides too short)
		REQUIRE(hh.getHelixE(2,2) == E_INF);
		REQUIRE(hh.getHelixLength1(2,2) == 0);
		REQUIRE(hh.getHelixLength2(2,2) == hh.getHelixLength1(2,2));

		// (5,3) no helix possible
		REQUIRE(hh.getHelixE(5,3) == E_INF);
		REQUIRE(hh.getHelixLength1(5,3) == 0);
		REQUIRE(hh.getHelixLength2(5,3) == hh.getHelixLength1(5,3));

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		//////////////////////
		Interaction interaction(r1,r2);
		hh.traceBackHelix(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 1);
		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 1);
		REQUIRE(interaction.basePairs.begin()->second == 5);

		REQUIRE(interaction.basePairs.rbegin()->first == 1);
		REQUIRE(interaction.basePairs.rbegin()->second == 5);

		// Case (5,5) - Possible but only 2 base pairs long (e.g no bp needs to be reported)
		//////////////////////

		interaction.clear();
		hh.traceBackHelix(interaction, 5, 5);

		REQUIRE(interaction.basePairs.size() == 0);

	}

	SECTION("Helix: Case 4 - No interaction possible", "[HelixHandlerNoBulgeMax]") {
		// Case 4 -NO HELIX POSSIBLE
		RnaSequence r1("r1", "AAAAAAA");
		RnaSequence r2("r2", "AAAAAAA");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 0, 999, 0, false);
		HelixHandlerNoBulgeMax hh(energy, hC);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   FILLHELIX  ////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		REQUIRE(hh.fillHelix(0, energy.size1() - 1, 0, energy.size2() - 1) == 0);

		// NO POSSIBLE HELICES
		// (2,2)
		REQUIRE(hh.getHelixE(2, 2) == E_INF);
		REQUIRE(hh.getHelixLength1(2, 2) == 0);
		REQUIRE(hh.getHelixLength2(2, 2) == hh.getHelixLength1(2, 2));

		// (0,3)
		REQUIRE(hh.getHelixE(0, 3) == E_INF);
		REQUIRE(hh.getHelixLength1(0, 3) == 0);
		REQUIRE(hh.getHelixLength2(0, 3) == hh.getHelixLength1(0, 3));

	}

	SECTION("Helix: Case 5 - Example from LimStackHeuristic test", "[HelixHandlerNoBulgeMax]") {
		// Case 5 (3 length helix at the end)
		RnaSequence r1("r1", "gggaaggg");
		RnaSequence r2("r2", "cccaaccc");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 0, 999, 0, false);
		HelixHandlerNoBulgeMax hh(energy, hC);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   FILLHELIX  ////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		REQUIRE(hh.fillHelix(0, energy.size1() - 1, 0, energy.size2() - 1) == 16);

		// Possible helices
		// (5,5)
		REQUIRE(E_2_Ekcal(hh.getHelixE(5, 5)) == -2);
		REQUIRE(hh.getHelixLength1(5, 5) == 3);
		REQUIRE(hh.getHelixLength2(5, 5) == hh.getHelixLength1(5, 5));

		// (0,0)
		REQUIRE(E_2_Ekcal(hh.getHelixE(0, 0)) == -2);
		REQUIRE(hh.getHelixLength1(0, 0) == 3);
		REQUIRE(hh.getHelixLength2(0, 0) == hh.getHelixLength1(0, 0));

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		//////////////////////
		Interaction interaction(r1,r2);

		hh.traceBackHelix(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 1);
		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 1);
		REQUIRE(interaction.basePairs.begin()->second == 6);

		REQUIRE(interaction.basePairs.rbegin()->first == 1);
		REQUIRE(interaction.basePairs.rbegin()->second == 6);

		// Case (5, 5)
		//////////////////////
		interaction.clear();

		hh.traceBackHelix(interaction, 5, 5);

		REQUIRE(interaction.basePairs.size() == 1);
		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 6);
		REQUIRE(interaction.basePairs.begin()->second == 1);

		REQUIRE(interaction.basePairs.rbegin()->first == 6);
		REQUIRE(interaction.basePairs.rbegin()->second == 1);
	}

	SECTION("Helix: Case 6 - Special case Helix+E_init() == Helix+IL+H", "[HelixHandlerNoBulgeMax]") {
		// Case 6
		RnaSequence r1("r1", "GGUUGAAUUACGACAG");
		RnaSequence r2("r2", "cugaaaaacauaacc");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 10, 0, 999, 0, false);
		HelixHandlerNoBulgeMax hh(energy, hC);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   FILLHELIX  ////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		REQUIRE(hh.fillHelix(0, energy.size1() - 1, 0, energy.size2() - 1) == 26);

		// (0,0)
		REQUIRE(E_2_Ekcal(hh.getHelixE(0,0)) == -4);
		REQUIRE(hh.getHelixLength1(0,0) == 5);
		REQUIRE(hh.getHelixLength2(0,0) == hh.getHelixLength1(0,0));

		// (3,5)
		REQUIRE(E_2_Ekcal(hh.getHelixE(3,5)) == -1);
		REQUIRE(hh.getHelixLength1(3,5) == 2);
		REQUIRE(hh.getHelixLength2(3,5) == hh.getHelixLength1(3,5));

		// (7,11)
		REQUIRE(E_2_Ekcal(hh.getHelixE(7,11)) == -2);
		REQUIRE(hh.getHelixLength1(7,11) == 3);
		REQUIRE(hh.getHelixLength2(7,11) == hh.getHelixLength1(7,11));

		// (13,12)
		REQUIRE(E_2_Ekcal(hh.getHelixE(13,12)) == -2);
		REQUIRE(hh.getHelixLength1(13,12) == 3);
		REQUIRE(hh.getHelixLength2(13,12) == hh.getHelixLength1(13,12));


		// Not working

		// (7,11)
		REQUIRE(hh.getHelixE(7,12) == E_INF);
		REQUIRE(hh.getHelixLength1(7,12) == 0);
		REQUIRE(hh.getHelixLength2(7,12) == hh.getHelixLength1(7,12));

		// (0,1)
		REQUIRE(hh.getHelixE(0,1) == E_INF);
		REQUIRE(hh.getHelixLength1(0,1) == 0);
		REQUIRE(hh.getHelixLength2(0,1) == hh.getHelixLength1(0,1));

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		//////////////////////
		Interaction interaction(r1,r2);

		hh.traceBackHelix(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 3);
		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 1);
		REQUIRE(interaction.basePairs.begin()->second == 13);

		REQUIRE(interaction.basePairs.rbegin()->first == 3);
		REQUIRE(interaction.basePairs.rbegin()->second == 11);

	}
}
