
#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/RnaSequence.h"
#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/HelixConstraint.h"
#include "IntaRNA/HelixHandlerUnpaired.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/Interaction.h"
#include "IntaRNA/ReverseAccessibility.h"


using namespace IntaRNA;

TEST_CASE( "HelixHandlerUnpaired", "[HelixHandlerUnpaired]") {


//	SECTION("getter", "[HelixHandlerUnpaired]") {
//
//		RnaSequence r1("r1", "GGGGG");
//		RnaSequence r2("r2", "CCCCC");
//		AccessibilityDisabled acc1(r1, 0, NULL);
//		AccessibilityDisabled acc2(r2, 0, NULL);
//		ReverseAccessibility racc(acc2);
//		InteractionEnergyBasePair energy(acc1, racc);
//
//		HelixConstraint hC(2, 10, 2, 999, 0, false);
//		HelixHandlerUnpaired hhU(energy, hC);
//
//		REQUIRE(&hhU.getInteractionEnergy() == &energy);
//		REQUIRE(&hhU.getConstraint() == &hC);
//
//	}
//
	SECTION("Helix: Case 1 - Everything is complementary", "[HelixHandlerUnpaired]") {

		// Case 1 Perfect sequence
		RnaSequence r1("r1", "GGGGG");
		RnaSequence r2("r2", "CCCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 2, 999, 0, false);
		HelixHandlerUnpaired hhU(energy, hC);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   FILLHELIX  ////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// When counting all non-inf values for helices -> 29
		// When only counting the optimal helices that are non-inf -> 16
		// (0,0) -> 4 ; (1,0) -> 4 ; (2,0) -> 3 ; (3,0) -> 2
		// (0,1) -> 4 ; (1,1) -> 4 ; (2,1) -> 3 ; (3,1) -> 2
		// (0,2) -> 3 ; (1,2) -> 3 ; (2,2) -> 3 ; (3,2) -> 2
		// (0,3) -> 2 ; (1,3) -> 2 ; (2,3) -> 2 ; (3,3) -> 2
		REQUIRE(hhU.fillHelix(0, energy.size1() - 1, 0, energy.size2() - 1) == 16);


		// All optimal combinations
		// (0,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(0, 0)) == -3);
		REQUIRE(hhU.getHelixLength1(0, 0) == 4);
		REQUIRE(hhU.getHelixLength2(0, 0) == hhU.getHelixLength1(0, 0));

		// (0,1)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(0, 1)) == -3);
		REQUIRE(hhU.getHelixLength1(0, 1) == 4);
		REQUIRE(hhU.getHelixLength2(0, 1) == hhU.getHelixLength1(0, 1));

		// (0,2)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(0, 2)) == -2);
		REQUIRE(hhU.getHelixLength1(0, 2) == 3);
		REQUIRE(hhU.getHelixLength2(0, 2) == hhU.getHelixLength1(0, 2));

		// (0,3)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(0, 3)) == -1);
		REQUIRE(hhU.getHelixLength1(0, 3) == 2);
		REQUIRE(hhU.getHelixLength2(0, 3) == hhU.getHelixLength1(0, 3));

		// (1,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(1, 0)) == -3);
		REQUIRE(hhU.getHelixLength1(1, 0) == 4);
		REQUIRE(hhU.getHelixLength2(1, 0) == hhU.getHelixLength1(1, 0));

		// (1,1)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(1, 1)) == -3);
		REQUIRE(hhU.getHelixLength1(1, 1) == 4);
		REQUIRE(hhU.getHelixLength2(1, 1) == hhU.getHelixLength1(1, 1));

		// (1,2)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(1, 2)) == -2);
		REQUIRE(hhU.getHelixLength1(1, 2) == 3);
		REQUIRE(hhU.getHelixLength2(1, 2) == hhU.getHelixLength1(1, 2));

		// (1,3)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(1, 3)) == -1);
		REQUIRE(hhU.getHelixLength1(1, 3) == 2);
		REQUIRE(hhU.getHelixLength2(1, 3) == hhU.getHelixLength1(1, 3));

		// (2,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(2, 0)) == -2);
		REQUIRE(hhU.getHelixLength1(2, 0) == 3);
		REQUIRE(hhU.getHelixLength2(2, 0) == hhU.getHelixLength1(2, 0));

		// (2,1)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(2, 1)) == -2);
		REQUIRE(hhU.getHelixLength1(2, 1) == 3);
		REQUIRE(hhU.getHelixLength2(2, 1) == hhU.getHelixLength1(2, 1));

		// (2,2)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(2, 2)) == -2);
		REQUIRE(hhU.getHelixLength1(2, 2) == 3);
		REQUIRE(hhU.getHelixLength2(2, 2) == hhU.getHelixLength1(2, 2));

		// (2,3)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(2, 3)) == -1);
		REQUIRE(hhU.getHelixLength1(2, 3) == 2);
		REQUIRE(hhU.getHelixLength2(2, 3) == hhU.getHelixLength1(2, 3));

		// (3,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(3, 0)) == -1);
		REQUIRE(hhU.getHelixLength1(3, 0) == 2);
		REQUIRE(hhU.getHelixLength2(3, 0) == hhU.getHelixLength1(3, 0));

		// (3,1)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(3, 1)) == -1);
		REQUIRE(hhU.getHelixLength1(3, 1) == 2);
		REQUIRE(hhU.getHelixLength2(3, 1) == hhU.getHelixLength1(3, 1));

		// (3,2)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(3, 2)) == -1);
		REQUIRE(hhU.getHelixLength1(3, 2) == 2);
		REQUIRE(hhU.getHelixLength2(3, 2) == hhU.getHelixLength1(3, 2));

		// (3,3)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(3, 3)) == -1);
		REQUIRE(hhU.getHelixLength1(3, 3) == 2);
		REQUIRE(hhU.getHelixLength2(3, 3) == hhU.getHelixLength1(3, 3));


		// Several cases that do not allow a helix

		// (4,4)
		REQUIRE(hhU.getHelixE(4, 4) == E_INF);
		REQUIRE(hhU.getHelixLength1(4, 4) == 0);
		REQUIRE(hhU.getHelixLength2(4, 4) == hhU.getHelixLength1(4, 4));

		// (4,3)
		REQUIRE(hhU.getHelixE(0, 4) == E_INF);
		REQUIRE(hhU.getHelixLength1(0, 4) == 0);
		REQUIRE(hhU.getHelixLength2(0, 4) == hhU.getHelixLength1(0, 4));


		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		//////////////////////
		Interaction interaction(r1, r2);

		hhU.traceBackHelix(interaction, 0, 0);


		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 1);
		REQUIRE(interaction.basePairs.begin()->second == 3);

		REQUIRE(interaction.basePairs.rbegin()->first == 2);
		REQUIRE(interaction.basePairs.rbegin()->second == 2);



		// Case (2,1)
		//////////////////////
		interaction.clear();

		hhU.traceBackHelix(interaction, 2, 1);

		REQUIRE(interaction.basePairs.size() == 1);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 3);
		REQUIRE(interaction.basePairs.begin()->second == 2);

		REQUIRE(interaction.basePairs.rbegin()->first == 3);
		REQUIRE(interaction.basePairs.rbegin()->second == 2);
	}

	SECTION("Helix: Case 2 - Sequence 1 contains an A", "[HelixHandlerUnpaired]") {
		// Case 2 - sequence containing an "A" to disrupt perfect stacking
		RnaSequence r1("r1", "GGGAGG");
		RnaSequence r2("r2", "CCCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 2, 999, 0, false);
		HelixHandlerUnpaired hhU(energy, hC);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   FILLHELIX  ////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		REQUIRE(hhU.fillHelix(0, energy.size1() - 1, 0, energy.size2() - 1) == 16);


		// All optimal combinations
		// (0,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(0, 0)) == -3);
		REQUIRE(hhU.getHelixLength1(0, 0) == 5);
		REQUIRE(hhU.getHelixLength2(0, 0) == 4);

		// (0,1)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(0, 1)) == -3);
		REQUIRE(hhU.getHelixLength1(0, 1) == 5);
		REQUIRE(hhU.getHelixLength2(0, 1) == 4);

		// (0,2)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(0, 2)) == -2);
		REQUIRE(hhU.getHelixLength1(0, 2) == 3);
		REQUIRE(hhU.getHelixLength2(0, 2) == hhU.getHelixLength1(0, 2));

		// (0,3)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(0, 3)) == -1);
		REQUIRE(hhU.getHelixLength1(0, 3) == 2);
		REQUIRE(hhU.getHelixLength2(0, 3) == hhU.getHelixLength1(0, 3));

		// (1,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(1, 0)) == -3);
		REQUIRE(hhU.getHelixLength1(1, 0) == 5);
		REQUIRE(hhU.getHelixLength2(1, 0) == 4);

		// (1,1)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(1, 1)) == -3);
		REQUIRE(hhU.getHelixLength1(1, 1) == 5);
		REQUIRE(hhU.getHelixLength2(1, 1) == 4);

		// (1,2)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(1, 2)) == -2);
		REQUIRE(hhU.getHelixLength1(1, 2) == 4);
		REQUIRE(hhU.getHelixLength2(1, 2) == 3);

		// (1,3)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(1, 3)) == -1);
		REQUIRE(hhU.getHelixLength1(1, 3) == 2);
		REQUIRE(hhU.getHelixLength2(1, 3) == hhU.getHelixLength1(1, 3));

		// (2,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(2, 0)) == -2);
		REQUIRE(hhU.getHelixLength1(2, 0) == 4);
		REQUIRE(hhU.getHelixLength2(2, 0) == 3);

		// (2,1)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(2, 1)) == -2);
		REQUIRE(hhU.getHelixLength1(2, 1) == 4);
		REQUIRE(hhU.getHelixLength2(2, 1) == 3);

		// (2,2)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(2, 2)) == -2);
		REQUIRE(hhU.getHelixLength1(2, 2) == 4);
		REQUIRE(hhU.getHelixLength2(2, 2) == 3);

		// (2,3)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(2, 3)) == -1);
		REQUIRE(hhU.getHelixLength1(2, 3) == 3);
		REQUIRE(hhU.getHelixLength2(2, 3) == 2);

		// (4,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(4, 0)) == -1);
		REQUIRE(hhU.getHelixLength1(4, 0) == 2);
		REQUIRE(hhU.getHelixLength2(4, 0) == hhU.getHelixLength1(4, 0));

		// (4,1)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(4, 1)) == -1);
		REQUIRE(hhU.getHelixLength1(4, 1) == 2);
		REQUIRE(hhU.getHelixLength2(4, 1) == hhU.getHelixLength1(4, 1));

		// (4,2)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(4, 2)) == -1);
		REQUIRE(hhU.getHelixLength1(4, 2) == 2);
		REQUIRE(hhU.getHelixLength2(4, 2) == hhU.getHelixLength1(4, 2));

		// (4,3)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(4, 3)) == -1);
		REQUIRE(hhU.getHelixLength1(4, 3) == 2);
		REQUIRE(hhU.getHelixLength2(4, 3) == hhU.getHelixLength1(4, 3));


		// Not viable cases
		// (3,0) Not complementary
		REQUIRE(hhU.getHelixE(3,0) == E_INF);
		REQUIRE(hhU.getHelixLength1(3,0) == 0);
		REQUIRE(hhU.getHelixLength2(3,0) == hhU.getHelixLength1(3,0));

		// (5,4) no Helix possible (both sides too short)
		REQUIRE(hhU.getHelixE(5,4) == E_INF);
		REQUIRE(hhU.getHelixLength1(5,4) == 0);
		REQUIRE(hhU.getHelixLength2(5,4) == hhU.getHelixLength1(5,4));

		// (5,3) no helix possible
		REQUIRE(hhU.getHelixE(5,3) == E_INF);
		REQUIRE(hhU.getHelixLength1(5,3) == 0);
		REQUIRE(hhU.getHelixLength2(5,3) == hhU.getHelixLength1(5,3));

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		//////////////////////
		Interaction interaction(r1,r2);

		hhU.traceBackHelix(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 2);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 1);
		REQUIRE(interaction.basePairs.begin()->second == 3);

		REQUIRE(interaction.basePairs.rbegin()->first == 2);
		REQUIRE(interaction.basePairs.rbegin()->second == 2);

		// Case (2,1)
		//////////////////////
		interaction.clear();

		hhU.traceBackHelix(interaction, 2, 1);

		REQUIRE(interaction.basePairs.size() == 1);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 4);
		REQUIRE(interaction.basePairs.begin()->second == 2);

		REQUIRE(interaction.basePairs.rbegin()->first == 4);
		REQUIRE(interaction.basePairs.rbegin()->second == 2);

	}

	SECTION("Helix: Case 3 - A 'wall' of A's disrupts the possible helices", "[HelixHandlerUnpaired]") {
		// Case 2 - sequence containing an "A"-wall to disrupt perfect stacking
		RnaSequence r1("r1", "GGGAAGG");
		RnaSequence r2("r2", "CCAACCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 2, 999, 0, false);
		HelixHandlerUnpaired hhU(energy, hC);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   FILLHELIX  ////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		REQUIRE(hhU.fillHelix(0, energy.size1() - 1, 0, energy.size2() - 1) == 15);

		REQUIRE_FALSE(energy.areComplementary(5,4));
		// All optimal combinations
		// (0,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(0, 0)) == -3);
		REQUIRE(hhU.getHelixLength1(0, 0) == 6);
		REQUIRE(hhU.getHelixLength2(0, 0) == 7);

		// (0,1)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(0, 1)) == -3);
		REQUIRE(hhU.getHelixLength1(0, 1) == 6);
		REQUIRE(hhU.getHelixLength2(0, 1) == 6);

		// (0,2)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(0, 2)) == -2);
		REQUIRE(hhU.getHelixLength1(0, 2) == 3);
		REQUIRE(hhU.getHelixLength2(0, 2) == 5);

		// (0,5)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(0, 5)) == -1);
		REQUIRE(hhU.getHelixLength1(0, 5) == 2);
		REQUIRE(hhU.getHelixLength2(0, 5) == hhU.getHelixLength1(0, 5));

		// (1,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(1, 0)) == -3);
		REQUIRE(hhU.getHelixLength1(1, 0) == 6);
		REQUIRE(hhU.getHelixLength2(1, 0) == 6);

		// (1,1)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(1, 1)) == -1);
		REQUIRE(hhU.getHelixLength1(1, 1) == 2);
		REQUIRE(hhU.getHelixLength2(1, 1) == hhU.getHelixLength1(1, 1));

		// (1,2)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(1, 2)) == -2);
		REQUIRE(hhU.getHelixLength1(1, 2) == 5);
		REQUIRE(hhU.getHelixLength2(1, 2) == 5);

		// (1,5)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(1, 5)) == -1);
		REQUIRE(hhU.getHelixLength1(1, 5) == 2);
		REQUIRE(hhU.getHelixLength2(1, 5) == hhU.getHelixLength1(1, 5));

		// (2,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(2, 0)) == -2);
		REQUIRE(hhU.getHelixLength1(2, 0) == 5);
		REQUIRE(hhU.getHelixLength2(2, 0) == 3);

		// (2,1)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(2, 1)) == -2);
		REQUIRE(hhU.getHelixLength1(2, 1) == 5);
		REQUIRE(hhU.getHelixLength2(2, 1) == 5);

		// (2,5)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(2, 5)) == -1);
		REQUIRE(hhU.getHelixLength1(2, 5) == 4);
		REQUIRE(hhU.getHelixLength2(2, 5) == 2);

		// (5,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(5, 0)) == -1);
		REQUIRE(hhU.getHelixLength1(5, 0) == 2);
		REQUIRE(hhU.getHelixLength2(5, 0) == hhU.getHelixLength1(5, 0));

		// (5,1)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(5, 1)) == -1);
		REQUIRE(hhU.getHelixLength1(5, 1) == 2);
		REQUIRE(hhU.getHelixLength2(5, 1) == hhU.getHelixLength1(5, 1));

		// (5,2)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(5, 2)) == -1);
		REQUIRE(hhU.getHelixLength1(5, 2) == 2);
		REQUIRE(hhU.getHelixLength2(5, 2) == 4);

		// (5,5)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(5, 5)) == -1);
		REQUIRE(hhU.getHelixLength1(5, 5) == 2);
		REQUIRE(hhU.getHelixLength2(5, 5) == hhU.getHelixLength1(5, 5));



		// Not viable cases
		// (2,2) Not possible. Too many unpaired bases.
		REQUIRE(hhU.getHelixE(2, 2) == E_INF);
		REQUIRE(hhU.getHelixLength1(2, 2) == 0);
		REQUIRE(hhU.getHelixLength2(2, 2) == 0);

		// (3,0) Not complementary
		REQUIRE(hhU.getHelixE(3,0) == E_INF);
		REQUIRE(hhU.getHelixLength1(3,0) == 0);
		REQUIRE(hhU.getHelixLength2(3,0) == hhU.getHelixLength1(3,0));

		// (5,4) no Helix possible (both sides too short)
		REQUIRE(hhU.getHelixE(2,2) == E_INF);
		REQUIRE(hhU.getHelixLength1(2,2) == 0);
		REQUIRE(hhU.getHelixLength2(2,2) == hhU.getHelixLength1(2,2));

		// (5,3) no helix possible
		REQUIRE(hhU.getHelixE(5,3) == E_INF);
		REQUIRE(hhU.getHelixLength1(5,3) == 0);
		REQUIRE(hhU.getHelixLength2(5,3) == hhU.getHelixLength1(5,3));

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		//////////////////////
		Interaction interaction(r1,r2);
		hhU.traceBackHelix(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 2);
		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 1);
		REQUIRE(interaction.basePairs.begin()->second == 4);

		REQUIRE(interaction.basePairs.rbegin()->first == 2);
		REQUIRE(interaction.basePairs.rbegin()->second == 1);

		// Case (2,1) - two unpaired bases in sequence 1
		//////////////////////

		interaction.clear();
		hhU.traceBackHelix(interaction, 2, 1);

		REQUIRE(interaction.basePairs.size() == 1);
		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 5);
		REQUIRE(interaction.basePairs.begin()->second == 4);

		REQUIRE(interaction.basePairs.rbegin()->first == 5);
		REQUIRE(interaction.basePairs.rbegin()->second == 4);

		// Case (5,5) - Possible but only 2 base pairs long (e.g no bp needs to be reported)
		//////////////////////

		interaction.clear();
		hhU.traceBackHelix(interaction, 5, 5);

		REQUIRE(interaction.basePairs.size() == 0);


		// Case (2,0)
		/////////////////////
		interaction.clear();
		hhU.traceBackHelix(interaction, 2, 0);

		REQUIRE(interaction.basePairs.size() == 1);
		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 5);
		REQUIRE(interaction.basePairs.begin()->second == 5);

		REQUIRE(interaction.basePairs.rbegin()->first == 5);
		REQUIRE(interaction.basePairs.rbegin()->second == 5);

	}

	SECTION("Helix: Case 4 - No interaction possible", "[HelixHandlerUnpaired]") {
		// Case 4 -NO HELIX POSSIBLE
		RnaSequence r1("r1", "AAAAAAA");
		RnaSequence r2("r2", "AAAAAAA");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 2, 999, 0, false);
		HelixHandlerUnpaired hhU(energy, hC);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   FILLHELIX  ////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		REQUIRE(hhU.fillHelix(0, energy.size1() - 1, 0, energy.size2() - 1) == 0);

		// NO POSSIBLE HELICES
		// (2,2)
		REQUIRE(hhU.getHelixE(2, 2) == E_INF);
		REQUIRE(hhU.getHelixLength1(2, 2) == 0);
		REQUIRE(hhU.getHelixLength2(2, 2) == hhU.getHelixLength1(2, 2));

		// (0,3)
		REQUIRE(hhU.getHelixE(0, 3) == E_INF);
		REQUIRE(hhU.getHelixLength1(0, 3) == 0);
		REQUIRE(hhU.getHelixLength2(0, 3) == hhU.getHelixLength1(0, 3));


	}

	SECTION("Helix: Case 5 - unpaired bases in sequence 1 ", "[HelixHandlerUnpaired]") {
		// Case 4 -NO HELIX POSSIBLE
		RnaSequence r1("r1", "GGAAGG");
		RnaSequence r2("r2", "CCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 2, 999, 0, false);
		HelixHandlerUnpaired hhU(energy, hC);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   FILLHELIX  ////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		REQUIRE(hhU.fillHelix(0, energy.size1() - 1, 0, energy.size2() - 1) == 9);

		// (0,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(0, 0)) == -3);
		REQUIRE(hhU.getHelixLength1(0, 0) == 6);
		REQUIRE(hhU.getHelixLength2(0, 0) == 4);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		//////////////////////
		Interaction interaction(r1,r2);
		hhU.traceBackHelix(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 2);
		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 1);
		REQUIRE(interaction.basePairs.begin()->second == 2);

		REQUIRE(interaction.basePairs.rbegin()->first == 4);
		REQUIRE(interaction.basePairs.rbegin()->second == 1);

	}

	SECTION("Helix: Case 6 - 2 base pairs ", "[HelixHandlerUnpaired]") {
		// Case 4 -NO HELIX POSSIBLE
		RnaSequence r1("r1", "GAG");
		RnaSequence r2("r2", "CAC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 2, 999, 0, false);
		HelixHandlerUnpaired hhU(energy, hC);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   FILLHELIX  ////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		REQUIRE(hhU.fillHelix(0, energy.size1() - 1, 0, energy.size2() - 1) == 1);

		// (0,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(0, 0)) == -1);
		REQUIRE(hhU.getHelixLength1(0, 0) == 3);
		REQUIRE(hhU.getHelixLength2(0, 0) == 3);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		//////////////////////
		Interaction interaction(r1,r2);
		hhU.traceBackHelix(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 0);

	}

	SECTION("Helix: Case 7 - unpaired bases after 2,2 ", "[HelixHandlerUnpaired]") {

		RnaSequence r1("r1", "GGGAGG");
		RnaSequence r2("r2", "CCACCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 5, 2, 999, 0, false);
		HelixHandlerUnpaired hhU(energy, hC);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   FILLHELIX  ////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		REQUIRE(hhU.fillHelix(0, energy.size1() - 1, 0, energy.size2() - 1) == 16);

		// (0,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(0, 0)) == -4);
		REQUIRE(hhU.getHelixLength1(0, 0) == 6);
		REQUIRE(hhU.getHelixLength2(0, 0) == 6);

		// (0,1)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(0, 1)) == -3);
		REQUIRE(hhU.getHelixLength1(0, 1) == 5);
		REQUIRE(hhU.getHelixLength2(0, 1) == 5);

		// (0,2)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(0, 2)) == -2);
		REQUIRE(hhU.getHelixLength1(0, 2) == 3);
		REQUIRE(hhU.getHelixLength2(0, 2) == 4);

		// (0,4)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(0, 4)) == -1);
		REQUIRE(hhU.getHelixLength1(0, 4) == 2);
		REQUIRE(hhU.getHelixLength2(0, 4) == 2);

		// (1,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(1, 0)) == -3);
		REQUIRE(hhU.getHelixLength1(1, 0) == 5);
		REQUIRE(hhU.getHelixLength2(1, 0) == 5);

		// (1,1)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(1, 1)) == -3);
		REQUIRE(hhU.getHelixLength1(1, 1) == 5);
		REQUIRE(hhU.getHelixLength2(1, 1) == 5);

		// (1,2)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(1, 2)) == -2);
		REQUIRE(hhU.getHelixLength1(1, 2) == 4);
		REQUIRE(hhU.getHelixLength2(1, 2) == 4);

		// (1,4)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(1, 4)) == -1);
		REQUIRE(hhU.getHelixLength1(1, 4) == 2);
		REQUIRE(hhU.getHelixLength2(1, 4) == 2);

		// (2,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(2, 0)) == -2);
		REQUIRE(hhU.getHelixLength1(2, 0) == 4);
		REQUIRE(hhU.getHelixLength2(2, 0) == 3);

		// (2,1)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(2, 1)) == -2);
		REQUIRE(hhU.getHelixLength1(2, 1) == 4);
		REQUIRE(hhU.getHelixLength2(2, 1) == 4);

		// (2,2)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(2, 2)) == -2);
		REQUIRE(hhU.getHelixLength1(2, 2) == 4);
		REQUIRE(hhU.getHelixLength2(2, 2) == 4);

		// (2,4)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(2, 4)) == -1);
		REQUIRE(hhU.getHelixLength1(2, 4) == 3);
		REQUIRE(hhU.getHelixLength2(2, 4) == 2);

		// (4,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(4, 0)) == -1);
		REQUIRE(hhU.getHelixLength1(4, 0) == 2);
		REQUIRE(hhU.getHelixLength2(4, 0) == 2);

		// (4,1)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(4, 1)) == -1);
		REQUIRE(hhU.getHelixLength1(4, 1) == 2);
		REQUIRE(hhU.getHelixLength2(4, 1) == 2);

		// (4,2)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(4, 2)) == -1);
		REQUIRE(hhU.getHelixLength1(4, 2) == 2);
		REQUIRE(hhU.getHelixLength2(4, 2) == 3);

		// (4,4)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(4, 4)) == -1);
		REQUIRE(hhU.getHelixLength1(4, 4) == 2);
		REQUIRE(hhU.getHelixLength2(4, 4) == 2);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (2,2)
		//////////////////////
		Interaction interaction(r1,r2);
		hhU.traceBackHelix(interaction, 2, 2);

		REQUIRE(interaction.basePairs.size() == 1);
		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 4);
		REQUIRE(interaction.basePairs.begin()->second == 1);

		REQUIRE(interaction.basePairs.rbegin()->first == 4);
		REQUIRE(interaction.basePairs.rbegin()->second == 1);

	}

	SECTION("Helix: Case 8 - unpaired bases, many unpaired bases", "[HelixHandlerUnpaired]") {
		RnaSequence r1("r1", "GAAGAGG");
		RnaSequence r2("r2", "CAACACC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 2, 999, 0, false);
		HelixHandlerUnpaired hhU(energy, hC);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   FILLHELIX  ////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		REQUIRE(hhU.fillHelix(0, energy.size1() - 1, 0, energy.size2() - 1) == 6);

		// (0,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(0, 0)) == -3);
		REQUIRE(hhU.getHelixLength1(0, 0) == 7);
		REQUIRE(hhU.getHelixLength2(0, 0) == 7);

		// (3,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(3, 0)) == -2);
		REQUIRE(hhU.getHelixLength1(3, 0) == 4);
		REQUIRE(hhU.getHelixLength2(3, 0) == 4);

		// (3,1)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(3, 1)) == -2);
		REQUIRE(hhU.getHelixLength1(3, 1) == 4);
		REQUIRE(hhU.getHelixLength2(3, 1) == 6);

		// (5,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(5, 0)) == -1);
		REQUIRE(hhU.getHelixLength1(5, 0) == 2);
		REQUIRE(hhU.getHelixLength2(5, 0) == 2);

		// (5,1)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(5, 1)) == -1);
		REQUIRE(hhU.getHelixLength1(5, 1) == 2);
		REQUIRE(hhU.getHelixLength2(5, 1) == 3);

		// (5,3)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(5, 3)) == -1);
		REQUIRE(hhU.getHelixLength1(5, 3) == 2);
		REQUIRE(hhU.getHelixLength2(5, 3) == 4);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (2,2)
		//////////////////////
		Interaction interaction(r1,r2);
		hhU.traceBackHelix(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 2);
		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 3);
		REQUIRE(interaction.basePairs.begin()->second == 5);

		REQUIRE(interaction.basePairs.rbegin()->first == 5);
		REQUIRE(interaction.basePairs.rbegin()->second == 3);

	}

	SECTION("Helix: Case 9 - unpaired bases, many unpaired bases --only allow 1 unpaired base", "[HelixHandlerUnpaired]") {
		RnaSequence r1("r1", "GAAGAGG");
		RnaSequence r2("r2", "CAACACC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 1, 999, 0, false);
		HelixHandlerUnpaired hhU(energy, hC);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   FILLHELIX  ////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		REQUIRE(hhU.fillHelix(0, energy.size1() - 1, 0, energy.size2() - 1) == 3);

		// (3,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(3, 0)) == -2);
		REQUIRE(hhU.getHelixLength1(3, 0) == 4);
		REQUIRE(hhU.getHelixLength2(3, 0) == 4);

		// (5,0)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(5, 0)) == -1);
		REQUIRE(hhU.getHelixLength1(5, 0) == 2);
		REQUIRE(hhU.getHelixLength2(5, 0) == 2);

		// (5,1)
		REQUIRE(E_2_Ekcal(hhU.getHelixE(5, 1)) == -1);
		REQUIRE(hhU.getHelixLength1(5, 1) == 2);
		REQUIRE(hhU.getHelixLength2(5, 1) == 3);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (5,0)
		//////////////////////
		Interaction interaction(r1,r2);
		hhU.traceBackHelix(interaction, 5, 0);

		REQUIRE(interaction.basePairs.size() == 0);

	}
}
