
#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/RnaSequence.h"
#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/HelixConstraint.h"
#include "IntaRNA/HelixHandlerUnpaired.h"
#include "IntaRNA/HelixHandlerIdxOffset.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/Interaction.h"
#include "IntaRNA/ReverseAccessibility.h"


using namespace IntaRNA;

TEST_CASE( "HelixHandlerIdxOffset for Unpaired", "[HelixHandlerIdxOffset]") {


	SECTION("getter", "[HelixHandlerIdxOffset]") {

		RnaSequence r1("r1", "GGGGG");
		RnaSequence r2("r2", "CCCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 10, 2, 999, 0, false);
		HelixHandlerIdxOffset hhIO(new HelixHandlerUnpaired(energy, hC));

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

	SECTION("Helix with Offset: Case 1 - Everything is complementary", "[HelixHandlerIdxOffset]") {

		// Case 1 Perfect sequence
		RnaSequence r1("r1", "AGGGGG");
		RnaSequence r2("r2", "CCCCCA");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 2, 999, 0, false);
		HelixHandlerIdxOffset hhIO(new HelixHandlerUnpaired(energy, hC));

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   FILLHELIX  ////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// When counting all non-inf values for helices -> 29
		// When only counting the optimal helices that are non-inf -> 16
		// (0,0) -> 4 ; (1,0) -> 4 ; (2,0) -> 3 ; (3,0) -> 2
		// (0,1) -> 4 ; (1,1) -> 4 ; (2,1) -> 3 ; (3,1) -> 2
		// (0,2) -> 3 ; (1,2) -> 3 ; (2,2) -> 3 ; (3,2) -> 2
		// (0,3) -> 2 ; (1,3) -> 2 ; (2,3) -> 2 ; (3,3) -> 2
		REQUIRE(hhIO.fillHelix(0, energy.size1() - 1, 0, energy.size2() - 1) == 16);

		// Set the offsets
		hhIO.setOffset1(1);
		hhIO.setOffset2(1);

		// All optimal combinations
		// (0,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(0, 0)) == -3);
		REQUIRE(hhIO.getHelixLength1(0, 0) == 4);
		REQUIRE(hhIO.getHelixLength2(0, 0) == hhIO.getHelixLength1(0, 0));

		// (0,1)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(0, 1)) == -3);
		REQUIRE(hhIO.getHelixLength1(0, 1) == 4);
		REQUIRE(hhIO.getHelixLength2(0, 1) == hhIO.getHelixLength1(0, 1));

		// (0,2)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(0, 2)) == -2);
		REQUIRE(hhIO.getHelixLength1(0, 2) == 3);
		REQUIRE(hhIO.getHelixLength2(0, 2) == hhIO.getHelixLength1(0, 2));

		// (0,3)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(0, 3)) == -1);
		REQUIRE(hhIO.getHelixLength1(0, 3) == 2);
		REQUIRE(hhIO.getHelixLength2(0, 3) == hhIO.getHelixLength1(0, 3));

		// (1,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(1, 0)) == -3);
		REQUIRE(hhIO.getHelixLength1(1, 0) == 4);
		REQUIRE(hhIO.getHelixLength2(1, 0) == hhIO.getHelixLength1(1, 0));

		// (1,1)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(1, 1)) == -3);
		REQUIRE(hhIO.getHelixLength1(1, 1) == 4);
		REQUIRE(hhIO.getHelixLength2(1, 1) == hhIO.getHelixLength1(1, 1));

		// (1,2)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(1, 2)) == -2);
		REQUIRE(hhIO.getHelixLength1(1, 2) == 3);
		REQUIRE(hhIO.getHelixLength2(1, 2) == hhIO.getHelixLength1(1, 2));

		// (1,3)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(1, 3)) == -1);
		REQUIRE(hhIO.getHelixLength1(1, 3) == 2);
		REQUIRE(hhIO.getHelixLength2(1, 3) == hhIO.getHelixLength1(1, 3));

		// (2,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(2, 0)) == -2);
		REQUIRE(hhIO.getHelixLength1(2, 0) == 3);
		REQUIRE(hhIO.getHelixLength2(2, 0) == hhIO.getHelixLength1(2, 0));

		// (2,1)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(2, 1)) == -2);
		REQUIRE(hhIO.getHelixLength1(2, 1) == 3);
		REQUIRE(hhIO.getHelixLength2(2, 1) == hhIO.getHelixLength1(2, 1));

		// (2,2)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(2, 2)) == -2);
		REQUIRE(hhIO.getHelixLength1(2, 2) == 3);
		REQUIRE(hhIO.getHelixLength2(2, 2) == hhIO.getHelixLength1(2, 2));

		// (2,3)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(2, 3)) == -1);
		REQUIRE(hhIO.getHelixLength1(2, 3) == 2);
		REQUIRE(hhIO.getHelixLength2(2, 3) == hhIO.getHelixLength1(2, 3));

		// (3,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(3, 0)) == -1);
		REQUIRE(hhIO.getHelixLength1(3, 0) == 2);
		REQUIRE(hhIO.getHelixLength2(3, 0) == hhIO.getHelixLength1(3, 0));

		// (3,1)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(3, 1)) == -1);
		REQUIRE(hhIO.getHelixLength1(3, 1) == 2);
		REQUIRE(hhIO.getHelixLength2(3, 1) == hhIO.getHelixLength1(3, 1));

		// (3,2)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(3, 2)) == -1);
		REQUIRE(hhIO.getHelixLength1(3, 2) == 2);
		REQUIRE(hhIO.getHelixLength2(3, 2) == hhIO.getHelixLength1(3, 2));

		// (3,3)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(3, 3)) == -1);
		REQUIRE(hhIO.getHelixLength1(3, 3) == 2);
		REQUIRE(hhIO.getHelixLength2(3, 3) == hhIO.getHelixLength1(3, 3));


		// Several cases that do not allow a helix

		// (4,4)
		REQUIRE(hhIO.getHelixE(4, 4) == E_INF);
		REQUIRE(hhIO.getHelixLength1(4, 4) == 0);
		REQUIRE(hhIO.getHelixLength2(4, 4) == hhIO.getHelixLength1(4, 4));

		// (4,3)
		REQUIRE(hhIO.getHelixE(0, 4) == E_INF);
		REQUIRE(hhIO.getHelixLength1(0, 4) == 0);
		REQUIRE(hhIO.getHelixLength2(0, 4) == hhIO.getHelixLength1(0, 4));


		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		//////////////////////
		Interaction interaction(r1, r2);

		hhIO.traceBackHelix(interaction, 0, 0);


		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 2);
		REQUIRE(interaction.basePairs.begin()->second == 3);

		REQUIRE(interaction.basePairs.rbegin()->first == 3);
		REQUIRE(interaction.basePairs.rbegin()->second == 2);



		// Case (2,1)
		//////////////////////
		interaction.clear();

		hhIO.traceBackHelix(interaction, 2, 1);

		REQUIRE(interaction.basePairs.size() == 1);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 4);
		REQUIRE(interaction.basePairs.begin()->second == 2);

		REQUIRE(interaction.basePairs.rbegin()->first == 4);
		REQUIRE(interaction.basePairs.rbegin()->second == 2);
	}

	SECTION("Helix with offset: Case 2 - Sequence 1 contains an A", "[HelixHandlerIdxOffset]") {
		// Case 2 - sequence containing an "A" to disrupt perfect stacking
		RnaSequence r1("r1", "AGGGAGG");
		RnaSequence r2("r2", "CCCCCA");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 2, 999, 0, false);
		HelixHandlerIdxOffset hhIO(new HelixHandlerUnpaired(energy, hC));

		// set offsets
		hhIO.setOffset1(1);
		hhIO.setOffset2(1);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   FILLHELIX  ////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		REQUIRE(hhIO.fillHelix(0, energy.size1()-hhIO.getOffset1() - 1, 0, energy.size2()-hhIO.getOffset2() - 1) == 16);


		// All optimal combinations
		// (0,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(0, 0)) == -3);
		REQUIRE(hhIO.getHelixLength1(0, 0) == 5);
		REQUIRE(hhIO.getHelixLength2(0, 0) == 4);

		// (0,1)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(0, 1)) == -3);
		REQUIRE(hhIO.getHelixLength1(0, 1) == 5);
		REQUIRE(hhIO.getHelixLength2(0, 1) == 4);

		// (0,2)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(0, 2)) == -2);
		REQUIRE(hhIO.getHelixLength1(0, 2) == 3);
		REQUIRE(hhIO.getHelixLength2(0, 2) == hhIO.getHelixLength1(0, 2));

		// (0,3)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(0, 3)) == -1);
		REQUIRE(hhIO.getHelixLength1(0, 3) == 2);
		REQUIRE(hhIO.getHelixLength2(0, 3) == hhIO.getHelixLength1(0, 3));

		// (1,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(1, 0)) == -3);
		REQUIRE(hhIO.getHelixLength1(1, 0) == 5);
		REQUIRE(hhIO.getHelixLength2(1, 0) == 4);

		// (1,1)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(1, 1)) == -3);
		REQUIRE(hhIO.getHelixLength1(1, 1) == 5);
		REQUIRE(hhIO.getHelixLength2(1, 1) == 4);

		// (1,2)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(1, 2)) == -2);
		REQUIRE(hhIO.getHelixLength1(1, 2) == 4);
		REQUIRE(hhIO.getHelixLength2(1, 2) == 3);

		// (1,3)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(1, 3)) == -1);
		REQUIRE(hhIO.getHelixLength1(1, 3) == 2);
		REQUIRE(hhIO.getHelixLength2(1, 3) == hhIO.getHelixLength1(1, 3));

		// (2,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(2, 0)) == -2);
		REQUIRE(hhIO.getHelixLength1(2, 0) == 4);
		REQUIRE(hhIO.getHelixLength2(2, 0) == 3);

		// (2,1)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(2, 1)) == -2);
		REQUIRE(hhIO.getHelixLength1(2, 1) == 4);
		REQUIRE(hhIO.getHelixLength2(2, 1) == 3);

		// (2,2)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(2, 2)) == -2);
		REQUIRE(hhIO.getHelixLength1(2, 2) == 4);
		REQUIRE(hhIO.getHelixLength2(2, 2) == 3);

		// (2,3)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(2, 3)) == -1);
		REQUIRE(hhIO.getHelixLength1(2, 3) == 3);
		REQUIRE(hhIO.getHelixLength2(2, 3) == 2);

		// (4,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(4, 0)) == -1);
		REQUIRE(hhIO.getHelixLength1(4, 0) == 2);
		REQUIRE(hhIO.getHelixLength2(4, 0) == hhIO.getHelixLength1(4, 0));

		// (4,1)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(4, 1)) == -1);
		REQUIRE(hhIO.getHelixLength1(4, 1) == 2);
		REQUIRE(hhIO.getHelixLength2(4, 1) == hhIO.getHelixLength1(4, 1));

		// (4,2)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(4, 2)) == -1);
		REQUIRE(hhIO.getHelixLength1(4, 2) == 2);
		REQUIRE(hhIO.getHelixLength2(4, 2) == hhIO.getHelixLength1(4, 2));

		// (4,3)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(4, 3)) == -1);
		REQUIRE(hhIO.getHelixLength1(4, 3) == 2);
		REQUIRE(hhIO.getHelixLength2(4, 3) == hhIO.getHelixLength1(4, 3));


		// Not viable cases
		// (3,0) Not complementary
		REQUIRE(hhIO.getHelixE(3,0) == E_INF);
		REQUIRE(hhIO.getHelixLength1(3,0) == 0);
		REQUIRE(hhIO.getHelixLength2(3,0) == hhIO.getHelixLength1(3,0));

		// (5,4) no Helix possible (both sides too short)
		REQUIRE(hhIO.getHelixE(5,4) == E_INF);
		REQUIRE(hhIO.getHelixLength1(5,4) == 0);
		REQUIRE(hhIO.getHelixLength2(5,4) == hhIO.getHelixLength1(5,4));

		// (5,3) no helix possible
		REQUIRE(hhIO.getHelixE(5,3) == E_INF);
		REQUIRE(hhIO.getHelixLength1(5,3) == 0);
		REQUIRE(hhIO.getHelixLength2(5,3) == hhIO.getHelixLength1(5,3));

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		//////////////////////
		Interaction interaction(r1,r2);

		hhIO.traceBackHelix(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 2);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 2);
		REQUIRE(interaction.basePairs.begin()->second == 3);

		REQUIRE(interaction.basePairs.rbegin()->first == 3);
		REQUIRE(interaction.basePairs.rbegin()->second == 2);

		// Case (2,1)
		//////////////////////
		interaction.clear();

		hhIO.traceBackHelix(interaction, 2, 1);

		REQUIRE(interaction.basePairs.size() == 1);

		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 5);
		REQUIRE(interaction.basePairs.begin()->second == 2);

		REQUIRE(interaction.basePairs.rbegin()->first == 5);
		REQUIRE(interaction.basePairs.rbegin()->second == 2);

	}

	SECTION("Helix with Offset: Case 3 - unpaired bases in sequence 1 ", "[HelixHandlerUnpaired]") {
		// Case 4 -NO HELIX POSSIBLE
		RnaSequence r1("r1", "AGGAAGG");
		RnaSequence r2("r2", "CCCCAA");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 2, 999, 0, false);
		HelixHandlerIdxOffset hhIO(new HelixHandlerUnpaired(energy, hC));

		// set offsets
		hhIO.setOffset1(1);
		hhIO.setOffset2(2);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   FILLHELIX  ////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		REQUIRE(hhIO.fillHelix(0, energy.size1()-hhIO.getOffset1() - 1, 0, energy.size2()-hhIO.getOffset2() - 1) == 9);


		// (0,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(0, 0)) == -3);
		REQUIRE(hhIO.getHelixLength1(0, 0) == 6);
		REQUIRE(hhIO.getHelixLength2(0, 0) == 4);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		//////////////////////
		Interaction interaction(r1,r2);
		hhIO.traceBackHelix(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 2);
		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 2);
		REQUIRE(interaction.basePairs.begin()->second == 2);

		REQUIRE(interaction.basePairs.rbegin()->first == 5);
		REQUIRE(interaction.basePairs.rbegin()->second == 1);

	}

	SECTION("Helix: Case 8 - unpaired bases, many unpaired bases", "[HelixHandlerUnpaired]") {
		RnaSequence r1("r1", "AGAAGAGG");
		RnaSequence r2("r2", "CAACACCAA");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 2, 999, 0, false);
		HelixHandlerIdxOffset hhIO(new HelixHandlerUnpaired(energy, hC));

		// set offsets
		hhIO.setOffset1(1);
		hhIO.setOffset2(2);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   FILLHELIX  ////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		REQUIRE(hhIO.fillHelix(0, energy.size1()-hhIO.getOffset1() - 1, 0, energy.size2()-hhIO.getOffset2() - 1) == 6);

		// (0,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(0, 0)) == -3);
		REQUIRE(hhIO.getHelixLength1(0, 0) == 7);
		REQUIRE(hhIO.getHelixLength2(0, 0) == 7);

		// (3,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(3, 0)) == -2);
		REQUIRE(hhIO.getHelixLength1(3, 0) == 4);
		REQUIRE(hhIO.getHelixLength2(3, 0) == 4);

		// (3,1)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(3, 1)) == -2);
		REQUIRE(hhIO.getHelixLength1(3, 1) == 4);
		REQUIRE(hhIO.getHelixLength2(3, 1) == 6);

		// (5,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(5, 0)) == -1);
		REQUIRE(hhIO.getHelixLength1(5, 0) == 2);
		REQUIRE(hhIO.getHelixLength2(5, 0) == 2);

		// (5,1)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(5, 1)) == -1);
		REQUIRE(hhIO.getHelixLength1(5, 1) == 2);
		REQUIRE(hhIO.getHelixLength2(5, 1) == 3);

		// (5,3)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(5, 3)) == -1);
		REQUIRE(hhIO.getHelixLength1(5, 3) == 2);
		REQUIRE(hhIO.getHelixLength2(5, 3) == 4);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (2,2)
		//////////////////////
		Interaction interaction(r1,r2);
		hhIO.traceBackHelix(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 2);
		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 4);
		REQUIRE(interaction.basePairs.begin()->second == 5);

		REQUIRE(interaction.basePairs.rbegin()->first == 6);
		REQUIRE(interaction.basePairs.rbegin()->second == 3);

	}

	SECTION("Helix: Case 9 - unpaired bases, many unpaired bases --only allow 1 unpaired base", "[HelixHandlerUnpaired]") {
		RnaSequence r1("r1", "AAGAAGAGG");
		RnaSequence r2("r2", "CAACACCAA");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 1, 999, 0, false);
		HelixHandlerIdxOffset hhIO(new HelixHandlerUnpaired(energy, hC));

		// set offsets
		hhIO.setOffset1(2);
		hhIO.setOffset2(2);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   FILLHELIX  ////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		REQUIRE(hhIO.fillHelix(0, energy.size1()-hhIO.getOffset1() - 1, 0, energy.size2()-hhIO.getOffset2() - 1) == 3);

		// (3,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(3, 0)) == -2);
		REQUIRE(hhIO.getHelixLength1(3, 0) == 4);
		REQUIRE(hhIO.getHelixLength2(3, 0) == 4);

		// (5,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(5, 0)) == -1);
		REQUIRE(hhIO.getHelixLength1(5, 0) == 2);
		REQUIRE(hhIO.getHelixLength2(5, 0) == 2);

		// (5,1)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(5, 1)) == -1);
		REQUIRE(hhIO.getHelixLength1(5, 1) == 2);
		REQUIRE(hhIO.getHelixLength2(5, 1) == 3);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (5,0)
		//////////////////////
		Interaction interaction(r1,r2);
		hhIO.traceBackHelix(interaction, 5, 0);

		REQUIRE(interaction.basePairs.size() == 0);

	}
}
