
#include "catch.hpp"

#undef NDEBUG
#define protected public

#include "IntaRNA/RnaSequence.h"
#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/HelixConstraint.h"
#include "IntaRNA/HelixHandlerNoBulgeMax.h"
#include "IntaRNA/HelixHandlerIdxOffset.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/Interaction.h"
#include "IntaRNA/ReverseAccessibility.h"


using namespace IntaRNA;

TEST_CASE( "HelixHandlerIdxOffset for NoBulgeMax", "[HelixHandlerIdxOffset]") {

	SECTION("getter", "[HelixHandlerIdxOffset]") {

		RnaSequence r1("r1", "GGGGA");
		RnaSequence r2("r2", "ACCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 10, 0, 999, 0, false);
		HelixHandlerIdxOffset hhIO(new HelixHandlerNoBulgeMax(energy, hC));

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

	SECTION("Helix with Offset: Case1 - offset 1", "[HelixHandlerIdxOffset]") {

		RnaSequence r1("r1", "GGGAAGG");
		RnaSequence r2("r2", "CCAACCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2,4, 0, 999, 0, false);
		HelixHandlerIdxOffset hhIO(new HelixHandlerNoBulgeMax(energy, hC));

		// Set the offsets
		hhIO.setOffset1(1);
		hhIO.setOffset2(1);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////   FILLHELIXOFFSET  /////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// possible helices
		REQUIRE(hhIO.fillHelix(0, energy.size1()-hhIO.getOffset1()-1, 0, energy.size2()-hhIO.getOffset2()-1)==4);

		// All Working
		// (0,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(0, 0)) == -1);
		REQUIRE(hhIO.getHelixLength1(0, 0) == 2);
		REQUIRE(hhIO.getHelixLength2(0, 0) == hhIO.getHelixLength1(0, 0));

		// (0,4)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(0, 4)) == -1);
		REQUIRE(hhIO.getHelixLength1(0, 4) == 2);
		REQUIRE(hhIO.getHelixLength2(0, 4) == hhIO.getHelixLength1(0, 4));

		// (0,4)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(4, 0)) == -1);
		REQUIRE(hhIO.getHelixLength1(4, 0) == 2);
		REQUIRE(hhIO.getHelixLength2(4, 0) == hhIO.getHelixLength1(4, 0));

		// (0,4)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(4, 4)) == -1);
		REQUIRE(hhIO.getHelixLength1(4, 4) == 2);
		REQUIRE(hhIO.getHelixLength2(4, 4) == hhIO.getHelixLength1(4, 4));

		// Not Working
		// (3,4)
		REQUIRE(hhIO.getHelixE(3, 4) == E_INF);
		REQUIRE(hhIO.getHelixLength1(3, 4) == 0);
		REQUIRE(hhIO.getHelixLength2(3, 4) == hhIO.getHelixLength1(3, 4));

		// (1,1)
		REQUIRE(hhIO.getHelixE(1, 1) == E_INF);
		REQUIRE(hhIO.getHelixLength1(1, 1) == 0);
		REQUIRE(hhIO.getHelixLength2(1, 1) == hhIO.getHelixLength1(1, 1));

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		Interaction interaction(r1,r2);
		hhIO.traceBackHelix(interaction,0,0);

		REQUIRE(interaction.basePairs.size() == 0);

	}

	SECTION("Helix with Offset: Case 2 - all complementary", "[HelixHandlerIdxOffset]") {

		RnaSequence r1("r1", "GGGGGG");
		RnaSequence r2("r2", "CCCCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2,4, 0, 999, 0, false);
		HelixHandlerIdxOffset hhIO(new HelixHandlerNoBulgeMax(energy, hC));

		// Set the offsets
		hhIO.setOffset1(1);
		hhIO.setOffset2(1);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////   FILLHELIXOFFSET  /////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// possible helices
		REQUIRE(hhIO.fillHelix(0, energy.size1()-hhIO.getOffset1()-1, 0, energy.size2()-hhIO.getOffset2()-1)==16);

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
		Interaction interaction(r1,r2);

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

	SECTION("Helix with Offset: Case 3 - no interaction possible", "[HelixHandlerIdxOffset]") {

		RnaSequence r1("r1", "AAAAAAA");
		RnaSequence r2("r2", "AAAAAAA");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2,4, 0, 999, 0, false);
		HelixHandlerIdxOffset hhIO(new HelixHandlerNoBulgeMax(energy, hC));

		// Set the offsets
		hhIO.setOffset1(1);
		hhIO.setOffset2(1);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////   FILLHELIXOFFSET  /////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// possible helices
		REQUIRE(hhIO.fillHelix(0, energy.size1()-hhIO.getOffset1()-1, 0, energy.size2()-hhIO.getOffset2()-1)==0);

		// Not Working

		// (0,0)
		REQUIRE(hhIO.getHelixE(0, 0) == E_INF);
		REQUIRE(hhIO.getHelixLength1(0, 0) == 0);
		REQUIRE(hhIO.getHelixLength2(0, 0) == hhIO.getHelixLength1(0, 0));

		// (1,1)
		REQUIRE(hhIO.getHelixE(1, 1) == E_INF);
		REQUIRE(hhIO.getHelixLength1(1, 1) == 0);
		REQUIRE(hhIO.getHelixLength2(1, 1) == hhIO.getHelixLength1(1, 1));

		// (3,4)
		REQUIRE(hhIO.getHelixE(3, 4) == E_INF);
		REQUIRE(hhIO.getHelixLength1(3, 4) == 0);
		REQUIRE(hhIO.getHelixLength2(3, 4) == hhIO.getHelixLength1(3, 4));

		// (1,1)
		REQUIRE(hhIO.getHelixE(1, 1) == E_INF);
		REQUIRE(hhIO.getHelixLength1(1, 1) == 0);
		REQUIRE(hhIO.getHelixLength2(1, 1) == hhIO.getHelixLength1(1, 1));


	}

	SECTION("Helix with Offset: Case 3 - A 'wall' of A's disrupts the possible helices", "[HelixHandlerIdxOffset]") {
		// Case 2 - sequence containing an "A"-wall to disrupt perfect stacking
		RnaSequence r1("r1", "GGGGGAAGG");
		RnaSequence r2("r2", "CCAACCCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2,4, 0, 999, 0, false);
		HelixHandlerIdxOffset hhIO(new HelixHandlerNoBulgeMax(energy, hC));

		// Set the offsets
		hhIO.setOffset1(2);
		hhIO.setOffset2(2);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   FILLHELIX  ////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// possible helices
		REQUIRE(hhIO.fillHelix(0, energy.size1()-hhIO.getOffset1()-1, 0, energy.size2()-hhIO.getOffset2()-1)==9);

		REQUIRE_FALSE(energy.areComplementary(5,4));
		// All optimal combinations
		// (0,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(0, 0)) == -2);
		REQUIRE(hhIO.getHelixLength1(0, 0) == 3);
		REQUIRE(hhIO.getHelixLength2(0, 0) == hhIO.getHelixLength1(0, 0));

		// (0,1)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(0, 1)) == -1);
		REQUIRE(hhIO.getHelixLength1(0, 1) == 2);
		REQUIRE(hhIO.getHelixLength2(0, 1) == hhIO.getHelixLength1(0, 1));

		// (0,5)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(0, 5)) == -1);
		REQUIRE(hhIO.getHelixLength1(0, 5) == 2);
		REQUIRE(hhIO.getHelixLength2(0, 5) == hhIO.getHelixLength1(0, 5));

		// (1,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(1, 0)) == -1);
		REQUIRE(hhIO.getHelixLength1(1, 0) == 2);
		REQUIRE(hhIO.getHelixLength2(1, 0) == hhIO.getHelixLength1(1, 0));

		// (1,1)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(1, 1)) == -1);
		REQUIRE(hhIO.getHelixLength1(1, 1) == 2);
		REQUIRE(hhIO.getHelixLength2(1, 1) == hhIO.getHelixLength1(1, 1));

		// (1,5)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(1, 5)) == -1);
		REQUIRE(hhIO.getHelixLength1(1, 5) == 2);
		REQUIRE(hhIO.getHelixLength2(1, 5) == hhIO.getHelixLength1(1, 5));

		// (5,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(5, 0)) == -1);
		REQUIRE(hhIO.getHelixLength1(5, 0) == 2);
		REQUIRE(hhIO.getHelixLength2(5, 0) == hhIO.getHelixLength1(5, 0));

		// (5,1)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(5, 1)) == -1);
		REQUIRE(hhIO.getHelixLength1(5, 1) == 2);
		REQUIRE(hhIO.getHelixLength2(5, 1) == hhIO.getHelixLength1(5, 1));

		// (5,5)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(5, 5)) == -1);
		REQUIRE(hhIO.getHelixLength1(5, 5) == 2);
		REQUIRE(hhIO.getHelixLength2(5, 5) == hhIO.getHelixLength1(5, 5));



		// Not viable cases
		// (3,0) Not complementary
		REQUIRE(hhIO.getHelixE(3,0) == E_INF);
		REQUIRE(hhIO.getHelixLength1(3,0) == 0);
		REQUIRE(hhIO.getHelixLength2(3,0) == hhIO.getHelixLength1(3,0));

		// (5,4) no Helix possible (both sides too short)
		REQUIRE(hhIO.getHelixE(2,2) == E_INF);
		REQUIRE(hhIO.getHelixLength1(2,2) == 0);
		REQUIRE(hhIO.getHelixLength2(2,2) == hhIO.getHelixLength1(2,2));

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

		REQUIRE(interaction.basePairs.size() == 1);
		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 3);
		REQUIRE(interaction.basePairs.begin()->second == 5);

		REQUIRE(interaction.basePairs.rbegin()->first == 3);
		REQUIRE(interaction.basePairs.rbegin()->second == 5);

		// Case (5,5) - Possible but only 2 base pairs long (e.g no bp needs to be reported)
		//////////////////////

		interaction.clear();
		hhIO.traceBackHelix(interaction, 5, 5);

		REQUIRE(interaction.basePairs.size() == 0);

	}

	SECTION("Helix with Offset: Case 4 - uneven offset", "[HelixHandlerIdxOffset]") {
		// Case 2 - sequence containing an "A"-wall to disrupt perfect stacking
		RnaSequence r1("r1", "GGGGG");
		RnaSequence r2("r2", "CCCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2,4, 0, 999, 0, false);
		HelixHandlerIdxOffset hhIO(new HelixHandlerNoBulgeMax(energy, hC));

		// Set the offsets
		hhIO.setOffset1(1);
		hhIO.setOffset2(2);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   FILLHELIX  ////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// possible helices
		REQUIRE(hhIO.fillHelix(0, energy.size1()-hhIO.getOffset1()-1, 0, energy.size2()-hhIO.getOffset2()-1)==6);

		// All optimal combinations
		// (0,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(0, 0)) == -2);
		REQUIRE(hhIO.getHelixLength1(0, 0) == 3);
		REQUIRE(hhIO.getHelixLength2(0, 0) == hhIO.getHelixLength1(0, 0));

		// (0,1)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(0, 1)) == -1);
		REQUIRE(hhIO.getHelixLength1(0, 1) == 2);
		REQUIRE(hhIO.getHelixLength2(0, 1) == hhIO.getHelixLength1(0, 1));

		// (1,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(1, 0)) == -2);
		REQUIRE(hhIO.getHelixLength1(1, 0) == 3);
		REQUIRE(hhIO.getHelixLength2(1, 0) == hhIO.getHelixLength1(1, 0));

		// (1,1)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(1, 1)) == -1);
		REQUIRE(hhIO.getHelixLength1(1, 1) == 2);
		REQUIRE(hhIO.getHelixLength2(1, 1) == hhIO.getHelixLength1(1, 1));

		// (2,0)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(2, 0)) == -1);
		REQUIRE(hhIO.getHelixLength1(2, 0) == 2);
		REQUIRE(hhIO.getHelixLength2(2, 0) == hhIO.getHelixLength1(2, 0));

		// (5,1)
		REQUIRE(E_2_Ekcal(hhIO.getHelixE(2, 1)) == -1);
		REQUIRE(hhIO.getHelixLength1(2, 1) == 2);
		REQUIRE(hhIO.getHelixLength2(2, 1) == hhIO.getHelixLength1(2, 1));

		// Not viable cases
		// (3,0) Not complementary
		REQUIRE(hhIO.getHelixE(3,0) == E_INF);
		REQUIRE(hhIO.getHelixLength1(3,0) == 0);
		REQUIRE(hhIO.getHelixLength2(3,0) == hhIO.getHelixLength1(3,0));

		// (5,4) no Helix possible (both sides too short)
		REQUIRE(hhIO.getHelixE(2,2) == E_INF);
		REQUIRE(hhIO.getHelixLength1(2,2) == 0);
		REQUIRE(hhIO.getHelixLength2(2,2) == hhIO.getHelixLength1(2,2));

		// (5,3) no helix possible
		REQUIRE(hhIO.getHelixE(0,2) == E_INF);
		REQUIRE(hhIO.getHelixLength1(0,2) == 0);
		REQUIRE(hhIO.getHelixLength2(0,2) == hhIO.getHelixLength1(0,2));

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		//////////////////////
		Interaction interaction(r1,r2);
		hhIO.traceBackHelix(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 1);
		// First / last base pair of helix
		REQUIRE(interaction.basePairs.begin()->first == 2);
		REQUIRE(interaction.basePairs.begin()->second == 1);

		REQUIRE(interaction.basePairs.rbegin()->first == 2);
		REQUIRE(interaction.basePairs.rbegin()->second == 1);

		// Case (5,5) - Possible but only 2 base pairs long (e.g no bp needs to be reported)
		//////////////////////

		interaction.clear();
		hhIO.traceBackHelix(interaction, 2, 1);

		REQUIRE(interaction.basePairs.size() == 0);

	}
}
