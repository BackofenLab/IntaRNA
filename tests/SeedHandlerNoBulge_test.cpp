
#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/SeedHandlerNoBulge.h"
#include "IntaRNA/SeedHandlerMfe.h"

using namespace IntaRNA;

TEST_CASE( "SeedHandlerNoBulge", "[SeedHandlerNoBulge]" ) {

	// setup easylogging++ stuff if not already done
	#include "testEasyLoggingSetup.icc"

	// CAAAACAAACAAAACAAACCAAAACAAAA
	// CCUUUUUUCUUUUCUUUCUUUUUUCUUUU

	RnaSequence rna1("test1","CAAAACAAACAAAACAAACCAAAACAAAA");
	RnaSequence rna2("test2","UUUUCUUUUUUCUUUCUUUUCUUUUUUCC"); // had to be reversed first!

	AccessibilityDisabled acc1(rna1,rna1.size(),NULL);
	AccessibilityDisabled acc2(rna2,rna2.size(),NULL);
	ReverseAccessibility rAcc2(acc2); // reverse reverse
	InteractionEnergyBasePair energy( acc1, rAcc2 );

//	LOG(DEBUG) <<"idx:       01234567890123456789012345678";
//	LOG(DEBUG) <<"s1 : "<<acc1.getSequence();
//	LOG(DEBUG) <<"s2 : "<<acc2.getSequence();
//	LOG(DEBUG) <<"r2 : "<<rAcc2.getSequence();

	// seed constraint for 3 bp
	SeedConstraint sConstr(3,0,0,0,0,10,0,IndexRangeList(),IndexRangeList(),"", false, false, false);

	SECTION( "compare with SeedHandlerMfe output" ) {

		SeedHandlerMfe shMfe(energy,sConstr);
		SeedHandlerNoBulge sh(energy,sConstr);

		// same number of seeds
		REQUIRE( sh.fillSeed(0,energy.size1()-1,0,energy.size2()-1)
					== shMfe.fillSeed(0,energy.size1()-1,0,energy.size2()-1) );

		// compare seed energies
		for (size_t i1=0; i1<energy.size1(); i1++) {
			for (size_t i2=0; i2<energy.size2(); i2++) {
				if (E_isINF(sh.getSeedE(i1,i2))) {
					REQUIRE( E_isINF(shMfe.getSeedE(i1,i2)) );
					REQUIRE( ! shMfe.isSeedBound(i1,i2)) ;
				} else {
//					LOG(DEBUG) <<"seed at "<<i1<<"-"<<i2;
					REQUIRE( E_equal( sh.getSeedE(i1,i2), shMfe.getSeedE(i1,i2) ) );
					REQUIRE( shMfe.isSeedBound(i1,i2)) ;
				}
			}
		}

	}

	SECTION( "no GU bp" ) {

		RnaSequence rna1("test1","GGGGGUUUU");
		AccessibilityDisabled acc1(rna1,rna1.size(),NULL);
		ReverseAccessibility rAcc2(acc1); // reverse reverse
		InteractionEnergyBasePair energy( acc1, rAcc2 );

		{
			// should find some seed
			SeedHandlerNoBulge sh(energy,sConstr);
			REQUIRE( sh.fillSeed(0,energy.size1()-1,0,energy.size2()-1) > 0 );
		}
		{
			// should find no seed
			SeedConstraint sConstr(3,0,0,0,0,10,0,IndexRangeList(),IndexRangeList(),"", true, false, false);
			SeedHandlerNoBulge sh(energy,sConstr);
			REQUIRE( sh.fillSeed(0,energy.size1()-1,0,energy.size2()-1) == 0 );
		}

	}

	SECTION( "maxE" ) {

		RnaSequence rna1("test1","GGGGGUUUU");
		AccessibilityDisabled acc1(rna1,rna1.size(),NULL);
		ReverseAccessibility rAcc2(acc1); // reverse reverse
		InteractionEnergyBasePair energy( acc1, rAcc2 );

		{
			// should find some seeds
			SeedConstraint sConstr(3,0,0,0,Ekcal_2_E(-1),10,0,IndexRangeList(),IndexRangeList(),"", false, false, false);
			SeedHandlerNoBulge sh(energy,sConstr);
			REQUIRE( sh.fillSeed(0,energy.size1()-1,0,energy.size2()-1) > 0 );
		}
		{
			// should find some seed
			SeedConstraint sConstr(3,0,0,0,Ekcal_2_E(-2),10,0,IndexRangeList(),IndexRangeList(),"", false, false, false);
			SeedHandlerNoBulge sh(energy,sConstr);
			REQUIRE( sh.fillSeed(0,energy.size1()-1,0,energy.size2()-1) > 0 );
		}
		{
			// should find no seed
			SeedConstraint sConstr(3,0,0,0,Ekcal_2_E(-3),10,0,IndexRangeList(),IndexRangeList(),"", false, false, false);
			SeedHandlerNoBulge sh(energy,sConstr);
			REQUIRE( sh.fillSeed(0,energy.size1()-1,0,energy.size2()-1) == 0 );
		}

	}

	SECTION( "maxEhybrid" ) {

		RnaSequence rna1("test1","GGGGGUUUU");
		AccessibilityDisabled acc1(rna1,rna1.size(),NULL);
		ReverseAccessibility rAcc2(acc1); // reverse reverse
		InteractionEnergyBasePair energy( acc1, rAcc2 );

		{
			// should find some seeds
			SeedConstraint sConstr(3,0,0,0,0,10,Ekcal_2_E(-1),IndexRangeList(),IndexRangeList(),"", false, false, false);
			SeedHandlerNoBulge sh(energy,sConstr);
			REQUIRE( sh.fillSeed(0,energy.size1()-1,0,energy.size2()-1) > 0 );
		}
		{
			// should find some seeds
			SeedConstraint sConstr(3,0,0,0,0,10,Ekcal_2_E(-2),IndexRangeList(),IndexRangeList(),"", false, false, false);
			SeedHandlerNoBulge sh(energy,sConstr);
			REQUIRE( sh.fillSeed(0,energy.size1()-1,0,energy.size2()-1) > 0 );
		}
		{
			// should find no seed
			SeedConstraint sConstr(3,0,0,0,0,10,Ekcal_2_E(-3),IndexRangeList(),IndexRangeList(),"", false, false, false);
			SeedHandlerNoBulge sh(energy,sConstr);
			REQUIRE( sh.fillSeed(0,energy.size1()-1,0,energy.size2()-1) == 0 );
		}

	}

	SECTION( "areLoopOverlapping" ) {

		RnaSequence rna1("test1","GGGGGGGUUUUUUU");
		AccessibilityDisabled acc1(rna1,rna1.size(),NULL);
		ReverseAccessibility rAcc2(acc1); // reverse reverse
		InteractionEnergyBasePair energy( acc1, rAcc2 );

		// should find some seeds
		SeedConstraint sConstr(3,0,0,0,0,10,Ekcal_2_E(-1),IndexRangeList(),IndexRangeList(),"", false, false, false);
		SeedHandlerNoBulge sh(energy,sConstr);
		REQUIRE( sh.fillSeed(0,energy.size1()-1,0,energy.size2()-1) > 0 );

		REQUIRE( sh.isSeedBound(1,1) );
		REQUIRE( sh.isSeedBound(2,2) );
		REQUIRE( sh.isSeedBound(2,3) );
		REQUIRE( sh.isSeedBound(3,2) );
		REQUIRE( sh.isSeedBound(3,3) );
		REQUIRE( sh.isSeedBound(5,5) );

		REQUIRE( sh.areLoopOverlapping(1,1,1,1) );
		REQUIRE( sh.areLoopOverlapping(1,1,2,2) );
		REQUIRE( sh.areLoopOverlapping(2,2,1,1) );
		REQUIRE_FALSE( sh.areLoopOverlapping(1,1,3,3) );
		REQUIRE_FALSE( sh.areLoopOverlapping(3,3,1,1) );
		REQUIRE_FALSE( sh.areLoopOverlapping(1,1,2,3) );
		REQUIRE_FALSE( sh.areLoopOverlapping(2,3,1,1) );
		REQUIRE_FALSE( sh.areLoopOverlapping(1,1,3,2) );
		REQUIRE_FALSE( sh.areLoopOverlapping(3,2,1,1) );
		REQUIRE_FALSE( sh.areLoopOverlapping(1,1,5,5) );
		REQUIRE_FALSE( sh.areLoopOverlapping(5,5,1,1) );

	}

}
