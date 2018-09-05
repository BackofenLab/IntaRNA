
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
	SeedConstraint sConstr(3,0,0,0,0,10,IndexRangeList(),IndexRangeList(),"");


	SECTION( "compare with SeedHandlerMfe output" ) {

		SeedHandlerMfe shNoBulge(energy,sConstr);
		SeedHandlerNoBulge shMfe(energy,sConstr);

		// same number of seeds
		REQUIRE( shMfe.fillSeed(0,energy.size1()-1,0,energy.size2()-1)
					== shNoBulge.fillSeed(0,energy.size1()-1,0,energy.size2()-1) );

		// compare seed energies
		for (size_t i1=0; i1<energy.size1(); i1++) {
			for (size_t i2=0; i2<energy.size2(); i2++) {
				if (E_isINF(shMfe.getSeedE(i1,i2))) {
					REQUIRE( E_isINF(shNoBulge.getSeedE(i1,i2)) );
				} else {
//					LOG(DEBUG) <<"seed at "<<i1<<"-"<<i2;
					REQUIRE( E_equal( shMfe.getSeedE(i1,i2), shNoBulge.getSeedE(i1,i2) ) );
				}
			}
		}

	}

}
