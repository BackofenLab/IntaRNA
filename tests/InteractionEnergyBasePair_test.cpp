
#include "catch.hpp"

#undef NDEBUG

#include <cmath>
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/AccessibilityDisabled.h"

using namespace IntaRNA;

TEST_CASE( "InteractionEnergyBasePair", "[InteractionEnergyBasePair]" ) {

	// setup easylogging++ stuff if not already done
	#include "testEasyLoggingSetup.icc"

	RnaSequence rna("test","ACGU");

	AccessibilityDisabled acc(rna,rna.size(),NULL);
	ReverseAccessibility rAcc(acc);

	size_t maxLoop1 = 1, maxLoop2 = 2;
	InteractionEnergyBasePair energy( acc, rAcc, maxLoop1, maxLoop2, true, 1, Ekcal_2_E(-1.0), 1);

	SECTION("data access") {
		// check
		REQUIRE( &energy.getAccessibility1() == &acc );
		REQUIRE( &energy.getAccessibility2().getAccessibilityOrigin() == &acc );
		REQUIRE( energy.getMaxInternalLoopSize1() == maxLoop1 );
		REQUIRE( energy.getMaxInternalLoopSize2() == maxLoop2 );

		REQUIRE( energy.getAccessibility1().getSequence().asString().at(0) == 'A' );
		REQUIRE( energy.getAccessibility2().getSequence().asString().at(0) == 'U' );

	}

	SECTION("energy computation" ) {

		// dangling end check
		REQUIRE( energy.getE_danglingLeft(0,0) == 0.0 );
		REQUIRE( energy.getE_danglingLeft(2,2) == 0.0 );
		REQUIRE( energy.getE_danglingRight(0,0) == 0.0 );
		REQUIRE( energy.getE_danglingRight(2,2) == 0.0 );

		// dangling end check
		REQUIRE( energy.getE_danglingLeft(0,0) == 0.0 );
		REQUIRE( energy.getE_danglingLeft(2,2) == 0.0 );
		REQUIRE( energy.getE_danglingRight(0,0) == 0.0 );
		REQUIRE( energy.getE_danglingRight(2,2) == 0.0 );

		// init only, but interior loop called
		REQUIRE( energy.getE_interLeft( 0,0, 0,0 ) > 0.0 );

		// base pairs possible
		REQUIRE( energy.getE_interLeft( 0,1, 0,1 ) < 0.0 );

		// base pairs overlapping
		REQUIRE_FALSE( energy.getE_interLeft( 0,0, 0,1 ) < 0.0 );
		REQUIRE_FALSE( energy.getE_interLeft( 0,1, 0,0 ) < 0.0 );

		// base pairs not possible
		REQUIRE_FALSE( energy.getE_interLeft( 0,0, 1,1 ) < 0.0 );
		REQUIRE_FALSE( energy.getE_interLeft( 0,1, 2,3 ) < 0.0 );
		REQUIRE_FALSE( energy.getE_interLeft( 0,1, 0,2 ) < 0.0 );
		REQUIRE_FALSE( energy.getE_interLeft( 0,2, 1,2 ) < 0.0 );

		// loop size exceeded
		REQUIRE_FALSE( energy.getE_interLeft( 0,3, 1,2 ) < 0.0 );
		REQUIRE_FALSE( energy.getE_interLeft( 1,2, 0,3 ) < 0.0 );
		REQUIRE_FALSE( energy.getE_interLeft( 0,3, 0,3 ) < 0.0 );

	}

  SECTION("ES computation") {
    REQUIRE( E_equal(energy.getES1(0, 3), Ekcal_2_E(-1.313186)) );
    REQUIRE( E_equal(energy.getES2(0, 3), Ekcal_2_E(-1.313186)) );
    REQUIRE( E_isINF(energy.getES1(0, 2)) );
    REQUIRE( E_isINF(energy.getES1(1, 2)) );
    REQUIRE( E_isINF(energy.getES2(0, 2)) );
    REQUIRE( E_isINF(energy.getES2(1, 2)) );
  }

}
