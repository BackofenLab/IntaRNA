
#include "catch.hpp"

#undef NDEBUG

#include "InteractionEnergyBasePair.h"
#include "AccessibilityDisabled.h"

TEST_CASE( "InteractionEnergyBasePair", "[InteractionEnergyBasePair]" ) {

	RnaSequence rna("test","ACGU");

	AccessibilityDisabled acc(rna);
	ReverseAccessibility rAcc(acc);

	size_t maxLoop1 = 1, maxLoop2 = 2;
	InteractionEnergyBasePair energy( acc, rAcc, maxLoop1, maxLoop2);

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
		REQUIRE( energy.getDanglingLeft(0,0) == 0.0 );
		REQUIRE( energy.getDanglingLeft(2,2) == 0.0 );
		REQUIRE( energy.getDanglingRight(0,0) == 0.0 );
		REQUIRE( energy.getDanglingRight(2,2) == 0.0 );

		// dangling end check
		REQUIRE( energy.getDanglingLeft(0,0) == 0.0 );
		REQUIRE( energy.getDanglingLeft(2,2) == 0.0 );
		REQUIRE( energy.getDanglingRight(0,0) == 0.0 );
		REQUIRE( energy.getDanglingRight(2,2) == 0.0 );

		// base pairs possible
		REQUIRE( energy.getInterLoopE( 0,0, 0,0 ) < 0.0 );
		REQUIRE( energy.getInterLoopE( 0,1, 0,1 ) < 0.0 );

		// base pairs overlapping
		REQUIRE_FALSE( energy.getInterLoopE( 0,0, 0,1 ) < 0.0 );
		REQUIRE_FALSE( energy.getInterLoopE( 0,1, 0,0 ) < 0.0 );

		// base pairs not possible
		REQUIRE_FALSE( energy.getInterLoopE( 0,0, 1,1 ) < 0.0 );
		REQUIRE_FALSE( energy.getInterLoopE( 0,1, 2,3 ) < 0.0 );
		REQUIRE_FALSE( energy.getInterLoopE( 0,1, 0,2 ) < 0.0 );
		REQUIRE_FALSE( energy.getInterLoopE( 0,2, 1,2 ) < 0.0 );

		// loop size exceeded
		REQUIRE_FALSE( energy.getInterLoopE( 0,3, 1,2 ) < 0.0 );
		REQUIRE_FALSE( energy.getInterLoopE( 1,2, 0,3 ) < 0.0 );
		REQUIRE_FALSE( energy.getInterLoopE( 0,3, 0,3 ) < 0.0 );

	}

}
