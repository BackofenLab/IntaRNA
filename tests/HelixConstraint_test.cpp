
#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/HelixConstraint.h"

using namespace IntaRNA;

TEST_CASE( "HelixConstraint", "[HelixConstraint]" ) {

	SECTION( "getter", "[HelixConstraint]" ) {

		size_t minBP= 2, maxBP = 10, maxUP=2, maxIL=2, maxED=0, maxE=0;
		bool noED=false;

		HelixConstraint hC( minBP, maxBP, maxUP, maxIL, maxED, maxE, noED);

		// check data access
		REQUIRE( hC.getMinBasePairs() == 2 );
		REQUIRE( hC.getMaxBasePairs() == 10 );
		REQUIRE( hC.getMaxUnpaired() == 2 );
		REQUIRE( hC.getMaxIL() == 2 );
		REQUIRE( hC.getMaxED() == 0 );
		REQUIRE( hC.getMaxE() == 0 );
		REQUIRE( !hC.noED() );
		REQUIRE( hC.getMaxLength1() == 28 );
		REQUIRE( hC.getMaxLength2() == 28 );
	}


}