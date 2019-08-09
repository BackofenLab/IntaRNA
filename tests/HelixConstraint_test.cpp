
#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/HelixConstraint.h"

using namespace IntaRNA;

TEST_CASE( "HelixConstraint", "[HelixConstraint]" ) {

	SECTION( "getter", "[HelixConstraint]" ) {

		size_t minBP= 2, maxBP = 10, maxIL=2, maxED=0, maxE=0;
		bool evaFullED=true;

		HelixConstraint hC( minBP, maxBP, maxIL, maxED, maxE, evaFullED);

		// check data access
		REQUIRE( hC.getMinBasePairs() == 2 );
		REQUIRE( hC.getMaxBasePairs() == 10 );
		REQUIRE( hC.getMaxIL() == 2 );
		REQUIRE( hC.getMaxED() == 0 );
		REQUIRE( hC.getMaxE() == 0 );
		REQUIRE( hC.evalFullE() == evaFullED );
		REQUIRE( hC.getMaxLength1() == 28 );
		REQUIRE( hC.getMaxLength2() == 28 );
	}


}
