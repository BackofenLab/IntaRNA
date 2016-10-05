

#include "catch.hpp"

#undef NDEBUG

#include "AccessibilityConstraint.h"

TEST_CASE( "AccessibilityConstraint", "[AccessibilityConstraint]" ) {


	SECTION("check empty construction") {
		AccessibilityConstraint c(10);
		REQUIRE( c.isEmpty() );
		REQUIRE_FALSE( c.isBlocked(0) );
		REQUIRE_FALSE( c.isAccessible(0) );
	}

	SECTION("check dot-bracket construction") {
		AccessibilityConstraint c("..bb..xx..bb");
		REQUIRE_FALSE( c.isEmpty() );
		REQUIRE_FALSE( c.isBlocked(0) );
		REQUIRE_FALSE( c.isAccessible(0) );
		REQUIRE( c.isBlocked(3) );
		REQUIRE( c.isAccessible(6) );
		REQUIRE( c.isBlocked(10) );
	}


}

