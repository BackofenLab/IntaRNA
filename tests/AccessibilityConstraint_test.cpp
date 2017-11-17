

#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/AccessibilityConstraint.h"

using namespace IntaRNA;

TEST_CASE( "AccessibilityConstraint", "[AccessibilityConstraint]" ) {

	// setup easylogging++ stuff if not already done
	#include "testEasyLoggingSetup.icc"


	SECTION("check empty construction") {
		AccessibilityConstraint c(10);
		REQUIRE( c.isEmpty() );
		REQUIRE_FALSE( c.isMarkedBlocked(0) );
		REQUIRE_FALSE( c.isMarkedAccessible(0) );
		REQUIRE( c.isAccessible(0) );
		REQUIRE( c.isUnconstrained(0) );
	}

	SECTION("check dot-bracket construction") {
		std::string constraint = "..bb..xxp.bb";
		AccessibilityConstraint c(constraint);

		REQUIRE_FALSE( c.isEmpty() );
		REQUIRE_FALSE( c.isMarkedBlocked(0) );
		REQUIRE_FALSE( c.isMarkedAccessible(0) );
		REQUIRE( c.isAccessible(0) );
		REQUIRE( c.isUnconstrained(0) );

		REQUIRE( c.isMarkedBlocked(3) );
		REQUIRE_FALSE( c.isAccessible(3) );
		REQUIRE_FALSE( c.isUnconstrained(3) );
		REQUIRE_FALSE( c.isMarkedAccessible(3) );
		REQUIRE_FALSE( c.isMarkedPaired(3) );

		REQUIRE( c.isMarkedAccessible(6) );
		REQUIRE( c.isAccessible(6) );
		REQUIRE_FALSE( c.isUnconstrained(6) );

		REQUIRE( c.isMarkedPaired(8) );
		REQUIRE_FALSE( c.isUnconstrained(8) );
		REQUIRE_FALSE( c.isAccessible(8) );
		REQUIRE_FALSE( c.isMarkedBlocked(8) );
		REQUIRE_FALSE( c.isMarkedAccessible(8) );

		REQUIRE( c.isMarkedBlocked(10) );
		REQUIRE_FALSE( c.isAccessible(10) );

		std::string vrnaStyle = "";
		for (size_t i=0; i<constraint.size(); i++) {
			vrnaStyle += toString(c.getVrnaDotBracket(i));
		}
		REQUIRE( vrnaStyle == "..xx..xx|.xx");
	}


}

