

#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/InteractionRange.h"

using namespace IntaRNA;

TEST_CASE( "InteractionRange", "[InteractionRange]" ) {

	// setup easylogging++ stuff if not already done
	#include "testEasyLoggingSetup.icc"

	RnaSequence r("test","AACCGGUU");


	SECTION("empty construction") {

		InteractionRange range( r, r );

		REQUIRE( range.isSane() );
		REQUIRE( range.s1 == &r );
		REQUIRE( range.s2 == &r );
		REQUIRE( range.r1.isAscending() );
		REQUIRE_FALSE( range.r1.isDescending() );
		REQUIRE( range.r2.isDescending() );
		REQUIRE_FALSE( range.r2.isAscending() );
		REQUIRE( range.r1.from < r.size() );
		REQUIRE( (range.r1.to < r.size() || range.r1.to == RnaSequence::lastPos) );
		REQUIRE( range.r2.from < r.size() );
		REQUIRE( (range.r2.to < r.size() || range.r2.to == RnaSequence::lastPos) );

	}

	// TODO init with interaction

	// TODO check operators

}
