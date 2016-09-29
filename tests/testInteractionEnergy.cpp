
#include "catch.hpp"

#undef NDEBUG

#include "InteractionEnergy.h"

TEST_CASE( "InteractionEnergy - isAllowedLoopRegion()", "[InteractionEnergy]" ) {

	RnaSequence rna("test","ACGUNA");

	// check (j-i) <= (w+1)
	REQUIRE( InteractionEnergy::isAllowedLoopRegion(rna, 0, 0, 0) );
	REQUIRE( InteractionEnergy::isAllowedLoopRegion(rna, 0, 1, 0) );
	REQUIRE_FALSE( InteractionEnergy::isAllowedLoopRegion(rna, 0, 2, 0) );
	REQUIRE( InteractionEnergy::isAllowedLoopRegion(rna, 1, 1, 0) );
	REQUIRE( InteractionEnergy::isAllowedLoopRegion(rna, 1, 2, 0) );
	REQUIRE_FALSE( InteractionEnergy::isAllowedLoopRegion(rna, 1, 3, 0) );

	// >= size()
	REQUIRE_FALSE( InteractionEnergy::isAllowedLoopRegion(rna, rna.size(), 0, rna.size()) );
	REQUIRE_FALSE( InteractionEnergy::isAllowedLoopRegion(rna, 0, rna.size(), rna.size()) );
	REQUIRE_FALSE( InteractionEnergy::isAllowedLoopRegion(rna, rna.size(), rna.size(), rna.size()) );
	// == N
	REQUIRE_FALSE( InteractionEnergy::isAllowedLoopRegion(rna, 0, 4, rna.size()) );
	REQUIRE_FALSE( InteractionEnergy::isAllowedLoopRegion(rna, 4, 5, rna.size()) );
	REQUIRE_FALSE( InteractionEnergy::isAllowedLoopRegion(rna, 4, 4, rna.size()) );
	// i > j
	REQUIRE_FALSE( InteractionEnergy::isAllowedLoopRegion(rna, 1, 0, 3) );

}
