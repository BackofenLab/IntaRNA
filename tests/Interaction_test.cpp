

#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/Interaction.h"
#include "IntaRNA/InteractionRange.h"

using namespace IntaRNA;

TEST_CASE( "Interaction", "[Interaction]" ) {

	// setup easylogging++ stuff if not already done
	#include "testEasyLoggingSetup.icc"

	RnaSequence r("test","AACCGGUU");

	SECTION("empty construction") {

		Interaction inter( r, r );

		REQUIRE( inter.isEmpty() );
		REQUIRE( inter.s1 == &r );
		REQUIRE( inter.s2 == &r );
		REQUIRE_FALSE( inter.isValid() );
	}

	SECTION("add base pairs and sort") {

		Interaction inter( r, r );

		inter.basePairs.push_back( Interaction::BasePair( 1, 6 ) );
		REQUIRE_FALSE( inter.isEmpty() );
		REQUIRE( inter.isValid() );
		REQUIRE( Interaction::dotBracket(inter) == "(&)");
		REQUIRE( Interaction::dotBar(inter) == "2|&7|");

		inter.basePairs.push_back( Interaction::BasePair( 0, 7 ) );
		REQUIRE_FALSE( inter.isEmpty() );
		REQUIRE_FALSE( inter.isValid() );

		inter.sort();
		REQUIRE( inter.isValid() );
		REQUIRE( Interaction::dotBracket(inter) == "((&))");
		REQUIRE( Interaction::dotBar(inter) == "1||&7||");

		inter.basePairs.push_back( Interaction::BasePair( 3, 4 ) );
		REQUIRE( inter.isValid() );
		REQUIRE( Interaction::dotBracket(inter) == "((.(&).))");
		REQUIRE( Interaction::dotBar(inter) == "1||.|&5|.||");

	}


	SECTION("init with interaction range") {

		InteractionRange range(r,r, IndexRange(0,2), IndexRange(7,5));
		REQUIRE( range.isSane() );

		Interaction inter(range);
		REQUIRE_FALSE( inter.isEmpty() );
		REQUIRE( inter.isValid() );
		REQUIRE( inter.basePairs.size() == 2 );
		REQUIRE( inter.basePairs.at(0).first == 0 );
		REQUIRE( inter.basePairs.at(0).second == 7 );
		REQUIRE( inter.basePairs.at(1).first == 2 );
		REQUIRE( inter.basePairs.at(1).second == 5 );

	}

	SECTION("assign interaction range") {

		InteractionRange range(r,r, IndexRange(0,2), IndexRange(7,5));
		REQUIRE( range.isSane() );

		// empty interaction
		Interaction inter(r,r);
		REQUIRE( inter.isEmpty() );

		// assignment operator
		inter = range;
		REQUIRE_FALSE( inter.isEmpty() );
		REQUIRE( inter.isValid() );
		REQUIRE( inter.basePairs.size() == 2 );
		REQUIRE( inter.basePairs.at(0).first == 0 );
		REQUIRE( inter.basePairs.at(0).second == 7 );
		REQUIRE( inter.basePairs.at(1).first == 2 );
		REQUIRE( inter.basePairs.at(1).second == 5 );

	}

}
