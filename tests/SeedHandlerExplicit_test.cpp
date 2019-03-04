
#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/SeedHandlerExplicit.h"
#include "IntaRNA/AccessibilityDisabled.h"

using namespace IntaRNA;

TEST_CASE( "SeedHandlerExplicit", "[SeedHandlerExplicit]" ) {

	// setup easylogging++ stuff if not already done
	#include "testEasyLoggingSetup.icc"

	RnaSequence rna("test","ACGUACGU");

	AccessibilityDisabled acc(rna,rna.size(),NULL);
	ReverseAccessibility rAcc(acc);
	InteractionEnergyBasePair energy( acc, rAcc );

	SECTION( "test SeedData construction" ) {

		SeedHandlerExplicit::SeedData seedEmpty("1||&3||", energy);

		// valid input
		seedEmpty = SeedHandlerExplicit::SeedData ("1|&8|", energy);
		REQUIRE( seedEmpty.start1 == 0 );
		REQUIRE( seedEmpty.start2 == 0 );
		REQUIRE( seedEmpty.dotBar1 == "|" );
		REQUIRE( seedEmpty.dotBar2 == "|" );
		REQUIRE( E_equal( seedEmpty.energy, 0 ) );
		REQUIRE( seedEmpty.isValid() );

		// valid input
		seedEmpty = SeedHandlerExplicit::SeedData ("1||&7||", energy);
		REQUIRE( seedEmpty.start1 == 0 );
		REQUIRE( seedEmpty.start2 == 0 );
		REQUIRE( seedEmpty.dotBar1 == "||" );
		REQUIRE( seedEmpty.dotBar2 == "||" );
		REQUIRE( E_equal( seedEmpty.energy, energy.getE_basePair() ) );
		REQUIRE( seedEmpty.isValid() );

		// valid input
		seedEmpty = SeedHandlerExplicit::SeedData ("1|..|&1|..|", energy);
		REQUIRE( seedEmpty.start1 == 0 );
		REQUIRE( seedEmpty.start2 == 4 );
		REQUIRE( seedEmpty.dotBar1 == "|..|" );
		REQUIRE( seedEmpty.dotBar2 == "|..|" );
		REQUIRE( E_equal( seedEmpty.energy, energy.getE_basePair() ) );
		REQUIRE( seedEmpty.isValid() );

		// invalid out of bounds
		seedEmpty = SeedHandlerExplicit::SeedData("1|&10|", energy);
		REQUIRE( seedEmpty.start1 == 0 );
		REQUIRE( seedEmpty.start2 == std::string::npos );
		REQUIRE( seedEmpty.dotBar1 == "|" );
		REQUIRE( seedEmpty.dotBar2 == "|" );
		REQUIRE( E_isINF( seedEmpty.energy ) );
		REQUIRE_FALSE( seedEmpty.isValid() );

		// invalid out of bounds
		seedEmpty = SeedHandlerExplicit::SeedData("10|&1|", energy);
		REQUIRE( seedEmpty.start1 == std::string::npos );
		REQUIRE( seedEmpty.start2 == 7 );
		REQUIRE( seedEmpty.dotBar1 == "|" );
		REQUIRE( seedEmpty.dotBar2 == "|" );
		REQUIRE( E_isINF( seedEmpty.energy ) );
		REQUIRE_FALSE( seedEmpty.isValid() );

	}

	SECTION( "test encoding checks" ) {

		// valid encodings
		REQUIRE( SeedHandlerExplicit::checkSeedEncoding( "1||&1||" ) == "" );
		REQUIRE( SeedHandlerExplicit::checkSeedEncoding( "12||&1||" ) == "" );
		REQUIRE( SeedHandlerExplicit::checkSeedEncoding( "1||&51||" ) == "" );
		REQUIRE( SeedHandlerExplicit::checkSeedEncoding( "1|...|&51||" ) == "" );
		REQUIRE( SeedHandlerExplicit::checkSeedEncoding( "1||&51|...|" ) == "" );
		REQUIRE( SeedHandlerExplicit::checkSeedEncoding( "1|||&51|..|.|" ) == "" );
		REQUIRE( SeedHandlerExplicit::checkSeedEncoding( "1|.|.|&51|..|.|" ) == "" );
		REQUIRE( SeedHandlerExplicit::checkSeedEncoding( "1|.|.|&51|||" ) == "" );

		// missing stuff
		REQUIRE( SeedHandlerExplicit::checkSeedEncoding( "" ) != "" );
		REQUIRE( SeedHandlerExplicit::checkSeedEncoding( "||||" ) != "" );
		REQUIRE( SeedHandlerExplicit::checkSeedEncoding( "1||||" ) != "" );
		REQUIRE( SeedHandlerExplicit::checkSeedEncoding( "||2||" ) != "" );
		REQUIRE( SeedHandlerExplicit::checkSeedEncoding( "1||2||" ) != "" );
		REQUIRE( SeedHandlerExplicit::checkSeedEncoding( "||&||" ) != "" );
		REQUIRE( SeedHandlerExplicit::checkSeedEncoding( "1||&||" ) != "" );
		REQUIRE( SeedHandlerExplicit::checkSeedEncoding( "||&2||" ) != "" );

		// wrong content
		REQUIRE( SeedHandlerExplicit::checkSeedEncoding( "2|[|&2||" ) != "" );
		REQUIRE( SeedHandlerExplicit::checkSeedEncoding( "a||&2||" ) != "" );
		REQUIRE( SeedHandlerExplicit::checkSeedEncoding( "3||&2||2" ) != "" );
		REQUIRE( SeedHandlerExplicit::checkSeedEncoding( "3||&-2||" ) != "" );
		REQUIRE( SeedHandlerExplicit::checkSeedEncoding( "-3||&2||" ) != "" );

		// base pair ends missing
		REQUIRE( SeedHandlerExplicit::checkSeedEncoding( "2.||&2||" ) != "" );
		REQUIRE( SeedHandlerExplicit::checkSeedEncoding( "2||.&2||" ) != "" );
		REQUIRE( SeedHandlerExplicit::checkSeedEncoding( "2||&2.||" ) != "" );
		REQUIRE( SeedHandlerExplicit::checkSeedEncoding( "2||&2||." ) != "" );

	}

}
