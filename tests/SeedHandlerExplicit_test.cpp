
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
		REQUIRE( SeedHandlerExplicit::checkSeedEncoding( "3||&-2||" ) == "" );
		REQUIRE( SeedHandlerExplicit::checkSeedEncoding( "-3||&2||" ) == "" );

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
		REQUIRE( SeedHandlerExplicit::checkSeedEncoding( "3||&--2||" ) != "" );
		REQUIRE( SeedHandlerExplicit::checkSeedEncoding( "--3||&2||" ) != "" );

		// base pair ends missing
		REQUIRE( SeedHandlerExplicit::checkSeedEncoding( "2.||&2||" ) != "" );
		REQUIRE( SeedHandlerExplicit::checkSeedEncoding( "2||.&2||" ) != "" );
		REQUIRE( SeedHandlerExplicit::checkSeedEncoding( "2||&2.||" ) != "" );
		REQUIRE( SeedHandlerExplicit::checkSeedEncoding( "2||&2||." ) != "" );

	}

	SECTION( "test seed enumeration" ) {

		RnaSequence r1("r1", "GGGGGG");  // 1GGGGGG6
		RnaSequence r2("r2", "CCCCCG");  // 6GCCCCC1
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ
		SeedConstraint sC(3,0,0,0,0
				, AccessibilityDisabled::ED_UPPER_BOUND, 0
				, IndexRangeList("")
				, IndexRangeList("")
				, "1||&3||,2|.|&1||"
				, false, false, false);

		// create instance (triggers parsing)
		SeedHandlerExplicit sh( energy, sC );

		// check parsing number
		REQUIRE( sh.fillSeed(0,5,0,5) == 2 );

		// check seed iteration
		size_t i1=0, i2=0;
		REQUIRE( sh.updateToNextSeed(i1,i2) );
		bool isSeed13 = (i1 == 0 && i2 == 2);
		bool isSeed21 = (i1 == 1 && i2 == 4);
		bool isInputSeed = isSeed13 || isSeed21;
		REQUIRE( isInputSeed );
		REQUIRE( sh.updateToNextSeed(i1,i2) );
		bool isTheOtherSeed = (isSeed13 && (i1 == 1 && i2 == 4) ) || (isSeed21 && (i1 == 0 && i2 == 2) );
		REQUIRE( isTheOtherSeed );
		REQUIRE_FALSE( sh.updateToNextSeed(i1,i2) );

		// check seed init
		i1=20; i2=0;
		REQUIRE( sh.updateToNextSeed(i1,i2) );
		isSeed13 = (i1 == 0 && i2 == 2);
		isSeed21 = (i1 == 1 && i2 == 4);
		isInputSeed = isSeed13 || isSeed21;
		REQUIRE( isInputSeed );

		// check seed init
		i1=0; i2=20;
		REQUIRE( sh.updateToNextSeed(i1,i2) );
		isSeed13 = (i1 == 0 && i2 == 2);
		isSeed21 = (i1 == 1 && i2 == 4);
		isInputSeed = isSeed13 || isSeed21;
		REQUIRE( isInputSeed );

		// check seed init
		i1=20; i2=20;
		REQUIRE( sh.updateToNextSeed(i1,i2) );
		isSeed13 = (i1 == 0 && i2 == 2);
		isSeed21 = (i1 == 1 && i2 == 4);
		isInputSeed = isSeed13 || isSeed21;
		REQUIRE( isInputSeed );

		//////////////////  seeds of subregion  ///////////////////////////

		i1=0; i2=0;
		REQUIRE_FALSE( sh.updateToNextSeed(i1,i2,3,5,0,5) );
		REQUIRE_FALSE( sh.updateToNextSeed(i1,i2,0,5,0,1) );

		i1=0; i2=0;
		REQUIRE( sh.updateToNextSeed(i1,i2,0,0,2,2) );
		REQUIRE_FALSE( sh.updateToNextSeed(i1,i2,0,0,2,2) );

	}

}
