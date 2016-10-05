
#include "catch.hpp"

#undef NDEBUG

#include "RnaSequence.h"

TEST_CASE( "RnaSequence", "[RNAsequence]" ) {


	SECTION( "sequence with non-ACGU" ) {
		bool exceptionRaised = false;
		try {
			RnaSequence rna("test","ACUGACerror");
		} catch (std::exception & ex) {
			exceptionRaised = true;
		}
		REQUIRE( exceptionRaised );
	}

	SECTION( "getter", "[RNAsequence]" ) {

		std::string id= "test", seq="AAAUUUGGGCCC";

		RnaSequence rna(id, seq);

		// check data access
		REQUIRE( rna.size() == seq.size() );
		REQUIRE( rna.getId() == id );
		REQUIRE( rna.asString() == seq );
		REQUIRE( rna.asCodes().size() == seq.size() );
	}


	SECTION( "areComplementary()" ) {

		std::string id= "test", seq="AAAUUUGGGCCC";

		RnaSequence rna(id, seq);

		// check complementarity
		REQUIRE( RnaSequence::areComplementary( rna, rna, 0, 3) );
		REQUIRE_FALSE( RnaSequence::areComplementary( rna, rna, 0, 6) );
		REQUIRE_FALSE( RnaSequence::areComplementary( rna, rna, 0, 9) );
		REQUIRE( RnaSequence::areComplementary( rna, rna, 3, 6) );
		REQUIRE_FALSE( RnaSequence::areComplementary( rna, rna, 3, 9) );
		REQUIRE( RnaSequence::areComplementary( rna, rna, 6, 9) );

	}

}
