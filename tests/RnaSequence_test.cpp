
#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/RnaSequence.h"

using namespace IntaRNA;

TEST_CASE( "RnaSequence", "[RNAsequence]" ) {

	// setup easylogging++ stuff if not already done
	#include "testEasyLoggingSetup.icc"


	SECTION( "sequence with non-ACGU" ) {
		bool exceptionRaised = ! RnaSequence::isValidSequenceIUPAC("ACUGAC_error");
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

	SECTION( "operator==()" ) {

		std::string id= "test", seq="AAAUUUGGGCCC";

		RnaSequence rna(id, seq);

		// check equality
		REQUIRE( rna == rna );
		REQUIRE( rna == RnaSequence(id,seq) );
		REQUIRE( RnaSequence(id,seq) == rna );
		REQUIRE_FALSE( rna == RnaSequence(id,seq+"A") );
		REQUIRE_FALSE( RnaSequence(id,seq+"A") == rna );
		REQUIRE_FALSE( rna == RnaSequence(id+"2",seq) );
		REQUIRE_FALSE( RnaSequence(id+"2",seq) == rna );
		REQUIRE_FALSE( rna == RnaSequence(id+"2",seq+"A") );
		REQUIRE_FALSE( RnaSequence(id+"2",seq+"A") == rna );
	}

}
