

#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/AccessibilityConstraint.h"

using namespace IntaRNA;

TEST_CASE( "AccessibilityConstraint", "[AccessibilityConstraint]" ) {

	// setup easylogging++ stuff if not already done
	#include "testEasyLoggingSetup.icc"

	const size_t seqLength = 12;

	SECTION("check empty construction") {
		AccessibilityConstraint c(seqLength,0,"","","");
		REQUIRE( c.isEmpty() );
		REQUIRE_FALSE( c.isMarkedBlocked(0) );
		REQUIRE_FALSE( c.isMarkedAccessible(0) );
		REQUIRE( c.isAccessible(0) );
		REQUIRE( c.isUnconstrained(0) );
	}

	SECTION("check regular expressions") {
		boost::regex regex;

		// general test
		regex = boost::regex(R"(\s)");
		REQUIRE_FALSE(boost::regex_match( "", regex, boost::match_perl ));
		REQUIRE(boost::regex_match( " ", regex, boost::match_perl ));

		// test dot-bracket alphabet
		regex = boost::regex("["+AccessibilityConstraint::dotBracketAlphabet+"]+");
		REQUIRE_FALSE(boost::regex_match( "", regex, boost::match_perl ));
		REQUIRE(boost::regex_match( AccessibilityConstraint::dotBracketAlphabet, regex, boost::match_perl ));

		// test dot-bracket alphabet
		regex = AccessibilityConstraint::regex;
		REQUIRE_FALSE(boost::regex_match( "", regex, boost::match_perl ));
		REQUIRE(boost::regex_match( AccessibilityConstraint::dotBracketAlphabet, regex, boost::match_perl ));
		REQUIRE(boost::regex_match( "b:3-4", regex, boost::match_perl ));
		REQUIRE(boost::regex_match( "b:3-4,11-12", regex, boost::match_perl ));
		REQUIRE(boost::regex_match( "x:11-12", regex, boost::match_perl ));
		REQUIRE(boost::regex_match( "x:3-4,11-12", regex, boost::match_perl ));
		REQUIRE(boost::regex_match( "p:11-12", regex, boost::match_perl ));
		REQUIRE(boost::regex_match( "p:3-4,11-12", regex, boost::match_perl ));
		REQUIRE(boost::regex_match( "b:3-4,11-12,x:7-8,p:9-9", regex, boost::match_perl ));
		REQUIRE(boost::regex_match( "b:3-4,11-12,x:7--8,p:9-9", regex, boost::match_perl ));

		REQUIRE_FALSE(boost::regex_match( "b:3-4d,11-12,x:7-8,p:9-9", regex, boost::match_perl ));
		REQUIRE_FALSE(boost::regex_match( "b:3-4,11-12,f:7-8,p:9-9", regex, boost::match_perl ));
		REQUIRE_FALSE(boost::regex_match( "b:3-4.0,11-12,x:7-8,p:9-9", regex, boost::match_perl ));
		REQUIRE_FALSE(boost::regex_match( "b:3-4,11,6-12,x:7-8,p:9-9", regex, boost::match_perl ));
		REQUIRE_FALSE(boost::regex_match( "b:,x:7-8,p:9-9", regex, boost::match_perl ));
	}

	SECTION("check dot-bracket construction") {
		RnaSequence rna("tmp",   "AAAAAAAAAAAA");
		std::string constraint = "..bb..xxp.bb";
		AccessibilityConstraint c(rna,constraint,0,"","","");

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


	SECTION("check range-based construction") {
		RnaSequence rna("tmp",   "AAAAAAAAAAAA");
		std::string constraint = "b:3-4,11-12,x:7-8,p:9-9";
		AccessibilityConstraint c(rna,constraint,0,"","","");

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
		for (size_t i=0; i<seqLength; i++) {
			vrnaStyle += toString(c.getVrnaDotBracket(i));
		}
		REQUIRE( vrnaStyle == "..xx..xx|.xx");
	}


}

