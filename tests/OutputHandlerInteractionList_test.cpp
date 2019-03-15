
#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/OutputHandlerInteractionList.h"
#include <stdexcept>

using namespace IntaRNA;

/**
 * Dummy class for testing the add result
 */
bool equalInteraction( const Interaction& i1 , const Interaction& i2 ) {
	return !(i1 < i2) && !(i2 < i1);
};

TEST_CASE( "OutputHandlerInteractionList", "[OutputHandlerInteractionList]" ) {

	// setup easylogging++ stuff if not already done
	#include "testEasyLoggingSetup.icc"

	RnaSequence r("test","AACCGGUU");

	OutputConstraint oc;

	SECTION("sorting") {

		Interaction i1(r,r);
		i1.basePairs.push_back( Interaction::BasePair(0,7));
		i1.basePairs.push_back( Interaction::BasePair(1,6));
		i1.energy = 0;
		REQUIRE( i1.isValid() );

		OutputHandlerInteractionList out(2);
		REQUIRE( out.empty() );

		out.add(i1,oc);
		REQUIRE( ! out.empty() );
		auto outIt = out.begin();
		REQUIRE( equalInteraction( **outIt, i1) );
		outIt++;
		REQUIRE( outIt == out.end() );

		Interaction i2(i1);
		i2.energy = -1;
		REQUIRE( i2.isValid() );

		// insert (should be new first element)
		out.add(i2,oc);
		REQUIRE( ! out.empty() );
		outIt = out.begin();
		REQUIRE( equalInteraction( **outIt, i2) );
		outIt++;
		REQUIRE( outIt != out.end() );
		REQUIRE( equalInteraction( **outIt, i1) );
		outIt++;
		REQUIRE( outIt == out.end() );

	}


	SECTION("double insertion") {

		Interaction i(r,r);
		i.basePairs.push_back( Interaction::BasePair(0,7));
		i.basePairs.push_back( Interaction::BasePair(1,6));
		i.energy = 0;
		REQUIRE( i.isValid() );

		OutputHandlerInteractionList out(2);
		REQUIRE( out.empty() );

		out.add(i,oc);
		REQUIRE( ! out.empty() );
		auto outIt = out.begin();
		REQUIRE( equalInteraction( **outIt, i) );
		outIt++;
		REQUIRE( outIt == out.end() );

		// insert a second time (should cause no insertion)
		out.add(i,oc);
		REQUIRE( ! out.empty() );
		outIt = out.begin();
		REQUIRE( equalInteraction( **outIt, i) );
		outIt++;
		REQUIRE( outIt == out.end() );

	}

}
