
#include "catch.hpp"

#undef NDEBUG

#include <boost/numeric/ublas/io.hpp>
#include "IntaRNA/AccessibilityVrna.h"

using namespace IntaRNA;



TEST_CASE("AccessibilityVrna", "[AccessibilityVrna]") {

	// setup easylogging++ stuff if not already done
	#include "testEasyLoggingSetup.icc"


  SECTION("ED test short") {

	const std::string seq = "GG";
	RnaSequence rna("test", seq);
	VrnaHandler vrnaHandler(37,"Turner04",false,false);
	AccessibilityVrna acc(rna, 2, NULL,vrnaHandler,2);

	REQUIRE( E_equal( acc.getED(0, 0), 0 ) );
	REQUIRE( E_equal( acc.getED(1, 1), 0 ) );
	REQUIRE( E_equal( acc.getED(0, 1), 0 ) );

  }
}
