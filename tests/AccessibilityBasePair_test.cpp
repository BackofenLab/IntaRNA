
#include "catch.hpp"

#undef NDEBUG

#include <boost/numeric/ublas/io.hpp>
#include "IntaRNA/AccessibilityBasePair.h"

using namespace IntaRNA;

const std::string seq = "gguccacguccaa";


TEST_CASE("AccessibilityBasePair", "[AccessibilityBasePair]") {

	// setup easylogging++ stuff if not already done
	#include "testEasyLoggingSetup.icc"

  RnaSequence rna("test", seq);

  SECTION("Basepair Energy difference") {
    AccessibilityBasePair acc(rna, 10, NULL);

//    std::cout <<"\n\n########\n";
//    for (int i =0; i<acc.getSequence().size(); i++) {
//    	std::cout <<(i);
//    	for (int j=0; j <acc.getSequence().size(); j++) {
//    		if (j<i) {
//    			std::cout <<" ---";
//    		} else {
//    			std::cout <<" "<<acc.getED(i,j);
//    		}
//    	}
//    	std::cout <<std::endl;
//    }
//    std::cout <<"########\n";

    REQUIRE( E_equal( acc.getED(0, 6), 452 ) );
    REQUIRE( E_equal( acc.getED(0, 12), 584 ) );
    REQUIRE( E_equal( acc.getED(1, 3), 237 ) );
    REQUIRE( E_equal( acc.getED(1, 5), 275 ) );
    REQUIRE( E_equal( acc.getED(3, 4), 51 ) );
    REQUIRE( E_equal( acc.getED(4, 7), 197 ) );
    REQUIRE( E_equal( acc.getED(4, 12), 584 ) );
    REQUIRE( E_equal( acc.getED(5, 5), 0) );
    REQUIRE( E_equal( acc.getED(5, 10), 362 ) );
    REQUIRE( E_equal( acc.getED(6, 6), 13 ) );
    REQUIRE( E_equal( acc.getED(10, 11), 55 ) );
    REQUIRE( E_equal( acc.getED(11, 12), 19 ) );

  }
}
