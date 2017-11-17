
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

    REQUIRE( std::abs( acc.getED(0, 6) - 4.52824 ) < 1e-4 );
    REQUIRE( std::abs( acc.getED(0, 12) - 5.8415) < 1e-4 );
    REQUIRE( std::abs( acc.getED(1, 3) - 2.37302 ) < 1e-4 );
    REQUIRE( std::abs( acc.getED(1, 5) - 2.75135 ) < 1e-4 );
    REQUIRE( std::abs( acc.getED(3, 4) - 0.519897 ) < 1e-4 );
    REQUIRE( std::abs( acc.getED(4, 7) - 1.978 ) < 1e-4 );
    REQUIRE( std::abs( acc.getED(4, 12) - 5.8415) < 1e-4 );
    REQUIRE( std::abs( acc.getED(5, 5) - 0.0) < 1e-4 );
    REQUIRE( std::abs( acc.getED(5, 10) - 3.62722 ) < 1e-4 );
    REQUIRE( std::abs( acc.getED(6, 6) - 0.131359 ) < 1e-4 );
    REQUIRE( std::abs( acc.getED(10, 11) - 0.559891 ) < 1e-4 );
    REQUIRE( std::abs( acc.getED(11, 12) - 0.191058 ) < 1e-4 );

  }
}
