
#include "catch.hpp"

#undef NDEBUG

#include <boost/numeric/ublas/io.hpp>
#include "IntaRNA/AccessibilityBasePair.h"

using namespace IntaRNA;

const std::string seq = "gguccacguccaa";


TEST_CASE("AccessibilityBasePair", "[AccessibilityBasePair]") {

  RnaSequence rna("test", seq);

  SECTION("Basepair Energy difference") {
    AccessibilityBasePair acc(rna, 10, NULL);

    REQUIRE( std::abs( acc.getED(0, 6) - 4.50986 ) < 1e-4 );
    REQUIRE( std::abs( acc.getED(0, 12) - 5.809142) < 1e-4 );
    REQUIRE( std::abs( acc.getED(1, 3) - 2.37515 ) < 1e-4 );
    REQUIRE( std::abs( acc.getED(1, 5) - 2.74887 ) < 1e-4 );
    REQUIRE( std::abs( acc.getED(3, 4) - 0.51919 ) < 1e-4 );
    REQUIRE( std::abs( acc.getED(4, 7) - 1.98050 ) < 1e-4 );
    REQUIRE( std::abs( acc.getED(4, 12) - 5.809142) < 1e-4 );
    REQUIRE( std::abs( acc.getED(5, 5) - 0.0) < 1e-4 );
    REQUIRE( std::abs( acc.getED(5, 10) - 3.6119 ) < 1e-4 );
    REQUIRE( std::abs( acc.getED(6, 6) - 0.13124 ) < 1e-4 );
    REQUIRE( std::abs( acc.getED(10, 11) - 0.56036 ) < 1e-4 );
    REQUIRE( std::abs( acc.getED(11, 12) - 0.19116 ) < 1e-4 );

  }
}
