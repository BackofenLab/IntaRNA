
#include "catch.hpp"

#undef NDEBUG

#include <boost/numeric/ublas/io.hpp>
#include "IntaRNA/AccessibilityBasePair.h"
#include <iostream>

const std::string seq = "gguccac";

#include <sstream>

TEST_CASE("AccessibilityBasePair", "[AccessibilityBasePair]") {

  RnaSequence rna("test", seq);

  SECTION("Basepair Energy difference") {
    AccessibilityBasePair acc(rna, 10, NULL);

    REQUIRE( abs( acc.getED(1, 5) - 2.7806 ) < 1e-4 );
    REQUIRE( abs( acc.getED(3, 4) - 0.8074 ) < 1e-4 );
    REQUIRE( abs( acc.getED(1, 3) - 2.2256 ) < 1e-4 );
    REQUIRE( abs( acc.getED(0, 6) - 4.0745 ) < 1e-4 );
    REQUIRE( abs( acc.getED(6, 6) - 0.3467 ) < 1e-4 );

  }
}
