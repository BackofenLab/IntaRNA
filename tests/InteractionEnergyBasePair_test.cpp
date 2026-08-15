
#include "catch.hpp"

#undef NDEBUG

#include <cmath>
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/AccessibilityDisabled.h"

using namespace IntaRNA;

namespace {

class CountingEdEnergy : public InteractionEnergyBasePair {
public:
	CountingEdEnergy( const Accessibility & accS1
			, const ReverseAccessibility & accS2 )
		: InteractionEnergyBasePair( accS1, accS2, 1, 1, false, 1
				, Ekcal_2_E(-1.0), 1, Ekcal_2_E(0.0), false )
		, ed1(Ekcal_2_E(0.1))
		, ed2(Ekcal_2_E(0.2))
		, ed1Calls(0)
		, ed2Calls(0)
	{}

	E_type getED1( const size_t, const size_t ) const override {
		++ed1Calls;
		return ed1;
	}

	E_type getED2( const size_t, const size_t ) const override {
		++ed2Calls;
		return ed2;
	}

	void setED( const E_type newEd1, const E_type newEd2 ) {
		ed1 = newEd1;
		ed2 = newEd2;
	}

	void resetCalls() const {
		ed1Calls = 0;
		ed2Calls = 0;
	}

	E_type ed1;
	E_type ed2;
	mutable size_t ed1Calls;
	mutable size_t ed2Calls;
};

}

TEST_CASE( "InteractionEnergyBasePair", "[InteractionEnergyBasePair]" ) {

	// setup easylogging++ stuff if not already done
	#include "testEasyLoggingSetup.icc"

	RnaSequence rna("test","ACGU");

	AccessibilityDisabled acc(rna,rna.size(),NULL);
	ReverseAccessibility rAcc(acc);

	size_t maxLoop1 = 1, maxLoop2 = 2;
	InteractionEnergyBasePair energy( acc, rAcc, maxLoop1, maxLoop2, true, 1, Ekcal_2_E(-1.0), 1);

	SECTION("data access") {
		// check
		REQUIRE( &energy.getAccessibility1() == &acc );
		REQUIRE( &energy.getAccessibility2().getAccessibilityOrigin() == &acc );
		REQUIRE( energy.getMaxInternalLoopSize1() == maxLoop1 );
		REQUIRE( energy.getMaxInternalLoopSize2() == maxLoop2 );

		REQUIRE( energy.getAccessibility1().getSequence().asString().at(0) == 'A' );
		REQUIRE( energy.getAccessibility2().getSequence().asString().at(0) == 'U' );

	}

	SECTION("energy computation" ) {

		// dangling end check
		REQUIRE( energy.getE_danglingLeft(0,0) == 0.0 );
		REQUIRE( energy.getE_danglingLeft(2,2) == 0.0 );
		REQUIRE( energy.getE_danglingRight(0,0) == 0.0 );
		REQUIRE( energy.getE_danglingRight(2,2) == 0.0 );

		// dangling end check
		REQUIRE( energy.getE_danglingLeft(0,0) == 0.0 );
		REQUIRE( energy.getE_danglingLeft(2,2) == 0.0 );
		REQUIRE( energy.getE_danglingRight(0,0) == 0.0 );
		REQUIRE( energy.getE_danglingRight(2,2) == 0.0 );

		// init only, but interior loop called
		REQUIRE( energy.getE_interLeft( 0,0, 0,0 ) > 0.0 );

		// base pairs possible
		REQUIRE( energy.getE_interLeft( 0,1, 0,1 ) < 0.0 );

		// base pairs overlapping
		REQUIRE_FALSE( energy.getE_interLeft( 0,0, 0,1 ) < 0.0 );
		REQUIRE_FALSE( energy.getE_interLeft( 0,1, 0,0 ) < 0.0 );

		// base pairs not possible
		REQUIRE_FALSE( energy.getE_interLeft( 0,0, 1,1 ) < 0.0 );
		REQUIRE_FALSE( energy.getE_interLeft( 0,1, 2,3 ) < 0.0 );
		REQUIRE_FALSE( energy.getE_interLeft( 0,1, 0,2 ) < 0.0 );
		REQUIRE_FALSE( energy.getE_interLeft( 0,2, 1,2 ) < 0.0 );

		// loop size exceeded
		REQUIRE_FALSE( energy.getE_interLeft( 0,3, 1,2 ) < 0.0 );
		REQUIRE_FALSE( energy.getE_interLeft( 1,2, 0,3 ) < 0.0 );
		REQUIRE_FALSE( energy.getE_interLeft( 0,3, 0,3 ) < 0.0 );

	}

  SECTION("ES computation") {
    // ES covers structures containing at least one intramolecular base pair.
    // ACGU has exactly one admissible pair of weight exp(1), hence ES=-1.
    REQUIRE( E_equal(energy.getES1(0, 3), Ekcal_2_E(-1.0)) );
    REQUIRE( E_equal(energy.getES2(0, 3), Ekcal_2_E(-1.0)) );
    REQUIRE( E_isINF(energy.getES1(0, 2)) );
    REQUIRE( E_isINF(energy.getES1(1, 2)) );
    REQUIRE( E_isINF(energy.getES2(0, 2)) );
    REQUIRE( E_isINF(energy.getES2(1, 2)) );
  }

	SECTION("overall energy reads each accessibility penalty once") {
		CountingEdEnergy countingEnergy(acc, rAcc);
		const E_type hybridE = Ekcal_2_E(-1.0);

		REQUIRE( countingEnergy.getE(0, 0, 0, 0, hybridE)
				== hybridE + Ekcal_2_E(0.1) + Ekcal_2_E(0.2) );
		REQUIRE( countingEnergy.ed1Calls == 1 );
		REQUIRE( countingEnergy.ed2Calls == 1 );

		countingEnergy.resetCalls();
		REQUIRE( E_isINF(countingEnergy.getE(0, 0, 0, 0, E_INF)) );
		REQUIRE( countingEnergy.ed1Calls == 0 );
		REQUIRE( countingEnergy.ed2Calls == 0 );

		countingEnergy.setED(Accessibility::ED_UPPER_BOUND, Ekcal_2_E(0.2));
		countingEnergy.resetCalls();
		REQUIRE( E_isINF(countingEnergy.getE(0, 0, 0, 0, hybridE)) );
		REQUIRE( countingEnergy.ed1Calls == 1 );
		REQUIRE( countingEnergy.ed2Calls == 0 );

		countingEnergy.setED(Ekcal_2_E(0.1), Accessibility::ED_UPPER_BOUND);
		countingEnergy.resetCalls();
		REQUIRE( E_isINF(countingEnergy.getE(0, 0, 0, 0, hybridE)) );
		REQUIRE( countingEnergy.ed1Calls == 1 );
		REQUIRE( countingEnergy.ed2Calls == 1 );
	}

}
