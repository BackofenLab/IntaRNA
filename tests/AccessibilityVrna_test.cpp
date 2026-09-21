
#include "catch.hpp"

#undef NDEBUG

#include <boost/numeric/ublas/io.hpp>
#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/AccessibilityVrna.h"
#include "IntaRNA/InteractionEnergyVrna.h"
#include "IntaRNA/ReverseAccessibility.h"

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


	SECTION("InteractionEnergyVrna respects accessibility base-pair span") {

		const std::string seq = "GGGGAAAACCCC";
		RnaSequence rna("test", seq);
		AccessibilityConstraint unrestrictedConstraint(rna.size(), 0, "", "", "");
		AccessibilityConstraint shortSpanConstraint(rna.size(), 3, "", "", "");
		AccessibilityDisabled unrestrictedAcc(rna, rna.size(), &unrestrictedConstraint);
		AccessibilityDisabled shortSpanAcc(rna, rna.size(), &shortSpanConstraint);
		ReverseAccessibility unrestrictedAccReversed(unrestrictedAcc);
		ReverseAccessibility shortSpanAccReversed(shortSpanAcc);
		VrnaHandler vrnaHandler(37, "Turner04", false, false);

		InteractionEnergyVrna unrestrictedEnergy(
				unrestrictedAcc, unrestrictedAccReversed, vrnaHandler, 16, 16, true);
		InteractionEnergyVrna shortSpanEnergy(
				shortSpanAcc, shortSpanAccReversed, vrnaHandler, 16, 16, true);

		REQUIRE( E_isNotINF(unrestrictedEnergy.getES1(0, rna.size()-1)) );
		REQUIRE( E_isINF(shortSpanEnergy.getES1(0, rna.size()-1)) );

		// Exercise both normal computeIntraEall() paths. Their ViennaRNA-owned
		// resources are released by the scoped owners in InteractionEnergyVrna.
		REQUIRE( E_isNotINF(unrestrictedEnergy.getEall1()) );
		REQUIRE( E_isNotINF(unrestrictedEnergy.getEall2()) );
	}
}
