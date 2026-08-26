
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

	SECTION("C++ string boundary preserves ViennaRNA constraints") {

		const std::string seq = "GGGGAAAACCCC";
		RnaSequence rna("test", seq);
		AccessibilityConstraint constraint(rna, "xp........px", 0, "", "", "");
		VrnaHandler vrnaHandler(37, "Turner04", false, false);
		AccessibilityVrna acc(rna, 1, &constraint, vrnaHandler, rna.size());

		// 'x' is passed at both C-string boundaries and forces unpaired bases.
		REQUIRE( E_equal(acc.getED(0, 0), 0) );
		REQUIRE( E_equal(acc.getED(rna.size()-1, rna.size()-1), 0) );

		// Exercise the same x|........|x bytes through computeES() and
		// computeIntraEall(), including the reversed accessibility path.
		AccessibilityDisabled disabledAcc(rna, rna.size(), &constraint);
		ReverseAccessibility disabledAccReversed(disabledAcc);
		InteractionEnergyVrna energy(
				disabledAcc, disabledAccReversed, vrnaHandler, 16, 16, true);
		REQUIRE( E_isNotINF(energy.getEall1()) );
		REQUIRE( E_isNotINF(energy.getEall2()) );
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
