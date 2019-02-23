
#include "catch.hpp"

#undef NDEBUG
#define protected public

#include "IntaRNA/RnaSequence.h"
#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/SeedHandlerMfe.h"
#include "IntaRNA/SeedHandlerIdxOffset.h"


using namespace IntaRNA;

TEST_CASE( "SeedHandlerMfe with offset", "[SeedHandlerIdxOffset]") {

	SECTION("SeedHandlerMfe: Case 1 - offset 1", "[SeedHandlerIdxOffset]") {
		RnaSequence r1("r1", "GGGGGG");
		RnaSequence r2("r2", "CCCCCG");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ
		SeedConstraint sC(3,0,0,0,0
				, AccessibilityDisabled::ED_UPPER_BOUND
				, IndexRangeList("")
				, IndexRangeList("")
				, "");

		SeedHandlerIdxOffset sHIO(new SeedHandlerMfe(energy, sC));

		sHIO.setOffset1(1);
		sHIO.setOffset2(1);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////   FILLSEED  ////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		REQUIRE(sHIO.fillSeed(0,energy.size1()-sHIO.getOffset1()-1, 0,energy.size2()-sHIO.getOffset2()-1) > 0);
	}
}