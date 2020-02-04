
#include "catch.hpp"

#undef NDEBUG
#define protected public

#include "IntaRNA/RnaSequence.h"
#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/SeedHandlerNoBulge.h"
#include "IntaRNA/SeedHandlerIdxOffset.h"

#include <set>

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
				, 0
				, IndexRangeList("")
				, IndexRangeList("")
				, ""
				, false, false, false );

		SeedHandlerNoBulge shOrig(energy, sC);
		SeedHandlerIdxOffset shIO( &shOrig, false );

		shIO.setOffset1(1);
		shIO.setOffset2(1);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////   FILLSEED  ////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// we omit 3 of 12 possible seeds (with i1==0)
		REQUIRE(shIO.fillSeed(0,energy.size1()-shIO.getOffset1()-1, 0,energy.size2()-shIO.getOffset2()-1) == 9);

		size_t i1=10, i2=10;

		// generate sorted list of seed starts without offset
		std::set< Interaction::BasePair > ioSeeds;
		while (shIO.updateToNextSeed(i1,i2)) {
			ioSeeds.insert( Interaction::BasePair(i1+shIO.getOffset1(),i2+shIO.getOffset2()) );
		}

		i1=10; i2=10;
		size_t count = 0;
		// check for list of seed starts without offset
		while (shOrig.updateToNextSeed(i1,i2, shIO.getOffset1(), 10, shIO.getOffset2(), 10)) {
			count++;
			REQUIRE( ioSeeds.find( Interaction::BasePair(i1,i2) ) != ioSeeds.end() );
		}
		REQUIRE( count == 9 );

	}
}
