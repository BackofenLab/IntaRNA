
#include "catch.hpp"

#undef NDEBUG
#define protected public

#include "IntaRNA/RnaSequence.h"
#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/SeedHandlerMfe.h"


using namespace IntaRNA;

TEST_CASE( "SeedHandlerMfe", "[SeedHandlerMfe]") {

	SECTION("SeedHandlerMfe: Case 1 - offset 1", "[SeedHandlerMfe]") {
		RnaSequence r1("r1", "GGGGG");
		RnaSequence r2("r2", "CCCCC");
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
				, false, false, false
				);

		SeedHandlerMfe sHM(energy, sC);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////   FILLSEED  ////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		REQUIRE(sHM.fillSeed(0,energy.size1()-1, 0,energy.size2()-1) > 0);
	}

	SECTION("SeedHandlerMfe: with lp", "[SeedHandlerMfe]") {
		RnaSequence r1("r1", "GGCGCGG");
		RnaSequence r2("r2", "CCCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ
		SeedConstraint sC(4,3,3,0,0
				, AccessibilityDisabled::ED_UPPER_BOUND
				, 0
				, IndexRangeList("")
				, IndexRangeList("")
				, ""
				, false, false, false
				);

		SeedHandlerMfe sHM(energy, sC);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////   FILLSEED  ////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		REQUIRE(sHM.fillSeed(0,energy.size1()-1, 0,energy.size2()-1) == 4); // getting only one seed per si, thus 2x2

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////   TRACEBACK  ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		size_t si1=9999999, si2=9999999;
		while( sHM.updateToNextSeed(si1,si2,0,r1.size()-1,0,r2.size()) ) {
			// prepare trace
			Interaction i(r1,r2);
			i.basePairs.push_back(energy.getBasePair(si1,si2));
			i.basePairs.push_back(energy.getBasePair(si1+sHM.getSeedLength1(si1,si2)-1,si2+sHM.getSeedLength2(si1,si2)-1));
			// trace base pairs
			sHM.traceBackSeed(i, si1, si2);
			// check number of traced base pairs
			REQUIRE(i.basePairs.size() == 4);
		}

	}

	SECTION("SeedHandlerMfe: no lp - 1", "[SeedHandlerMfe]") {
		RnaSequence r1("r1", "GCGCGG");
		RnaSequence r2("r2", "CCCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		{
			size_t bp = 2;
			// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ
			SeedConstraint sC(bp,3,3,0,0
					, AccessibilityDisabled::ED_UPPER_BOUND
					, 0
					, IndexRangeList("")
					, IndexRangeList("")
					, ""
					, false, false, true
					);

			SeedHandlerMfe sHM(energy, sC);

			//////////////////////////////////////////////////////////////////////////////////////////////////////////
			/////////////////////////////////////////////   FILLSEED  ////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////////////////////////
			REQUIRE(sHM.fillSeed(0,energy.size1()-1, 0,energy.size2()-1) == 4);
		}

		for (size_t bp = 3; bp <= 4; bp++) {
			// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ
			SeedConstraint sC(bp,3,3,0,0
					, AccessibilityDisabled::ED_UPPER_BOUND
					, 0
					, IndexRangeList("")
					, IndexRangeList("")
					, ""
					, false, false, true
					);

			SeedHandlerMfe sHM(energy, sC);

			//////////////////////////////////////////////////////////////////////////////////////////////////////////
			/////////////////////////////////////////////   FILLSEED  ////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////////////////////////
			REQUIRE(sHM.fillSeed(0,energy.size1()-1, 0,energy.size2()-1) == 0);


			//////////////////////////////////////////////////////////////////////////////////////////////////////////
			/////////////////////////////////////////////   TRACEBACK  ///////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////////////////////////
			size_t si1=9999999, si2=9999999;
			while( sHM.updateToNextSeed(si1,si2,0,r1.size()-1,0,r2.size()) ) {
				// prepare trace
				Interaction i(r1,r2);
				i.basePairs.push_back(energy.getBasePair(si1,si2));
				i.basePairs.push_back(energy.getBasePair(si1+sHM.getSeedLength1(si1,si2)-1,si2+sHM.getSeedLength2(si1,si2)-1));
				// trace base pairs
				sHM.traceBackSeed(i, si1, si2);
				// check number of traced base pairs
				REQUIRE(i.basePairs.size() == bp);
			}

		}
	}

	SECTION("SeedHandlerMfe: no lp - 2", "[SeedHandlerMfe]") {
		RnaSequence r1("r1", "GGCGCGG");
		RnaSequence r2("r2", "CCCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		size_t seedBP=4;
		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ
		SeedConstraint sC(seedBP,3,3,0,0
				, AccessibilityDisabled::ED_UPPER_BOUND
				, 0
				, IndexRangeList("")
				, IndexRangeList("")
				, ""
				, false, false, true
				);

		SeedHandlerMfe sHM(energy, sC);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////   FILLSEED  ////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		REQUIRE(sHM.fillSeed(0,energy.size1()-1, 0,energy.size2()-1) == 2);


		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////   TRACEBACK  ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		size_t si1=9999999, si2=9999999;
		while( sHM.updateToNextSeed(si1,si2,0,r1.size()-1,0,r2.size()) ) {
			// prepare trace
			Interaction i(r1,r2);
			i.basePairs.push_back(energy.getBasePair(si1,si2));
			i.basePairs.push_back(energy.getBasePair(si1+sHM.getSeedLength1(si1,si2)-1,si2+sHM.getSeedLength2(si1,si2)-1));
			// trace base pairs
			sHM.traceBackSeed(i, si1, si2);
			// check number of traced base pairs
			REQUIRE(i.basePairs.size() == seedBP);
		}
	}

	SECTION("reused handler matches fresh handlers after range shape changes") {
		RnaSequence r1("r1", "GGGGGGGG");
		RnaSequence r2("r2", "CCCCCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);
		SeedConstraint sC(3, 2, 2, 2, 0,
				AccessibilityDisabled::ED_UPPER_BOUND, 0,
				IndexRangeList(""), IndexRangeList(""), "",
				false, false, false);
		SeedHandlerMfe reused(energy, sC);

		auto compareWithFresh = [&](const IndexRange & range1,
				const IndexRange & range2) {
			CAPTURE(range1.from, range1.to, range2.from, range2.to);
			SeedHandlerMfe fresh(energy, sC);
			REQUIRE(reused.fillSeed(range1.from, range1.to,
					range2.from, range2.to)
					== fresh.fillSeed(range1.from, range1.to,
							range2.from, range2.to));
			for (size_t i1 = range1.from; i1 <= range1.to; ++i1) {
				for (size_t i2 = range2.from; i2 <= range2.to; ++i2) {
					CAPTURE(i1, i2);
					REQUIRE(reused.getSeedE(i1, i2)
							== fresh.getSeedE(i1, i2));
					REQUIRE(reused.getSeedLength1(i1, i2)
							== fresh.getSeedLength1(i1, i2));
					REQUIRE(reused.getSeedLength2(i1, i2)
							== fresh.getSeedLength2(i1, i2));
				}
			}
		};

		compareWithFresh(IndexRange(0, 5), IndexRange(0, 4)); // 6 x 5
		compareWithFresh(IndexRange(1, 6), IndexRange(1, 5)); // same shape
		compareWithFresh(IndexRange(2, 5), IndexRange(2, 4)); // smaller: 4 x 3
		compareWithFresh(IndexRange(3, 5), IndexRange(1, 4)); // same product: 3 x 4
	}


}
