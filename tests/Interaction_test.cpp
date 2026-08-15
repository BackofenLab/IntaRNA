

#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/Interaction.h"
#include "IntaRNA/InteractionRange.h"

using namespace IntaRNA;

TEST_CASE( "Interaction", "[Interaction]" ) {

	// setup easylogging++ stuff if not already done
	#include "testEasyLoggingSetup.icc"

	RnaSequence r("test","AACCGGUU");

	SECTION("empty construction") {

		Interaction inter( r, r );

		REQUIRE( inter.isEmpty() );
		REQUIRE( inter.s1 == &r );
		REQUIRE( inter.s2 == &r );
		REQUIRE_FALSE( inter.isValid() );
	}

	SECTION("add base pairs and sort") {

		Interaction inter( r, r );

		inter.basePairs.push_back( Interaction::BasePair( 1, 6 ) );
		REQUIRE_FALSE( inter.isEmpty() );
		REQUIRE( inter.isValid() );
		REQUIRE( Interaction::dotBracket(inter) == "(&)");
		REQUIRE( Interaction::dotBar(inter) == "2|&7|");

		inter.basePairs.push_back( Interaction::BasePair( 0, 7 ) );
		REQUIRE_FALSE( inter.isEmpty() );
		REQUIRE_FALSE( inter.isValid() );

		inter.sort();
		REQUIRE( inter.isValid() );
		REQUIRE( Interaction::dotBracket(inter) == "((&))");
		REQUIRE( Interaction::dotBar(inter) == "1||&7||");

		inter.basePairs.push_back( Interaction::BasePair( 3, 4 ) );
		REQUIRE( inter.isValid() );
		REQUIRE( Interaction::dotBracket(inter) == "((.(&).))");
		REQUIRE( Interaction::dotBar(inter) == "1||.|&5|.||");

	}


	SECTION("init with interaction range") {

		InteractionRange range(r,r, IndexRange(0,2), IndexRange(7,5));
		REQUIRE( range.isSane() );

		Interaction inter(range);
		REQUIRE_FALSE( inter.isEmpty() );
		REQUIRE( inter.isValid() );
		REQUIRE( inter.basePairs.size() == 2 );
		REQUIRE( inter.basePairs.at(0).first == 0 );
		REQUIRE( inter.basePairs.at(0).second == 7 );
		REQUIRE( inter.basePairs.at(1).first == 2 );
		REQUIRE( inter.basePairs.at(1).second == 5 );

	}

	SECTION("assign interaction range") {

		InteractionRange range(r,r, IndexRange(0,2), IndexRange(7,5));
		REQUIRE( range.isSane() );

		// empty interaction
		Interaction inter(r,r);
		REQUIRE( inter.isEmpty() );

		// assignment operator
		inter = range;
		REQUIRE_FALSE( inter.isEmpty() );
		REQUIRE( inter.isValid() );
		REQUIRE( inter.basePairs.size() == 2 );
		REQUIRE( inter.basePairs.at(0).first == 0 );
		REQUIRE( inter.basePairs.at(0).second == 7 );
		REQUIRE( inter.basePairs.at(1).first == 2 );
		REQUIRE( inter.basePairs.at(1).second == 5 );

	}

	SECTION("self assignment preserves owned state") {

		Interaction inter(r,r);
		inter.basePairs.push_back( Interaction::BasePair(0,7) );
		inter.basePairs.push_back( Interaction::BasePair(1,6) );
		inter.energy = Ekcal_2_E(-2.0);
		inter.seed = new Interaction::SeedSet();
		inter.seed->insert( Interaction::Seed(
				Interaction::BasePair(0,7), Interaction::BasePair(1,6), inter.energy) );

		inter = inter;

		REQUIRE( inter.basePairs.size() == 2 );
		REQUIRE( inter.basePairs.front() == Interaction::BasePair(0,7) );
		REQUIRE( inter.basePairs.back() == Interaction::BasePair(1,6) );
		REQUIRE( inter.energy == Ekcal_2_E(-2.0) );
		REQUIRE( inter.seed != NULL );
		REQUIRE( inter.seed->size() == 1 );
	}

	SECTION("seeded and unseeded interactions compare safely") {

		Interaction unseeded(r,r);
		unseeded.basePairs.push_back( Interaction::BasePair(0,7) );
		unseeded.energy = Ekcal_2_E(-1.0);

		Interaction seeded(unseeded);
		seeded.seed = new Interaction::SeedSet();
		seeded.seed->insert( Interaction::Seed(
				Interaction::BasePair(0,7), Interaction::BasePair(0,7), seeded.energy) );

		REQUIRE_FALSE( unseeded == seeded );
		REQUIRE_FALSE( seeded == unseeded );
	}

	SECTION("seed ordering preserves distinct equal-energy ranges") {

		const E_type seedEnergy = -100;
		const Interaction::Seed seed(
				Interaction::BasePair(1,6), Interaction::BasePair(3,4), seedEnergy );
		const Interaction::Seed differentSeq2Right(
				Interaction::BasePair(1,7), Interaction::BasePair(3,4), seedEnergy );
		const Interaction::Seed differentSeq1Right(
				Interaction::BasePair(1,6), Interaction::BasePair(4,4), seedEnergy );
		const Interaction::Seed differentSeq2Left(
				Interaction::BasePair(1,6), Interaction::BasePair(3,3), seedEnergy );
		Interaction::SeedSet seeds;

		REQUIRE( seeds.insert(seed).second );
		REQUIRE( seeds.insert(differentSeq2Right).second );
		REQUIRE( seeds.insert(differentSeq1Right).second );
		REQUIRE( seeds.insert(differentSeq2Left).second );
		REQUIRE( seeds.size() == 4 );

		// Only an exact duplicate is equivalent in the ordering.
		REQUIRE_FALSE( seeds.insert(seed).second );
		REQUIRE( seeds.size() == 4 );
	}

}
