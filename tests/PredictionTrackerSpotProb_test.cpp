
#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/PredictionTrackerSpotProb.h"
#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/ReverseAccessibility.h"
#include "IntaRNA/InteractionEnergyBasePair.h"

#include <stdexcept>

using namespace IntaRNA;


TEST_CASE( "PredictionTrackerSpotProb", "[PredictionTrackerSpotProb]" ) {

	// setup easylogging++ stuff if not already done
	#include "testEasyLoggingSetup.icc"

	// setup dummy data
	RnaSequence r1("r1","AACCG");
	RnaSequence r2("r2","AGUUNNNN");
	AccessibilityDisabled acc1(r1, 0, NULL);
	AccessibilityDisabled acc2(r2, 0, NULL);
	ReverseAccessibility racc( acc2 );
	InteractionEnergyBasePair energy( acc1, racc );

	SECTION("empty tracking") {
		// output streams
		std::stringstream out;

		// create and directly destroy
		PredictionTrackerSpotProb * tracker = new PredictionTrackerSpotProb( energy, "", out);
		// destroy to flush output
		delete tracker; tracker = NULL;

//		LOG(DEBUG)<<"empty tracking: "<<out.str();

		// check output
		REQUIRE( out.str() == "spot;probability\n0&0;0\n" );
	}

	SECTION("empty tracking - no update") {
		// output streams
		std::stringstream out;

		// create and directly destroy
		PredictionTrackerSpotProb * tracker = new PredictionTrackerSpotProb( energy, "1&1", out);
		// destroy to flush output
		delete tracker; tracker = NULL;

//		LOG(DEBUG)<<"empty tracking: "<<out.str();

		// check output
		REQUIRE( out.str() == "spot;probability\n0&0;0\n1&1;0\n" );
	}

	SECTION("single update - no match") {
		// output streams
		std::stringstream out;

		// create
		PredictionTrackerSpotProb * tracker = new PredictionTrackerSpotProb( energy, "1&1", out);
		// add range
//		AACCG
//		01234567
//		NNNNUUGA.AGUUNNNN
		tracker->updateOptimumCalled( 0,1, 4,5, 2.0 );
		// destroy to flush output
		delete tracker; tracker = NULL;

		// check output
		REQUIRE( out.str() == "spot;probability\n0&0;1\n1&1;0\n" );
	}


	SECTION("single update - match") {
		// output streams
		std::stringstream out;

		// create
		PredictionTrackerSpotProb * tracker = new PredictionTrackerSpotProb( energy, "1&4", out);
		// add range
		//		AACCG
		//		01234567
		//		NNNNUUGA.AGUUNNNN
		tracker->updateOptimumCalled( 0,1, 4,5, 2.0 );
		// destroy to flush output
		delete tracker; tracker = NULL;

		// check output
		REQUIRE( out.str() == "spot;probability\n0&0;0\n1&4;1\n" );
	}


	SECTION("single update - mix") {
		// output streams
		std::stringstream out;

		// create
		PredictionTrackerSpotProb * tracker = new PredictionTrackerSpotProb( energy, "1&3,5&4", out);
		// add range
		//		AACCG
		//		01234567 01234567
		//		NNNNUUGA.AGUUNNNN
		tracker->updateOptimumCalled( 0,1, 4,5, 2.0 ); // 1&3
		// destroy to flush output
		delete tracker; tracker = NULL;

		// check output
		REQUIRE( out.str() == "spot;probability\n0&0;0\n1&3;1\n5&4;0\n" );
	}

	SECTION("double update - mix") {
		// output streams
		std::stringstream out;

		// create
		PredictionTrackerSpotProb * tracker = new PredictionTrackerSpotProb( energy, "1&3,5&4", out);
		// add range
		//		AACCG
		//		01234567 01234567
		//		NNNNUUGA.AGUUNNNN
		tracker->updateOptimumCalled( 0,1, 4,5, 2.0 ); // 1&3
		tracker->updateOptimumCalled( 4,4, 5,5, 2.0 ); // none
		// destroy to flush output
		delete tracker; tracker = NULL;

		// check output
		REQUIRE( out.str() == "spot;probability\n0&0;0.5\n1&3;0.5\n5&4;0\n" );
	}

}
