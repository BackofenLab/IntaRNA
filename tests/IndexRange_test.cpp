


#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/IndexRange.h"

using namespace IntaRNA;

#include <sstream>

TEST_CASE( "IndexRange", "[IndexRange]" ) {

	// setup easylogging++ stuff if not already done
	#include "testEasyLoggingSetup.icc"


	IndexRange range;

	SECTION("check default construction") {
		REQUIRE( range.isAscending() );
		REQUIRE( range.from == 0 );
		REQUIRE( range.to > range.from );
		REQUIRE_FALSE( range.isDescending() );
		REQUIRE( range.to > std::numeric_limits<size_t>::max()-1 );
	}

	SECTION("check ordering") {

		IndexRange r2;

		// equal
		REQUIRE_FALSE( range < r2 );

		// equal from
		r2.to = 10;
		REQUIRE_FALSE( range < r2 );
		REQUIRE( r2 < range );
		REQUIRE( r2.isAscending() );
		REQUIRE_FALSE( r2.isDescending() );

		// from ordered
		r2.from = 10;
		REQUIRE( range < r2 );
		REQUIRE( r2.isDescending() );
		REQUIRE( r2.isAscending() );
	}


	SECTION("check string en-/decoding") {

		// write to string
		std::stringstream s;
		range = (4,8);
		s <<range;

		// parse from string
		IndexRange r2(s.str());

		// check
		REQUIRE( range == r2 );
	}
	
	SECTION("check overlapping Windows") {
		
		// first example: starting with index 0, cutting into windows works out exactly
		size_t windowWidth = 10;
		size_t windowsOverlap = 3;
		IndexRange r3 = IndexRange(0, 23);
		IndexRange rw1 = IndexRange(0, 9);
		IndexRange rw2 = IndexRange(7, 16);
		IndexRange rw3 = IndexRange(14, 23);
		std::vector<IndexRange> windows = r3.overlappingWindows(windowWidth, windowsOverlap);
		
		REQUIRE ( windows.size() == 3 );
		REQUIRE ( windows[0] == rw1 );
		REQUIRE ( windows[1] == rw2 );
		REQUIRE ( windows[2] == rw3 );
		
		// second example: starting with index 123, the last window is smaller than windowWidth
		IndexRange r4 = IndexRange(123, 151);
		IndexRange rw4 = IndexRange(123, 132);
		IndexRange rw5 = IndexRange(130, 139);
		IndexRange rw6 = IndexRange(137, 146);
		IndexRange rw7 = IndexRange(144, 151);
		windows = r4.overlappingWindows(windowWidth, windowsOverlap);
		
		REQUIRE ( windows.size() == 4 );
		REQUIRE ( windows[0] == rw4 );
		REQUIRE ( windows[1] == rw5 );
		REQUIRE ( windows[2] == rw6 );
		REQUIRE ( windows[3] == rw7 );
		
		// third example: one window is enough for the whole IndexRange
		IndexRange r5 = IndexRange(453, 458);
		windows = r5.overlappingWindows(windowWidth, windowsOverlap);
		
		REQUIRE ( windows.size() == 1 );
		REQUIRE ( windows[0] == r5 );
	}
	
	
}

