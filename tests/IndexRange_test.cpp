


#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/IndexRange.h"

using namespace IntaRNA;

#include <sstream>

TEST_CASE( "IndexRange", "[IndexRange]" ) {


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

}

