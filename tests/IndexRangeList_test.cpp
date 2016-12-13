


#include "catch.hpp"

#undef NDEBUG

#include "IndexRangeList.h"

TEST_CASE( "IndexRangeList", "[IndexRangeList]" ) {


	IndexRangeList rangeList;

	SECTION("check construction and empty list") {
		REQUIRE( rangeList.size() == 0 );
		REQUIRE( rangeList.empty() );
		REQUIRE( rangeList.begin() == rangeList.end() );
		REQUIRE( rangeList.rbegin() == rangeList.rend() );
		REQUIRE_FALSE( rangeList.covers( 3 ) );
		REQUIRE_FALSE( rangeList.overlaps( IndexRange(2,3) ) );
	}

	SECTION("check push_back()") {

		IndexRange r5_10(5,10);

		rangeList.push_back(r5_10);
		REQUIRE( *rangeList.begin() == r5_10 );
		REQUIRE( *rangeList.rbegin() == r5_10 );
		REQUIRE( rangeList.size() == 1 );
		REQUIRE_FALSE( rangeList.empty() );
		REQUIRE( rangeList.begin() != rangeList.end() );
		REQUIRE( rangeList.rbegin() != rangeList.rend() );
		REQUIRE_FALSE( rangeList.covers( 3 ) );
		REQUIRE_FALSE( rangeList.overlaps( IndexRange(2,3) ) );
		REQUIRE( rangeList.covers( 8 ) );
		REQUIRE( rangeList.overlaps( IndexRange(2,8) ) );
		REQUIRE( rangeList.overlaps( IndexRange(6,8) ) );
		REQUIRE( rangeList.overlaps( IndexRange(6,12) ) );

		IndexRange r15_20(15,20);

		rangeList.push_back(r15_20);
		REQUIRE( *rangeList.begin() == r5_10 );
		REQUIRE( *rangeList.begin() != r15_20 );
		REQUIRE( *rangeList.rbegin() == r15_20 );
		REQUIRE( rangeList.size() == 2 );
		REQUIRE_FALSE( rangeList.empty() );
		REQUIRE_FALSE( rangeList.covers( 3 ) );
		REQUIRE_FALSE( rangeList.covers( 13 ) );
		REQUIRE_FALSE( rangeList.overlaps( IndexRange(12,13) ) );
		REQUIRE_FALSE( rangeList.overlaps( IndexRange(28,30) ) );
		REQUIRE( rangeList.covers( 8 ) );
		REQUIRE( rangeList.covers( 18 ) );
		REQUIRE( rangeList.overlaps( IndexRange(2,8) ) );
		REQUIRE( rangeList.overlaps( IndexRange(6,8) ) );
		REQUIRE( rangeList.overlaps( IndexRange(6,12) ) );
		REQUIRE( rangeList.overlaps( IndexRange(6,18) ) );
		REQUIRE( rangeList.overlaps( IndexRange(2,28) ) );
		REQUIRE( rangeList.overlaps( IndexRange(2,28) ) );
		REQUIRE( rangeList.overlaps( IndexRange(12,28) ) );

	}

	SECTION("check insert()") {

		IndexRange r5_10(5,10);

		IndexRangeList::iterator it = rangeList.insert(r5_10);
		REQUIRE( rangeList.begin() == it );
		REQUIRE( *it == r5_10 );
		REQUIRE( *rangeList.begin() == r5_10 );
		REQUIRE( *rangeList.rbegin() == r5_10 );
		REQUIRE( rangeList.size() == 1 );
		REQUIRE_FALSE( rangeList.empty() );
		REQUIRE( rangeList.begin() != rangeList.end() );
		REQUIRE( rangeList.rbegin() != rangeList.rend() );

		IndexRange r1_2(1,2);

		it = rangeList.insert(r1_2);
		REQUIRE( rangeList.begin() == it );
		REQUIRE( *it == r1_2 );
		REQUIRE( *rangeList.begin() == r1_2 );
		REQUIRE( *rangeList.rbegin() == r5_10 );
		REQUIRE( rangeList.size() == 2 );
		REQUIRE_FALSE( rangeList.empty() );

		IndexRange r3_4(3,4);

		it = rangeList.insert(r3_4);
		REQUIRE( rangeList.begin() != it );
		REQUIRE( *it == r3_4 );
		REQUIRE( *rangeList.begin() == r1_2 );
		REQUIRE( *(++rangeList.begin()) == r3_4 );
		REQUIRE( *rangeList.rbegin() == r5_10 );
		REQUIRE( rangeList.size() == 4 );

	}


}

