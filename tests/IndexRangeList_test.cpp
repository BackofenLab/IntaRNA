


#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/IndexRangeList.h"

using namespace IntaRNA;


TEST_CASE( "IndexRangeList", "[IndexRangeList]" ) {

	// setup easylogging++ stuff if not already done
	#include "testEasyLoggingSetup.icc"

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
		REQUIRE( rangeList.covers( 8, 10 ) );
		REQUIRE( rangeList.covers( 5, 10 ) );
		REQUIRE( rangeList.covers( 5, 6 ) );
		REQUIRE_FALSE( rangeList.covers( 5, 12 ) );
		REQUIRE_FALSE( rangeList.covers( 1, 12 ) );
		REQUIRE_FALSE( rangeList.covers( 1, 8 ) );
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
		REQUIRE( rangeList.covers( 8 ) );
		REQUIRE( rangeList.covers( 8, 10 ) );
		REQUIRE( rangeList.covers( 5, 10 ) );
		REQUIRE( rangeList.covers( 5, 6 ) );
		REQUIRE( rangeList.covers( 16, 18 ) );
		REQUIRE_FALSE( rangeList.covers( 12, 13 ) );
		REQUIRE_FALSE( rangeList.covers( 12, 16 ) );
		REQUIRE_FALSE( rangeList.covers( 10, 18 ) );
		REQUIRE( rangeList.overlaps( IndexRange(2,8) ) );
		REQUIRE( rangeList.overlaps( IndexRange(5,8) ) );
		REQUIRE( rangeList.overlaps( IndexRange(6,8) ) );
		REQUIRE( rangeList.overlaps( IndexRange(6,10) ) );
		REQUIRE( rangeList.overlaps( IndexRange(6,12) ) );
		REQUIRE( rangeList.overlaps( IndexRange(5,18) ) );
		REQUIRE( rangeList.overlaps( IndexRange(6,18) ) );
		REQUIRE( rangeList.overlaps( IndexRange(2,28) ) );
		REQUIRE( rangeList.overlaps( IndexRange(2,28) ) );
		REQUIRE( rangeList.overlaps( IndexRange(12,28) ) );

	}

	SECTION("check insert()") {

		rangeList.clear();

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
		REQUIRE( rangeList.size() == 3 );

	}


	SECTION("check string en-/decoding") {

		// write to string
		std::stringstream s;
		// parse from string
		IndexRangeList list2(s.str());
		REQUIRE( rangeList == list2 );

		rangeList.insert(IndexRange(4,8));
		s.str(""); s <<rangeList;
		list2=IndexRangeList(s.str());
		REQUIRE( rangeList == list2 );

		rangeList.insert(IndexRange(1,2));
		s.str(""); s <<rangeList;
		list2=IndexRangeList(s.str());
		REQUIRE( rangeList == list2 );

		rangeList.insert(IndexRange(10,10));
		s.str(""); s <<rangeList;
		list2=IndexRangeList(s.str());
		REQUIRE( rangeList == list2 );

	}

	SECTION("check shift") {

		// create list
		rangeList.insert(IndexRange(1,2));
		rangeList.insert(IndexRange(4,8));
		rangeList.insert(IndexRange(10,10));

		// shift to left
		IndexRangeList l5 = rangeList.shift(-5,10);
		IndexRangeList r1 = rangeList.shift(+1,10);

		// write to string for comparison
		std::stringstream s;
		s.str(""); s <<rangeList;
		REQUIRE( s.str() == "1-2,4-8,10-10" );
		s.str(""); s <<l5;
		REQUIRE( s.str() == "0-3,5-5" );
		s.str(""); s<<r1;
		REQUIRE( s.str() == "2-3,5-9" );

	}

	SECTION("check get") {

		// create list
		rangeList.insert(IndexRange(1,2));
		rangeList.insert(IndexRange(4,8));
		rangeList.insert(IndexRange(10,10));

		REQUIRE( rangeList.get(0) == IndexRange(1,2) );
		REQUIRE( rangeList.get(1) == IndexRange(4,8) );
		REQUIRE( rangeList.get(2) == IndexRange(10,10) );

	}

	SECTION("check reverse") {

		// create list
		rangeList.insert(IndexRange(1,2));
		rangeList.insert(IndexRange(4,8));
		rangeList.insert(IndexRange(10,10));

		rangeList.reverseInplace(11);
		REQUIRE( toString(rangeList) == "0-0,2-6,8-9" );

	}

}

