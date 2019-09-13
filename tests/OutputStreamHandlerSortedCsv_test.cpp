
#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/OutputStreamHandlerSortedCsv.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/AccessibilityDisabled.h"

#include <stdexcept>

using namespace IntaRNA;

TEST_CASE( "OutputStreamHandlerSortedCsv", "[OutputStreamHandlerSortedCsv]" ) {

	// setup easylogging++ stuff if not already done
	#include "testEasyLoggingSetup.icc"

	const std::string csvUnsorted =
			"id1,start1,E,id2\n"
			"eee,45|000,+5,cc\n"
			"bbb,00|105,-5,ee\n"
			"aaa,7|0005,20,bb\n"
			"ege,005|23,00,aa\n"
			;
	const std::string csvSorted0 =
			"id1,start1,E,id2\n"
			"aaa,7|0005,20,bb\n"
			"bbb,00|105,-5,ee\n"
			"eee,45|000,+5,cc\n"
			"ege,005|23,00,aa\n"
			;
	const std::string csvSorted1 =
			"id1,start1,E,id2\n"
			"bbb,00|105,-5,ee\n"
			"ege,005|23,00,aa\n"
			"aaa,7|0005,20,bb\n"
			"eee,45|000,+5,cc\n"
			;
	const std::string csvSorted2 =
			"id1,start1,E,id2\n"
			"bbb,00|105,-5,ee\n"
			"ege,005|23,00,aa\n"
			"eee,45|000,+5,cc\n"
			"aaa,7|0005,20,bb\n"
			;
	const std::string csvSorted3 =
			"id1,start1,E,id2\n"
			"ege,005|23,00,aa\n"
			"aaa,7|0005,20,bb\n"
			"eee,45|000,+5,cc\n"
			"bbb,00|105,-5,ee\n"
			;

	SECTION("sort id1") {
		std::stringstream outStream;
		{
			OutputStreamHandlerSortedCsv oshSorted( new OutputStreamHandler(&outStream), 0, true, ",", true);
			oshSorted.getOutStream() << csvUnsorted;
		} // ensure oshSorted is destroyed
		REQUIRE( outStream.str() == csvSorted0 );
	}

	SECTION("sort start1 list") {
		std::stringstream outStream;
		{
			OutputStreamHandlerSortedCsv oshSorted( new OutputStreamHandler(&outStream), 1, false, ",", true,"|");
			oshSorted.getOutStream() << csvUnsorted;
		} // ensure oshSorted is destroyed
		REQUIRE( outStream.str() == csvSorted1 );
	}

	SECTION("sort E") {
		std::stringstream outStream;
		{
			OutputStreamHandlerSortedCsv oshSorted( new OutputStreamHandler(&outStream), 2, false, ",", true);
			oshSorted.getOutStream() << csvUnsorted;
		} // ensure oshSorted is destroyed
		REQUIRE( outStream.str() == csvSorted2 );
	}

	SECTION("sort E pseudo list") {
		std::stringstream outStream;
		{
			OutputStreamHandlerSortedCsv oshSorted( new OutputStreamHandler(&outStream), 2, false, ",", true, "|");
			oshSorted.getOutStream() << csvUnsorted;
		} // ensure oshSorted is destroyed
		REQUIRE( outStream.str() == csvSorted2 );
	}

	SECTION("sort id2") {
		std::stringstream outStream;
		{
			OutputStreamHandlerSortedCsv oshSorted( new OutputStreamHandler(&outStream), 3, true, ",", true);
			oshSorted.getOutStream() << csvUnsorted;
		} // ensure oshSorted is destroyed
		REQUIRE( outStream.str() == csvSorted3 );
	}

}
