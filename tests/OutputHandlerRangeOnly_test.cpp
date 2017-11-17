
#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/OutputHandlerRangeOnly.h"
#include <stdexcept>

using namespace IntaRNA;

/**
 * Dummy class for testing the add result
 */
class RangeStore : public OutputHandler {
public:

	InteractionRange lastAdded;

	RangeStore( const RnaSequence& s ) : lastAdded(s,s) {}

	void add( const Interaction& i ) {
		// reset, ignoring the given interaction
		lastAdded = InteractionRange(*(i.s1),*(i.s2));
	}
	void add( const InteractionRange& r ) {
		// store the range
		lastAdded = r;
	}
};

TEST_CASE( "OutputHandlerRangeOnly", "[OutputHandlerRangeOnly]" ) {

	RnaSequence r("test","AACCGGUU");

	RangeStore succOut(r);

	OutputHandlerRangeOnly out(succOut);

	Interaction i(r,r);
	i.basePairs.push_back( Interaction::BasePair(0,7));
	i.basePairs.push_back( Interaction::BasePair(1,6));
	REQUIRE( i.isValid() );

	InteractionRange ir(i);

	SECTION("add interaction") {
		// ensure difference before add
		REQUIRE( ((succOut.lastAdded < ir)  ||  (ir < succOut.lastAdded)));
		out.add(i);
		// check equivalence after add
		REQUIRE( (!(succOut.lastAdded < ir)  &&  !(ir < succOut.lastAdded)));
	}

	SECTION("add interaction range") {
		// ensure difference before add
		REQUIRE( ((succOut.lastAdded < ir)  ||  (ir < succOut.lastAdded)));
		out.add(ir);
		// check equivalence after add
		REQUIRE( (!(succOut.lastAdded < ir)  &&  !(ir < succOut.lastAdded)));
	}

}
