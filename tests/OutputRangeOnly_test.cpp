
#include "catch.hpp"

#undef NDEBUG

#include "OutputHandlerRangeOnly.h"
#include <stdexcept>


class RangeStore : public OutputHandler {
public:

	InteractionRange lastAdded;

	RangeStore( const RnaSequence& s ) : lastAdded(s,s) {}

	void add( const Interaction& i ) {
		// reset
		lastAdded = InteractionRange(*(i.s1),*(i.s2));
	}
	void add( const InteractionRange& r ) {
		// store
		lastAdded = r;
	}
};

TEST_CASE( "OutputHandlerRangeOnly", "[OutputHandlerRangeOnly]" ) {

	RnaSequence r("test","AACCGGUU");

	RangeStore succOut(r);

	OutputHandlerRangeOnly out(succOut);

	Interaction i(r,r);
	i.addInteraction(0,7);
	i.addInteraction(1,6);
	REQUIRE( i.isValid() );

	InteractionRange ir(i);

	SECTION("add interaction") {
		out.add(i);
		REQUIRE( (!(succOut.lastAdded < ir)  &&  !(ir < succOut.lastAdded)));
	}

	SECTION("add interaction range") {
		out.add(ir);
		REQUIRE( (!(succOut.lastAdded < ir)  &&  !(ir < succOut.lastAdded)));
	}

}
