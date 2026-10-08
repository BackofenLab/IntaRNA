#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/OutputHandlerHub.h"
#include "IntaRNA/RnaSequence.h"

using namespace IntaRNA;

namespace {

class CountingOutputHandler : public OutputHandler {
public:
	explicit CountingOutputHandler(const OutputConstraint & constraint)
	 : OutputHandler(constraint)
	{}

	void add(const Interaction &) override
	{
		++reportedInteractions;
	}
};

} // namespace

TEST_CASE("OutputHandlerHub reports the largest child count", "[OutputHandlerHub]") {

	#include "testEasyLoggingSetup.icc"

	OutputConstraint constraint;
	CountingOutputHandler first(constraint);
	CountingOutputHandler second(constraint);
	OutputHandlerHub hub(constraint, false);

	REQUIRE(hub.reported() == 0);
	hub.addOutputHandler(&first);
	hub.addOutputHandler(&second);

	RnaSequence sequence("sequence", "GG");
	Interaction interaction(sequence, sequence);
	hub.add(interaction);

	REQUIRE(first.reported() == 1);
	REQUIRE(second.reported() == 1);
	REQUIRE(hub.reported() == 1);

	second.add(interaction);
	REQUIRE(first.reported() == 1);
	REQUIRE(second.reported() == 2);
	REQUIRE(hub.reported() == 2);
}

TEST_CASE("exact-zero partition status follows copied hubs and new children", "[OutputHandlerHub]") {
	OutputConstraint constraint;
	CountingOutputHandler output(constraint),other(constraint);
	OutputHandlerHub hub(constraint,false);hub.addOutputHandler(&output);hub.setExactPartition(true);
	hub.incrementZ(1e-100);
	REQUIRE(hub.hasNonzeroZ());
	OutputHandlerHub copy(hub);REQUIRE(copy.usesExactPartition());REQUIRE(copy.hasNonzeroZ());
	copy.addOutputHandler(&other);REQUIRE(other.usesExactPartition());
	OutputHandlerHub assigned(constraint,false);assigned=hub;REQUIRE(assigned.usesExactPartition());
}
