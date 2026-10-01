#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/InteractionEnergyVrna.h"
#include "IntaRNA/OutputHandlerInteractionList.h"
#include "IntaRNA/PredictorEvalOnly.h"
#include "IntaRNA/PredictorMfe2d.h"

#include <cmath>
#include <limits>
#include <tuple>
#include <vector>

using namespace IntaRNA;

namespace {

class EvalAccessibility : public Accessibility {
public:
	EvalAccessibility(const RnaSequence & sequence, const E_type perBase);
	E_type getED(const size_t from, const size_t to) const override;
private:
	E_type perBase;
};

inline EvalAccessibility::EvalAccessibility(const RnaSequence & sequence, const E_type perBase)
	: Accessibility(sequence, 0, NULL), perBase(perBase)
{}

inline E_type EvalAccessibility::getED(const size_t from, const size_t to) const
{
	checkIndices(from, to);
	return (to-from+1)*perBase;
}

using EvalUpdate = std::tuple<size_t, size_t, size_t, size_t, E_type>;

class EvalTracker : public PredictionTracker {
public:
	EvalTracker(std::vector<EvalUpdate> & updates);
	void updateOptimumCalled(const size_t i1, const size_t j1,
			const size_t i2, const size_t j2, const E_type energy) override;
private:
	std::vector<EvalUpdate> & updates;
};

inline EvalTracker::EvalTracker(std::vector<EvalUpdate> & updates) : updates(updates) {}

inline void EvalTracker::updateOptimumCalled(const size_t i1, const size_t j1,
		const size_t i2, const size_t j2, const E_type energy)
{
	updates.emplace_back(i1, j1, i2, j2, energy);
}

} // namespace

TEST_CASE("Evaluation parses hybridDB in original sequence coordinates", "[PredictorEvalOnly]")
{
	#include "testEasyLoggingSetup.icc"
	RnaSequence target("t", "GGAGG", -2), query("q", "CCACC", 8);
	auto input = PredictorEvalOnly::parseInteractions("-2||.||&8||.||:-1|&9|", target, query);
	REQUIRE(input.size() == 2);
	const Interaction::PairingVec expected = {{0,4}, {1,3}, {3,1}, {4,0}};
	REQUIRE(input[0].basePairs == expected);
	REQUIRE(Interaction::dotBar(input[0]) == "-2||.||&8||.||");
	REQUIRE(input[1].basePairs.front() == Interaction::BasePair(1,1));
	REQUIRE_THROWS_AS(PredictorEvalOnly::parseInteractions("0|&8|", target, query), std::invalid_argument);

	RnaSequence g("g", "GGGG"), c("c", "CCCC");
	auto full = PredictorEvalOnly::parseInteractions("1.||.&1.||.", g, c);
	REQUIRE(Interaction::dotBar(full[0], true) == "1.||.&1.||.");
	for (const auto & bad : {"", "1|||", "1|&1|&1|", "1|&1||", "1...&1...",
			"1x|&1||", "1(&1)", "0|&1|", "-1|&1|", "4||&1||", "1|&5|",
			"9223372036854775808|&1|", "-9223372036854775809|&1|", "1|&1|:", ":1|&1|", "1|&1|::1|&1|"}) {
		CAPTURE(bad);
		REQUIRE_THROWS_AS(PredictorEvalOnly::parseInteractions(bad, g, c), std::invalid_argument);
	}
	REQUIRE_THROWS_AS(PredictorEvalOnly::parseInteractions("1|&1|", g, g), std::invalid_argument);
	RnaSequence n("n", "NCCC");
	REQUIRE_THROWS_AS(PredictorEvalOnly::parseInteractions("1|&1|", g, n), std::invalid_argument);
	RnaSequence zero("zero", "GGGG", 0);
	REQUIRE(PredictorEvalOnly::parseInteractions("0|&1|", zero, c)[0].basePairs.front().first == 0);
}

TEST_CASE("Evaluation reports energies, structures and restricted ensemble without filtering", "[PredictorEvalOnly]")
{
	#include "testEasyLoggingSetup.icc"
	RnaSequence target("t", "GGGGGG"), query("q", "UUUUUUU");
	EvalAccessibility acc1(target, 20), acc2(query, 30);
	ReverseAccessibility reversed(acc2);
	InteractionEnergyBasePair energy(acc1, reversed, 20, 20, false, 1, -100, 3, 500);
	// Deliberately exclude these positive-energy, lonely GU pairs with filters.
	OutputConstraint constraints(0, OutputConstraint::OVERLAP_NONE, -1000, 0, true, true, true, true, false, 0);
	OutputHandlerInteractionList output(constraints, 10);
	std::vector<EvalUpdate> updates;
	auto input = PredictorEvalOnly::parseInteractions("2||.||&2|..|||:1|&1|:2||.||&2|..|||", target, query);
	input[0].energy = -12345;
	input[0].setSeedRange(input[0].basePairs.front(), input[0].basePairs.back(), -999);
	PredictorEvalOnly predictor(energy, output, new EvalTracker(updates), input);
	// Ranges do not trim or shift the supplied coordinates.
	predictor.predict(IndexRange(0,0), IndexRange(0,0));
	REQUIRE(output.reported() == 2);
	auto it = output.begin();
	REQUIRE((*it)->energy == 380); // -4 bp + 5*0.2 ED1 + 6*0.3 ED2 + 5 shift
	REQUIRE((*it)->basePairs == input[0].basePairs);
	REQUIRE((*it)->seed == NULL);
	REQUIRE((*++it)->energy == 450);
	REQUIRE(updates.size() == 2);
	REQUIRE(updates[0] == EvalUpdate(1,5,0,5,380));
	const double partition = std::exp(-3.8) + std::exp(-4.5);
	REQUIRE(static_cast<double>(predictor.getZall()) == Approx(partition));
	REQUIRE(static_cast<double>(output.getZ()) == Approx(partition));
	predictor.predict();
	REQUIRE(static_cast<double>(predictor.getZall()) == Approx(partition));
	REQUIRE(static_cast<double>(output.getZ()) == Approx(2*partition));
	REQUIRE(updates.size() == 4);
}

TEST_CASE("Evaluation validates API inputs and reports unavailable energies before output", "[PredictorEvalOnly]")
{
	#include "testEasyLoggingSetup.icc"
	RnaSequence target("t", "GGGG"), query("q", "CCCC");
	AccessibilityDisabled acc1(target, 0, NULL), acc2(query, 0, NULL);
	ReverseAccessibility reversed(acc2);
	InteractionEnergyBasePair energy(acc1, reversed, 0, 0);
	OutputHandlerInteractionList output(OutputConstraint(), 10);
	REQUIRE_THROWS_AS(PredictorEvalOnly(energy, output, NULL, {}), std::invalid_argument);
	Interaction invalid(target, query);
	REQUIRE_THROWS_AS(PredictorEvalOnly(energy, output, NULL, {invalid}), std::invalid_argument);
	invalid.basePairs = {{0,0}, {1,1}};
	REQUIRE_THROWS_AS(PredictorEvalOnly(energy, output, NULL, {invalid}), std::invalid_argument);
	invalid.basePairs = {{4,0}};
	REQUIRE_THROWS_AS(PredictorEvalOnly(energy, output, NULL, {invalid}), std::invalid_argument);
	Interaction mismatch(query, target);
	mismatch.basePairs = {{0,0}};
	REQUIRE_THROWS_AS(PredictorEvalOnly(energy, output, NULL, {mismatch}), std::invalid_argument);
	auto input = PredictorEvalOnly::parseInteractions("1|&1|:1|..|&1|..|", target, query);
	PredictorEvalOnly predictor(energy, output, NULL, input);
	REQUIRE_THROWS_AS(predictor.predict(), std::runtime_error);
	REQUIRE(output.empty());
	REQUIRE(output.getZ() == 0);

	AccessibilityDisabled narrow(target, 1, NULL);
	InteractionEnergyBasePair limited(narrow, reversed, 10, 10);
	PredictorEvalOnly inaccessible(limited, output, NULL, input);
	REQUIRE_THROWS_AS(inaccessible.predict(), std::runtime_error);
	REQUIRE(output.empty());
}

TEST_CASE("Evaluation reproduces ViennaRNA prediction energies and contributions", "[PredictorEvalOnly]")
{
	#include "testEasyLoggingSetup.icc"
	RnaSequence target("t", "AGCGACGCA"), query("q", "UGCGUCGCU");
	EvalAccessibility acc1(target, 7), acc2(query, 13);
	ReverseAccessibility reversed(acc2);
	VrnaHandler vrna;
	for (const bool dangles : {false, true}) {
		InteractionEnergyVrna energy(acc1, reversed, vrna, 20, 20, false, 37, dangles);
		OutputConstraint constraints(20, OutputConstraint::OVERLAP_BOTH, E_INF, E_INF);
		OutputHandlerInteractionList predicted(constraints, 20), evaluated(constraints, 20);
		PredictorMfe2d search(energy, predicted, NULL);
		search.predict();
		REQUIRE_FALSE(predicted.empty());
		std::vector<Interaction> input;
		for (const auto * interaction : predicted) input.push_back(*interaction);
		PredictorEvalOnly evaluation(energy, evaluated, NULL, input);
		evaluation.predict();
		REQUIRE(evaluated.reported() == predicted.reported());
		auto result = evaluated.begin();
		for (const auto * reference : predicted) {
			REQUIRE((*result)->basePairs == reference->basePairs);
			REQUIRE((*result)->energy == reference->energy);
			const auto parts = energy.getE_contributions(**result);
			REQUIRE(parts.init + parts.loops + parts.ED1 + parts.ED2 + parts.dangleLeft
					+ parts.dangleRight + parts.endLeft + parts.endRight + parts.energyAdd == reference->energy);
			++result;
		}
	}
}
