// Standalone installed-library example; see README.md#bpProb.
#include <IntaRNA/AccessibilityDisabled.h>
#include <IntaRNA/BasePairProbabilityWriter.h>
#include <IntaRNA/InteractionEnergyBasePair.h>
#include <IntaRNA/OutputHandlerInteractionList.h>
#include <IntaRNA/PredictorMfeEns2dSeedExtension.h>
#include <IntaRNA/SeedHandlerNoBulge.h>
#include <iostream>

INITIALIZE_EASYLOGGINGPP

int main() {
	using namespace IntaRNA;
	RnaSequence target("target", "GGAGGG"), query("query", "CCCC");
	AccessibilityDisabled targetAcc(target, 0, nullptr), queryAcc(query, 0, nullptr);
	ReverseAccessibility reversed(queryAcc);
	InteractionEnergyBasePair energy(targetAcc, reversed);
	SeedConstraint seeds(2, 0, 0, 0, E_INF, Accessibility::ED_UPPER_BOUND,
			E_INF, IndexRangeList(), IndexRangeList(), "", false, false, false);
	OutputConstraint constraints(1, OutputConstraint::OVERLAP_BOTH, E_INF, E_INF,
			false, false, false, true, false); // needZall=true; needBPs=false
	OutputHandlerInteractionList output(constraints, 1);
	BasePairProbabilities probabilities(target.size(), query.size());
	PredictorMfeEns2dSeedExtension predictor(energy, output, nullptr,
			new SeedHandlerNoBulge(energy, seeds), &probabilities);
	predictor.predict();
	probabilities.finalize();
	BasePairProbabilityWriter::write(std::cout, probabilities, target, query);
}
