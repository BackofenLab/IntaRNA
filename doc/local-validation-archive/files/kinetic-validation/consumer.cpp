#include <IntaRNA/PredictorSeedExtensionKinetic.h>
#include <IntaRNA/AccessibilityDisabled.h>
#include <IntaRNA/InteractionEnergyBasePair.h>
#include <IntaRNA/OutputHandlerInteractionList.h>
#include <IntaRNA/SeedHandlerNoBulge.h>

INITIALIZE_EASYLOGGINGPP

int main() {
    using namespace IntaRNA;
    RnaSequence target("target", "GGGG"), query("query", "CCCC");
    AccessibilityDisabled a(target, 0, nullptr), b(query, 0, nullptr);
    ReverseAccessibility reverse(b);
    InteractionEnergyBasePair energy(a, reverse, 1, 1);
    SeedConstraint seeds(2, 0, 0, 0, E_INF, Accessibility::ED_UPPER_BOUND,
            E_INF, IndexRangeList(""), IndexRangeList(""), "", false, false, true);
    OutputHandlerInteractionList output(OutputConstraint(), 1);
    PredictorSeedExtensionKinetic predictor(energy, output, nullptr,
            new SeedHandlerNoBulge(energy, seeds));
    predictor.predict();
    return output.empty() || (*output.begin())->energy != -400;
}
