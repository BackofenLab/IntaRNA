#include <IntaRNA/PredictorEvalOnly.h>
#include <IntaRNA/AccessibilityDisabled.h>
#include <IntaRNA/InteractionEnergyBasePair.h>
#include <IntaRNA/OutputHandlerInteractionList.h>
INITIALIZE_EASYLOGGINGPP
int main() {
    using namespace IntaRNA;
    RnaSequence t("t", "GG"), q("q", "CC");
    AccessibilityDisabled a(t, 0, nullptr), b(q, 0, nullptr);
    ReverseAccessibility reversed(b);
    InteractionEnergyBasePair energy(a, reversed);
    OutputHandlerInteractionList output(OutputConstraint(), 10);
    PredictorEvalOnly eval(energy, output, nullptr,
        PredictorEvalOnly::parseInteractions("1||&1||", t, q));
    eval.predict();
    return output.reported() == 1 && (*output.begin())->energy == -200 ? 0 : 1;
}
