#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/OutputHandlerInteractionList.h"
#include "IntaRNA/PredictorMfeEns2dHeuristic.h"
#include "IntaRNA/ReverseAccessibility.h"
#include "IntaRNA/RnaSequence.h"
INITIALIZE_EASYLOGGINGPP
using namespace IntaRNA;
int main() {
 RnaSequence target("target", "CC"), query("query", "GG");
 AccessibilityDisabled targetAcc(target,0,nullptr), queryAcc(query,0,nullptr);
 ReverseAccessibility reverse(queryAcc);
 InteractionEnergyBasePair energy(targetAcc,reverse);
 OutputConstraint constraint(2,OutputConstraint::OVERLAP_NONE,0,10000);
 OutputHandlerInteractionList output(constraint,2);
 PredictorMfeEns2dHeuristic predictor(energy,output,nullptr);
 predictor.predict();
}
