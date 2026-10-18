#include <IntaRNA/AccessibilityDisabled.h>
#include <IntaRNA/InteractionEnergyBasePair.h>
#include <IntaRNA/OutputHandlerInteractionList.h>
#include <IntaRNA/PredictorMfeEns2dHeuristicSeedExtension.h>
#include <IntaRNA/SeedHandlerNoBulge.h>
INITIALIZE_EASYLOGGINGPP
int main() {
 using namespace IntaRNA;
 RnaSequence t("t","GGGGGGGG"), q("q","CCCC");
 AccessibilityDisabled at(t,4,nullptr), aq(q,4,nullptr);
 ReverseAccessibility ar(aq); InteractionEnergyBasePair e(at,ar);
 OutputConstraint oc(1,OutputConstraint::OVERLAP_BOTH,E_INF,E_INF,false,false,false,true,false);
 OutputHandlerInteractionList out(oc,1);
 SeedConstraint sc(2,0,0,0,E_INF,Accessibility::ED_UPPER_BOUND,E_INF,IndexRangeList(),IndexRangeList(),"",false,false,false);
 PredictorMfeEns2dHeuristicSeedExtension p(e,out,nullptr,new SeedHandlerNoBulge(e,sc));
 p.predict(); return out.empty() || p.getZall()<=0;
}
