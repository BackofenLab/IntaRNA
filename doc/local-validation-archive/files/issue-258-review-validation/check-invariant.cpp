#include <IntaRNA/StackSeedDomain.h>
#include <IntaRNA/SeedHandlerNoBulge.h>
#include <IntaRNA/AccessibilityDisabled.h>
#include <IntaRNA/InteractionEnergyBasePair.h>
#include <iostream>
INITIALIZE_EASYLOGGINGPP

class ObservedSeeds : public IntaRNA::SeedHandlerNoBulge {
public:
    using SeedHandlerNoBulge::SeedHandlerNoBulge;
    mutable size_t secondLengthCalls=0;
    bool mismatch=false;
    size_t getSeedLength2(size_t i,size_t j) const override {
        ++secondLengthCalls;
        return SeedHandlerNoBulge::getSeedLength2(i,j)+(mismatch?1:0);
    }
};

int main(int argc,char **) {
    using namespace IntaRNA;
    RnaSequence target("target","GGG"),query("query","CCC");
    AccessibilityDisabled at(target,0,nullptr),aq(query,0,nullptr);
    ReverseAccessibility ar(aq);
    InteractionEnergyBasePair energy(at,ar);
    SeedConstraint constraint(2,0,0,0,E_INF,Accessibility::ED_UPPER_BOUND,E_INF,
        IndexRangeList(),IndexRangeList(),"",false,false,false);
    ObservedSeeds handler(energy,constraint);
    handler.fillSeed(0,2,0,2);
    handler.secondLengthCalls=0;
    handler.mismatch=argc>1;
    StackSeedDomain domain(handler,3,3,3,3);
    std::cout << "length=" << domain.maxSeedLength()
              << " secondLengthCalls=" << handler.secondLengthCalls << '\n';
    return domain.maxSeedLength()==2?0:1;
}
