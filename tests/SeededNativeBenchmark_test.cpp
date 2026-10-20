#include "catch.hpp"
#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/AccessibilityVrna.h"
#include "IntaRNA/InteractionEnergyVrna.h"
#include "IntaRNA/SeedHandlerNoBulge.h"
#include "IntaRNA/PredictorMfeEns2dSeedExtension.h"
#include "IntaRNA/OutputHandlerInteractionList.h"
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <memory>

using namespace IntaRNA;
namespace {
using Clock=std::chrono::steady_clock;
double elapsed(Clock::time_point start) { return std::chrono::duration<double,std::milli>(Clock::now()-start).count(); }
class TimedStackSeeds:public SeedHandlerNoBulge {
	double & ms;
public:
	TimedStackSeeds(const InteractionEnergy & e,const SeedConstraint & c,double & ms):SeedHandlerNoBulge(e,c),ms(ms) {}
	size_t fillSeed(size_t i,size_t j,size_t k,size_t l) override {
		auto start=Clock::now();auto count=SeedHandlerNoBulge::fillSeed(i,j,k,l);ms=elapsed(start);return count;
	}
};
}
TEST_CASE("native seeded preprocessing partition and outside benchmark", "[.][SeededNativeBenchmark]") {
	#include "testEasyLoggingSetup.icc"
	for(std::string name:{"none","sparse100","sparse400","dense","OxyS-fhlA"}) {
		const char * requested=std::getenv("INTARNA_BENCH_CASE");
		if(requested && name!=requested) continue;
		const size_t n=name=="sparse400"?400:100;
		std::string target(n,'A'),query(n,'A');
		if(name=="sparse100" || name=="sparse400") { target.replace(n/2,8,8,'G');query.replace(n/2,8,8,'C'); }
		if(name=="dense") { target=std::string(40,'G');query=std::string(40,'C'); }
		if(name=="OxyS-fhlA") { target="AGUUAGUCAAUGACCUUUUGCACCGCUUUGCGGUGCUUUCCUGGAACAACAAAAUGUCAUAUACACCGAUGAGUGAUCUCGGACAACAAGGGUUGUUCGACAUCACUCGGAC";query="GAAACGGAGCGGCACCUCUUUUAACCCUUGAAGUCACUGCCCGUUUCGAGAGUUUCUCAACUCGAAUAACUAAAGCCAACGUGAACUUUUGCGGAUCUCCAGGAUCCG"; }
		RnaSequence t("target",target),q("query",query);
		VrnaHandler vrna;
		auto start=Clock::now();
		std::unique_ptr<Accessibility> at,aq;
		if(name=="OxyS-fhlA") {
			at=std::make_unique<AccessibilityVrna>(t,30,nullptr,vrna,0);
			aq=std::make_unique<AccessibilityVrna>(q,30,nullptr,vrna,0);
		} else {
			at=std::make_unique<AccessibilityDisabled>(t,30,nullptr);
			aq=std::make_unique<AccessibilityDisabled>(q,30,nullptr);
		}
		const double accMs=elapsed(start);
		ReverseAccessibility ar(*aq);
		InteractionEnergyVrna energy(*at,ar,vrna,3,3);
		SeedConstraint sc(4,0,0,0,0,Accessibility::ED_UPPER_BOUND,E_INF,IndexRangeList(),IndexRangeList(),"",false,false,false);
		Z_type reference=0;
		for(bool pairs:{false,true}) {
			const char * mode=std::getenv("INTARNA_BENCH_OUTSIDE");
			if(mode && pairs!=(std::string(mode)=="1")) continue;
			OutputConstraint oc(0,OutputConstraint::OVERLAP_BOTH,E_INF,E_INF,false,false,false,true,false);
			OutputHandlerInteractionList output(oc,0);
			std::unique_ptr<BasePairProbabilities> result;
			start=Clock::now();
			if(pairs) result=std::make_unique<BasePairProbabilities>(t.size(),q.size());
			double seedMs=0;
			PredictorMfeEns2dSeedExtension p(energy,output,nullptr,new TimedStackSeeds(energy,sc,seedMs),result.get());
			p.predict();if(result) result->finalize();
			const double ms=elapsed(start);
			if(!pairs) reference=p.getZall();else if(!mode) REQUIRE(p.getZall()==reference);
			std::cout<<name<<","<<pairs<<","<<accMs<<","<<seedMs<<","<<ms-seedMs<<","<<ms+accMs<<","<<p.getZall()<<"\n";
		}
	}
}
