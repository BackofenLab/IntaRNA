#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/OutputHandlerCsv.h"
#include "IntaRNA/OutputHandlerInteractionList.h"
#include "IntaRNA/OutputHandlerMfe.h"
#include "IntaRNA/SeedHandlerExplicit.h"

#include <sstream>

using namespace IntaRNA;

TEST_CASE("OutputHandlerMfe", "[OutputHandlerMfe]") {

	// setup easylogging++ stuff if not already done
	#include "testEasyLoggingSetup.icc"

	RnaSequence target("target","GGGG");
	RnaSequence query("query","CCCC");
	AccessibilityDisabled targetAcc(target,target.size(),NULL);
	AccessibilityDisabled queryAcc(query,query.size(),NULL);
	ReverseAccessibility queryAccReversed(queryAcc);
	InteractionEnergyBasePair energy(targetAcc,queryAccReversed,16,16,false,
			1,Ekcal_2_E(-1.0),1);
	OutputConstraint constraint(10,OutputConstraint::OVERLAP_BOTH,
			Ekcal_2_E(100.0),E_INF,false,false,false,false,true);

	// Keep a raw observer while the adapter owns the final collector.
	OutputHandlerInteractionList * stored =
			new OutputHandlerInteractionList(constraint,10);
	OutputHandlerMfe output(stored,energy);
	Interaction site(target,query);
	// Include predictor-traced interior pairs and seed metadata. The adapter must
	// normalize these back to boundaries before its representative traceback.
	for (size_t i=0; i<target.size(); ++i)
		site.basePairs.push_back(Interaction::BasePair(i,query.size()-i-1));
	site.setSeedRange(site.basePairs.at(1),site.basePairs.at(2),Ekcal_2_E(-2.0));
	site.energy = Ekcal_2_E(-5.25);

	output.add(site);
	REQUIRE(output.reported() == 1);
	REQUIRE(stored->begin() != stored->end());
	const Interaction & representative = **stored->begin();
	REQUIRE(representative.basePairs.size() == 4);
	REQUIRE(representative.seed == NULL);
	REQUIRE(E_equal(representative.energy,Ekcal_2_E(-4.0)));
	REQUIRE(E_equal(representative.ensembleEnergy,site.energy));
	REQUIRE(E_equal(representative.getRankingEnergy(),site.energy));

	std::ostringstream csv;
	OutputHandlerCsv csvOutput(constraint,csv,energy,
			OutputHandlerCsv::string2list("E,Eens"),";",true);
	csvOutput.add(representative);
	REQUIRE(csv.str() == "E;Eens\n-4;-5.25\n");

	SeedConstraint explicitConstraint(2,0,0,0,0,
			AccessibilityDisabled::ED_UPPER_BOUND,0,IndexRangeList(""),
			IndexRangeList(""),"2|&2|",false,false,false);
	OutputHandlerInteractionList * explicitStored =
			new OutputHandlerInteractionList(constraint,1);
	OutputHandlerMfe explicitOutput(explicitStored,energy,
			new SeedHandlerExplicit(energy,explicitConstraint));
	Interaction singleton(target,query);
	singleton.basePairs.push_back(Interaction::BasePair(1,1));
	singleton.energy = Ekcal_2_E(-1.5);
	explicitOutput.add(singleton);
	REQUIRE(explicitStored->begin() != explicitStored->end());
	const Interaction & explicitRepresentative = **explicitStored->begin();
	REQUIRE(explicitRepresentative.basePairs.size() == 1);
	REQUIRE(explicitRepresentative.seed != NULL);
	REQUIRE(E_equal(explicitRepresentative.energy,Ekcal_2_E(-1.0)));
	REQUIRE(E_equal(explicitRepresentative.ensembleEnergy,singleton.energy));
}
