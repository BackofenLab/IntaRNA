
#include "InteractionEnergyVrna.h"

#include <cassert>
#include <set>



////////////////////////////////////////////////////////////////////////////

InteractionEnergyVrna::InteractionEnergyVrna(
		const Accessibility & accS1
		, const ReverseAccessibility & accS2
		, VrnaHandler &vrnaHandler
		, const size_t maxInternalLoopSize1
		, const size_t maxInternalLoopSize2
	)
 :
	InteractionEnergy(accS1, accS2, maxInternalLoopSize1, maxInternalLoopSize2)
// get final VRNA folding parameters
	, foldModel( vrnaHandler.getModel() )
	, foldParams( vrna_params( &foldModel ) )
	, RT(vrnaHandler.getRT())
	, bpCG( BP_pair[RnaSequence::getCodeForChar('C')][RnaSequence::getCodeForChar('G')] )
	, bpGC( BP_pair[RnaSequence::getCodeForChar('G')][RnaSequence::getCodeForChar('C')] )
{
	vrna_md_defaults_reset( &foldModel );
}

////////////////////////////////////////////////////////////////////////////

InteractionEnergyVrna::~InteractionEnergyVrna()
{
	// garbage collection
	if (foldParams != NULL) {
		free(foldParams);
		foldParams = NULL;
	}
}


////////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergyVrna::
getBestE_interLoop() const
{
	// TODO maybe setup member variable (init=E_INF) with lazy initialization

	// get all possible base pair codes handled
	std::set<int> basePairCodes;
	for (int i=0; i<NBASES; i++) {
	for (int j=0; j<NBASES; j++) {
		basePairCodes.insert( BP_pair[i][j] );
	}
	}
	// get minimal energy for any base pair code combination
	E_type minStackingE = E_INF;

	for (std::set<int>::const_iterator p1=basePairCodes.begin(); p1!=basePairCodes.end(); p1++) {
	for (std::set<int>::const_iterator p2=basePairCodes.begin(); p2!=basePairCodes.end(); p2++) {
		minStackingE = std::min( minStackingE
				, (E_type)E_IntLoop(	0	// unpaired region 1
						, 0	// unpaired region 2
						, *p1	// type BP (i1,i2)
						, *p2	// type BP (j2,j1)
						, 0
						, 0
						, 0
						, 0
						, foldParams)
					// correct from dcal/mol to kcal/mol
					/ (E_type)100.0
				);
	}
	}

	return minStackingE;
}


////////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergyVrna::
getBestE_dangling() const
{
	// TODO maybe setup member variable (init=E_INF) with lazy initialization

	// get all possible base pair codes handled
	std::set<int> basePairCodes;
	for (int i=0; i<NBASES; i++) {
	for (int j=0; j<NBASES; j++) {
		basePairCodes.insert( BP_pair[i][j] );
	}
	}
	// get minimal energy for any base pair code combination
	E_type minDangleE = std::numeric_limits<E_type>::infinity();

	// get codes for the sequence alphabet
	const RnaSequence::CodeSeq_type alphabet = RnaSequence::getCodeForString(RnaSequence::SequenceAlphabet);

	// get minimal
	for (std::set<int>::const_iterator p1=basePairCodes.begin(); p1!=basePairCodes.end(); p1++) {
	for (size_t i=0; i<alphabet.size(); i++) {
	for (size_t j=0; j<alphabet.size(); j++) {
		minDangleE = std::min( minDangleE
				, (E_type)E_Stem( *p1
					  , alphabet.at(i)
					  , alphabet.at(j)
					  , 1 // is an external loop
					  , foldParams
					  )
					// correct from dcal/mol to kcal/mol
					/ (E_type)100.0
				);
	}
	}
	}

	return minDangleE;
}

////////////////////////////////////////////////////////////////////////////

