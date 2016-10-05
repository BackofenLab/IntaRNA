
#include "InteractionEnergyVrna.h"

#include <cassert>

extern "C" {
	#include <ViennaRNA/utils.h>
	#include <ViennaRNA/fold_vars.h>
	#include <ViennaRNA/params.h>
	#include <ViennaRNA/pair_mat.h>
	#include <ViennaRNA/loop_energies.h>
}

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
	, foldParams( vrna_params( &vrnaHandler.getModel() ) )
	, RT(vrnaHandler.getRT())
{
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
getInterLoopE( const size_t i1, const size_t j1, const size_t i2, const size_t j2 ) const
{
	// if valid internal loop
	if ( isValidInternalLoop(i1,j1,i2,j2) ) {
		if (i1 == j1 ) {
			assert( i2==j2 );
			// no internal loop --> duplex initialization
			return getBestInitEnergy();
		} else {
			assert( i2!=j2 );
			// Vienna RNA : compute internal loop / stacking energy for base pair [i1,i2]
			return (E_type)E_IntLoop(	(int)j1-i1-1	// unpaired region 1
								, (int)j2-i2-1	// unpaired region 2
								, BP_pair[accS1.getSequence().asCodes().at(i1)][accS2.getSequence().asCodes().at(i2)]	// type BP (i1,i2)
								, BP_pair[accS2.getSequence().asCodes().at(j2)][accS1.getSequence().asCodes().at(j1)]	// type BP (j2,j1)
								, accS1.getSequence().asCodes().at(i1+1)
								, accS2.getSequence().asCodes().at(i2+1)
								, accS1.getSequence().asCodes().at(j1-1)
								, accS2.getSequence().asCodes().at(j2-1)
								, foldParams)
					// correct from dcal/mol to kcal/mol
					/ (E_type)100.0
					;
		}
	} else {
		return E_INF;
	}
}

////////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergyVrna::
getDanglingLeft( const size_t i1, const size_t i2 ) const
{
	// Vienna RNA : dangling end contribution
	return (E_type) E_Stem( BP_pair[accS1.getSequence().asCodes().at(i1)][accS2.getSequence().asCodes().at(i2)]
							  , ( i1==0 ? -1 : accS1.getSequence().asCodes().at(i1-1) )
							  , ( i2==0 ? -1 : accS2.getSequence().asCodes().at(i2-1) )
							  , 1 // is an external loop
							  , foldParams
							  )
					// correct from dcal/mol to kcal/mol
							  /(E_type)100.0;
}

////////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergyVrna::
getDanglingRight( const size_t j1, const size_t j2 ) const
{
	// Vienna RNA : dangling end contribution (reverse base pair to be sequence end conform)
	return (E_type) E_Stem( BP_pair[accS2.getSequence().asCodes().at(j2)][accS1.getSequence().asCodes().at(j1)]
							  , ( j2+1==accS2.getSequence().size() ? -1 : accS2.getSequence().asCodes().at(j2+1) )
							  , ( j1+1==accS1.getSequence().size() ? -1 : accS1.getSequence().asCodes().at(j1+1) )
							  , 1 // is an external loop
							  , foldParams
							  )
					// correct from dcal/mol to kcal/mol
							  /(E_type)100.0;
}

////////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergyVrna::
getRT() const
{
	return RT;
}

////////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergyVrna::
getBestStackingEnergy() const
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
getBestInitEnergy() const
{
	return (E_type)foldParams->DuplexInit/(E_type)100.0;
}

////////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergyVrna::
getBestDangleEnergy() const
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

