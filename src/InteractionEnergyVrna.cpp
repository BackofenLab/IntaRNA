
#include "InteractionEnergyVrna.h"

#include <cassert>

extern "C" {
	#include <ViennaRNA/utils.h>
	#include <ViennaRNA/fold_vars.h>
	#include <ViennaRNA/params.h>
	#include <ViennaRNA/pair_mat.h>
	#include <ViennaRNA/loop_energies.h>
}


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
			// no internal loop --> no energy
			return (E_type)0.0;
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
		return E_MAX;
	}
}

////////////////////////////////////////////////////////////////////////////

E_type
InteractionEnergyVrna::
getDanglingLeft( const size_t i1, const size_t i2 ) const
{
	// Vienna RNA : dangling end contribution
	return (E_type) E_Stem( BP_pair[accS1.getSequence().asCodes().at(i1)][accS2.getSequence().asCodes().at(i2)]
							  , -1+(int)i1
							  , -1+(int)i2
							  , 1 // is an external loop  : TODO check if return value higher than ML-dangles (should be)
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
							  , ( accS2.getSequence().size()-j2-1 == 0 ? -1 : (int)j2+1)  // check if last position in S2
							  , ( accS1.getSequence().size()-j1-1 == 0 ? -1 : (int)j1+1)  // check if last position in S1
							  , 1 // is an external loop  : TODO check if return value higher than ML-dangles (should be)
							  , foldParams
							  )
					// correct from dcal/mol to kcal/mol
							  /(E_type)100.0;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

