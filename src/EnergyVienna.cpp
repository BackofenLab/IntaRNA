/*
 * EnergyVienna.cpp
 *
 *  Created on: 27.06.2014
 *      Author: Mmann
 */

#include "EnergyVienna.h"

#include <cassert>

extern "C" {
	#include "ViennaRNA/utils.h"
	#include "ViennaRNA/fold_vars.h"
	#include "ViennaRNA/params.h"
	#include "ViennaRNA/pair_mat.h"
	#include <ViennaRNA/loop_energies.h>
}


////////////////////////////////////////////////////////////////////////////

EnergyVienna::EnergyVienna(
		const Accessibility & accS1
		, const ReverseAccessibility & accS2
		, const size_t maxInternalLoopSize1
		, const size_t maxInternalLoopSize2
		, const T_type T
	)
 :
	Energy(accS1, accS2, maxInternalLoopSize1, maxInternalLoopSize2)
	, modelDetails()
	, foldParams(NULL)
{
	// Vienna RNA : get default model details
	set_model_details(&modelDetails);

	// set OUR model details
	modelDetails.dangles = 2;
	modelDetails.noLP = 0;
	modelDetails.noGU = 0;
	modelDetails.noGUclosure = 0;
	modelDetails.circ = 0;
	modelDetails.gquad = 0;

	// Vienna RNA : get final folding parameters
	foldParams = get_scaled_parameters((double)T, modelDetails);

}

////////////////////////////////////////////////////////////////////////////

EnergyVienna::~EnergyVienna()
{
	// garbage collection
	if (foldParams != NULL) {
		delete foldParams;
		foldParams = NULL;
	}
}


////////////////////////////////////////////////////////////////////////////

E_type
EnergyVienna::
getInterLoopE( const size_t i1, const size_t j1, const size_t i2, const size_t j2 )
{
	// if region indices and loop sizes are ok
	if (Energy::isAllowedLoopRegion(accS1.getSequence(), i1, j1, maxInternalLoopSize1)
		&& Energy::isAllowedLoopRegion(accS2.getSequence(), i2, j2, maxInternalLoopSize2)
	) {
		if (i1 == j1 ) {
			assert( i2==j2 );
			// no internal loop --> no energy
			return (E_type)0.0;
		} else {
			assert( i2!=j2 );
			// Vienna RNA : compute internal loop / stacking energy for base pair [i1,i2]
			return (E_type)E_IntLoop(	(int)j1-i1-1	// unpaired region 1
								, (int)j2-i2-1	// unpaired region 2
								, BP_pair[accS1.getSequence().asCodes().at(i1)][accS1.getSequence().asCodes().at(i2)]	// type BP (i1,i2)
								, BP_pair[accS1.getSequence().asCodes().at(j2)][accS1.getSequence().asCodes().at(j1)]	// type BP (j2,j1)
								, accS1.getSequence().asCodes().at(i1+1)
								, accS2.getSequence().asCodes().at(i2+1)
								, accS1.getSequence().asCodes().at(j1-1)
								, accS2.getSequence().asCodes().at(j2-1)
								, foldParams)
					// correct from dcal/mol to kcal/mol
					/ (E_type)10.0
					;
		}
	} else {
		return E_MAX;
	}
}

////////////////////////////////////////////////////////////////////////////

E_type
EnergyVienna::
getDanglingLeft( const size_t i1, const size_t i2 )
{
	// Vienna RNA : dangling end contribution
	return (E_type) E_Stem( BP_pair[accS1.getSequence().asCodes().at(i1)][accS1.getSequence().asCodes().at(i2)]
							  , -1+(int)i1
							  , -1+(int)i2
							  , 1 // is an external loop  : TODO check if return value higher than ML-dangles (should be)
							  , foldParams
							  );
}

////////////////////////////////////////////////////////////////////////////

E_type
EnergyVienna::
getDanglingRight( const size_t j1, const size_t j2 )
{
	// Vienna RNA : dangling end contribution (reverse base pair to be sequence end conform)
	return (E_type) E_Stem( BP_pair[accS1.getSequence().asCodes().at(j2)][accS1.getSequence().asCodes().at(j1)]
							  , ( accS2.getSequence().size()-j2-1 == 0 ? -1 : (int)j2+1)  // check if last position in S2
							  , ( accS1.getSequence().size()-j1-1 == 0 ? -1 : (int)j1+1)  // check if last position in S1
							  , 1 // is an external loop  : TODO check if return value higher than ML-dangles (should be)
							  , foldParams
							  );
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

