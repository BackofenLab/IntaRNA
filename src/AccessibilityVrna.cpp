/*
 * AccessibilityVienna.cpp
 *
 *  Created on: 25.06.2014
 *      Author: Mmann
 */

#include "AccessibilityVrna.h"

#include <cmath>
#include <cassert>
#include <algorithm>
#include <limits>

extern "C" {
	#include <ViennaRNA/part_func.h>
	#include <ViennaRNA/fold.h>
}


/////////////////////////////////////////////////////////////////////////////

AccessibilityVrna::AccessibilityVrna(
			const RnaSequence& seq
			, const VrnaHandler & vrnaHandler
			, const size_t maxLength
			, const std::string & accConstraint
			, std::ostream * log
		)
 :
	Accessibility( seq, maxLength, accConstraint ),
	edValues( getSequence().size(), getSequence().size() ),
	modelDetails( vrnaHandler.getModel() ),
	partFoldParams( NULL )
{


	// get scaling factor to avoid math problems in partition function computation
	const double pfScale = getPfScale( seq, modelDetails );

	// Vienna RNA : get final partition function folding parameters
	partFoldParams = get_boltzmann_factors( modelDetails.temperature, modelDetails.betaScale, modelDetails, pfScale);


	// compute free energy of whole structure ensemble
	E_type E_all = calc_ensemble_free_energy(-1, -1);

	const std::string prefix =" ED("+getSequence().getId()+") computation: " ;
	int count=0;
	const int seq_len = (int)getSequence().size();
	// compute ED values for _all_ regions [i,j]
	for(int i=0; i<seq_len; i++)
	{
		bool regionUnconstrained = getAccConstraint().at(i) == '.';
		// compute only for region lengths (j-i+1) <= maxLength
		for(int j=i; j<std::min(seq_len,i+(int)maxLength); j++)
		{
			// extend knowledge about "unconstrainedness" for the region
			regionUnconstrained = regionUnconstrained && getAccConstraint().at(j) == '.';
			// check if unconstrained within region (i,j)
			if (regionUnconstrained) {
				// compute ED value = E(unstructured in [i,j]) - E_all
				edValues (i,j) = (calc_ensemble_free_energy(i,j) - E_all);
			} else {
				// region covers constrained elements --> set to upper bound
				edValues (i,j) = ED_UPPER_BOUND;
			}

			// Output progress if needed
			if( log != NULL)
			{
				count++;
				// give already processed part in percent (overwriting last)
				*log <<"\r" << prefix << count*200/(seq_len*(seq_len+1)) << "%" ;
				log->flush();
			}
		}
	}

	// Output progress if needed
	if( log != NULL)
	{
		// give already processed part in percent (overwriting last)
		*log <<"\r"<< prefix << "completed\n";
		log->flush();
	}

}

/////////////////////////////////////////////////////////////////////////////

AccessibilityVrna::~AccessibilityVrna()
{
	// garbage collection
	if (partFoldParams != NULL) {
		delete partFoldParams;
		partFoldParams = NULL;
	}

}

/////////////////////////////////////////////////////////////////////////////

E_type
AccessibilityVrna::
getED( const size_t from, const size_t to ) const
{
	// input range check
	checkIndices(from,to);

	if ((to-from+1) <= maxLength) {
		// return according ED value from the precomputed matrix
		return edValues (from,to);
	} else {
		// region length exceeds maximally allowed length -> no value
		return ED_UPPER_BOUND;
	}
}

///////////////////////////////////////////////////////////////////////////////


E_type
AccessibilityVrna::
calc_ensemble_free_energy( const int start_unfold, const int end_unfold )
{
#ifdef DEBUG
	if (start_unfold >= 0 && end_unfold >= 0) {
		checkIndices((size_t)start_unfold, (size_t)end_unfold);
	} else {
		if (start_unfold != -1 || end_unfold != -1) {
			throw std::runtime_error("AccessibilityVienna::calc_ensemble_free_energy : range ["+toString(start_unfold)+","+toString(end_unfold)+"] not allowed");
		}
	}
#endif


	// get sequence length
	int len = (int)getSequence().size();

	// generate structure constraint
	// ('.' = 'unconstrained' and 'x' = 'unstructured/unpaired')
	char c_structure[len+1];
	c_structure[len] = '\0';
	if (start_unfold < 0) {
		for (int i=0; i<len; i++) {
			// copy constraint from global accessibility constraint
			c_structure[i] = getAccConstraint().at(i);
		}
	} else {
		int i=0;
		// first : constrained by global accessibility constraint
		for (i=0; i<start_unfold; i++) {
			c_structure[i] = getAccConstraint().at(i);
		}
		// followed by unpaired region for interaction
		for (; i<=end_unfold; i++) {
			c_structure[i] = 'x';
		}
		// at the end : constrained by global accessibility constraint
		for (; i<len; i++) {
			c_structure[i] = getAccConstraint().at(i);
		}

	}

	// Vienna RNA : get free energy of structure ensemble via partition function
	// while applying the structure constraint for the unstructured region
	const double energy = pf_fold_par(	getSequence().asString().c_str()
										, c_structure
										, partFoldParams
										, 0
										, (start_unfold==-1?0:1)
										, 0
									);

	// memory cleanup
	free_pf_arrays();

	return (E_type)energy;
}

///////////////////////////////////////////////////////////////////////////////

double
AccessibilityVrna::
getPfScale( const RnaSequence & seq, const model_detailsT & modelDetails )
{

	// get sequence length
	int len = (int)seq.size();

	// generate structure constraint
	// ('.' = 'unconstrained' and 'x' = 'unstructured/unpaired')
	char c_structure[len+1];
	c_structure[len] = '\0';
	for (int i=0; i<len; i++) {
		// unconstrained
		c_structure[i] = '.';
	}

	// Vienna RNA : get final folding parameters
	paramT * foldParams = get_scaled_parameters( modelDetails.temperature, modelDetails);

	// Vienna RNA : get mfe value
	char structure[len+1];
	strncpy(structure, c_structure, len+1);
	const double min_en = fold_par(	seq.asString().c_str()
									, structure
									, foldParams
									, 0
									, 0
								);
	// memory cleanup
	free_arrays();

	// garbage collection
	delete foldParams;


	// compute a scaling factor to avoid overflow in partition function
	// computation using the constants defined in the Vienna RNA package
	// get RT = GASCONST(R in cal/mol) * (T(in Celsius) + K0(Celsius2KelvinShift)) / 1000 (cal->kcal)
	const double RT = GASCONST*(modelDetails.temperature+K0)/1000.0;
	// arbitrary factor to get an energy slightly below the mfe
	double sfact = modelDetails.sfact;

	// get scaling factor
	return std::exp(-(sfact*min_en)/RT/(double)len);

}

///////////////////////////////////////////////////////////////////////////////

