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

// constraint-based ED filling
extern "C" {
	#include <ViennaRNA/part_func.h>
	#include <ViennaRNA/fold.h>
}

// RNAup-like ED filling
extern "C" {
	#include <ViennaRNA/fold_vars.h>
	#include <ViennaRNA/fold.h>
	#include <ViennaRNA/part_func.h>
	#include <ViennaRNA/part_func_up.h>
	#include <ViennaRNA/utils.h>
}

// RNAplfold-like ED filling
extern "C" {
	#include <ViennaRNA/fold_vars.h>
	#include <ViennaRNA/fold.h>
	#include <ViennaRNA/part_func.h>
	#include <ViennaRNA/LPfold.h>
	#include <ViennaRNA/structure_utils.h>
	#include <ViennaRNA/utils.h>
}


/////////////////////////////////////////////////////////////////////////////

AccessibilityVrna::AccessibilityVrna(
			const RnaSequence& seq
			, VrnaHandler & vrnaHandler
			, const size_t maxLength
			, const std::string & accConstraint
			, std::ostream * log
		)
 :
	Accessibility( seq, maxLength, accConstraint ),
	edValues( getSequence().size(), getSequence().size(), 0, getMaxLength() )
{

	fillByConstraints(vrnaHandler,log);
	(*log) <<"\n################### constrained-based ##################\n" <<(*this) <<std::endl;
	fillByRNAup(vrnaHandler,log);
	(*log) <<"\n################### RNAup-based ##################\n" <<(*this) <<std::endl;
	fillByRNAplfold(vrnaHandler,log);
	(*log) <<"\n################### RNAplfold-based ##################\n" <<(*this) <<std::endl;

	(*log) <<"\n################### done ##################\n" <<std::endl;
}

/////////////////////////////////////////////////////////////////////////////

AccessibilityVrna::~AccessibilityVrna()
{

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
calc_ensemble_free_energy( const int start_unfold, const int end_unfold, vrna_exp_param_s * partFoldParams )
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
getPfScale( const RnaSequence & seq, VrnaHandler & vrnaHandler )
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
	paramT * foldParams = get_scaled_parameters( vrnaHandler.getModel().temperature, vrnaHandler.getModel());

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
	free(foldParams);

	// compute a scaling factor to avoid overflow in partition function
	return std::exp(-(vrnaHandler.getModel().sfact*min_en)/ vrnaHandler.getRT() /(double)len);

}

///////////////////////////////////////////////////////////////////////////////


void
AccessibilityVrna::
fillByConstraints( VrnaHandler &vrnaHandler, std::ostream * log )
{

	// get scaling factor to avoid math problems in partition function computation
	const double pfScale = getPfScale( seq, vrnaHandler );

	// Vienna RNA : get final partition function folding parameters
	vrna_exp_param_s* partFoldParams
		= get_boltzmann_factors( vrnaHandler.getModel().temperature
								, vrnaHandler.getModel().betaScale
								, vrnaHandler.getModel()
								, pfScale);


	// compute free energy of whole structure ensemble
	E_type E_all = calc_ensemble_free_energy(-1, -1, partFoldParams);

//	const std::string prefix =" ED("+getSequence().getId()+") computation: " ;
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
				edValues (i,j) = (calc_ensemble_free_energy(i,j, partFoldParams) - E_all);
			} else {
				// region covers constrained elements --> set to upper bound
				edValues (i,j) = ED_UPPER_BOUND;
			}

//			// Output progress if needed
//			if( log != NULL)
//			{
//				count++;
//				// give already processed part in percent (overwriting last)
//				*log <<"\r" << prefix << count*200/(seq_len*(seq_len+1)) << "%" ;
//				log->flush();
//			}
		}
	}

//	// Output progress if needed
//	if( log != NULL)
//	{
//		// give already processed part in percent (overwriting last)
//		*log <<"\r"<< prefix << "completed\n";
//		log->flush();
//	}

}

///////////////////////////////////////////////////////////////////////////////


void
AccessibilityVrna::
fillByRNAplfold( VrnaHandler &vrnaHandler, std::ostream * log )
{

	char * sequence = (char *) vrna_alloc(sizeof(char) * (seq.size() + 1));
	char * structure = (char *) vrna_alloc(sizeof(char) * (seq.size() + 1));
	for (int i=0; i<seq.size(); i++) {
		sequence[i] = seq.asString().at(i);
		structure[i] = this->getAccConstraint().at(i);
	}
	sequence[seq.size()] = structure[seq.size()] = '\0';

	int length = seq.size();

	double ** pup = NULL;
    pup       =(double **)  vrna_alloc((length+1)*sizeof(double *));
    pup[0]    =(double *)   vrna_alloc(sizeof(double)); /*I only need entry 0*/
    pup[0][0] = std::min((int)(maxLength>0?maxLength:length),length); // length of unpaired stretch in window

    vrna_plist_t * dpp = NULL; // ?? whatever..

    vrna_exp_param_t * pf_parameters = vrna_exp_params(& vrnaHandler.getModel());

	// call folding and unpaired prob calculation
    vrna_plist_t * pl = pfl_fold_par(sequence
    		, length // winsize
    		, length // base pair distance
    		, 0.0 // printing probability cut-off
    		, pup // the unpaired probabilities to be filled
    		, &dpp
    		, NULL // pUfp
    		, NULL // spup
    		, pf_parameters
    );

    const double RT = vrnaHandler.getRT();

    // copy data
    for (int j=1; j<=length; j++) {
		bool regionUnconstrained = getAccConstraint().at(j-1) == '.';
    	for (int i=std::max(1,j-(int)pup[0][0]);i<=j; i++) {
			// extend knowledge about "unconstrainedness" for the region
			regionUnconstrained = regionUnconstrained && getAccConstraint().at(i-1) == '.';
			// check if unconstrained within region (i,j)
			if (regionUnconstrained) {
				// compute overall unpaired probability
				double prob_unpaired = pup[j][j-i+1];
				// compute ED value = E(unstructured in [i,j]) - E_all
				edValues (i-1,j-1) = -RT*std::log(prob_unpaired);

			} else {
				// region covers constrained elements --> set to upper bound
				edValues (i-1,j-1) = ED_UPPER_BOUND;
			}

    	}
    }


    // garbage collection
    if (dpp!=NULL) free(dpp);
    free(pf_parameters);
    for (int i=1; i<=length; i++) {
    	free(pup[i]);
    }
    free(pup[0]);
    free(pup);
    free(structure);
    free(sequence);

}

///////////////////////////////////////////////////////////////////////////////


void
AccessibilityVrna::
fillByRNAup( VrnaHandler &vrnaHandler, std::ostream * log )
{

		// make model parameters accessible in VRNA2-API
		vrna_md_defaults_reset( &(vrnaHandler.getModel()) );
	    fold_constrained=1;
	    tetra_loop=1;
	    noLonelyPairs = 0;
	    noGU = 0;
	    no_closingGU = 0;
	    energy_set = 0;

		////////  RNAup-like (VRNA2-API) unpaired probability calculation  ///////

		char * sequence = (char *) vrna_alloc(sizeof(char) * (seq.size() + 1));
		char * structure = (char *) vrna_alloc(sizeof(char) * (seq.size() + 1));
		for (int i=0; i<seq.size(); i++) {
			sequence[i] = seq.asString().at(i);
			structure[i] = this->getAccConstraint().at(i);
		}
		sequence[seq.size()] = structure[seq.size()] = '\0';

		// get used RT from vienna package
		const double RT = vrnaHandler.getRT();

		// copy folding constraint
		// get mfe for this sequence (for scaling reasons)
	    const double min_energy = fold( sequence, structure );
	    // overwrite structure again with constraint since it now holds the mfe structure
	    for (int i=0; i<seq.size(); i++) { structure[i] = this->getAccConstraint().at(i); }
	    // compute the partition function scaling value
	    const double pf_scale = std::exp(-(vrnaHandler.getModel().sfact*min_energy)/RT/seq.size());
	    // compute partition function matrices to prepare unpaired probability calculation
	    const double ensemble_energy = pf_fold(sequence, structure);
	    // compute unpaired probabilities (cast-hack due to non-conform VRNA interface)
	    pu_contrib* unstr_out = pf_unstru( sequence, (int)maxLength);
	    // free folding data and Q matrices
	    free_pf_arrays();
	    free(structure);
	    free(sequence);

		// compute ED values for _all_ regions [i,j]
	    for (int i=(int)seq.size(); i>0; i--) {
			bool regionUnconstrained = getAccConstraint().at(i-1) == '.';
			// compute only for region lengths (j-i+1) <= maxLength
			for(int j=i; j<=std::min((int)seq.size(),i+(int)maxLength)-1;j++)
			{
				// extend knowledge about "unconstrainedness" for the region
				regionUnconstrained = regionUnconstrained && getAccConstraint().at(j-1) == '.';
				// check if unconstrained within region (i,j)
				if (regionUnconstrained) {
					// compute overall unpaired probability
					double prob_unpaired =
							unstr_out->H[i][j-i]+
							unstr_out->I[i][j-i]+
							unstr_out->M[i][j-i]+
							unstr_out->E[i][j-i];
					// compute ED value = E(unstructured in [i,j]) - E_all
					edValues (i-1,j-1) = -RT*std::log(prob_unpaired);

				} else {
					// region covers constrained elements --> set to upper bound
					edValues (i-1,j-1) = ED_UPPER_BOUND;
				}
			}
		}

	    // free data
	    free_pu_contrib( unstr_out );

//		// Output progress if needed
//		if( log != NULL)
//		{
//			const std::string prefix =" ED("+getSequence().getId()+") computation: " ;
//			// give already processed part in percent (overwriting last)
//			*log <<"\r"<< prefix << "completed\n";
//			log->flush();
//		}

}

/////////////////////////////////////////////////////////////////////////////

