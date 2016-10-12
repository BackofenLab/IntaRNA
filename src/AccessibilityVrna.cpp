
#include "AccessibilityVrna.h"

#include <cmath>
#include <cassert>
#include <algorithm>
#include <limits>
#include <stdexcept>

// constraint-based ED filling
extern "C" {
	#include <ViennaRNA/part_func.h>
	#include <ViennaRNA/fold.h>
	#include <ViennaRNA/model.h>
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
			, const size_t plFoldW
			, const size_t plFoldL
			, const AccessibilityConstraint * const accConstraint
			, std::ostream * log
		)
 :
	Accessibility( seq, maxLength, accConstraint ),
	edValues( getSequence().size(), getSequence().size(), 0, getMaxLength() )
{

	// check if constraint given
	if ( ! getAccConstraint().isEmpty() ) {
		if (plFoldW > 0 || plFoldL > 0) {
			throw std::runtime_error("AccessibilityVrna() : cannot compute accessibility with sliding window when constraints given");
		}
		fillByRNAup(vrnaHandler,log);
//		(*log) <<"\n################### RNAup-based ##################\n" <<(*this) <<std::endl;
		// old inefficient intarna-1 ED value computation
//		fillByConstraints(vrnaHandler,log);
//		(*log) <<"\n################### constrained-based ##################\n" <<(*this) <<std::endl;
	} else {
		fillByRNAplfold(vrnaHandler, plFoldW, plFoldL, log);
//		(*log) <<"\n################### RNAplfold-based ##################\n" <<(*this) <<std::endl;
	}
//	(*log) <<"\n################### done ##################\n" <<std::endl;

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

	if ((to-from+1) <= getMaxLength()) {
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
#ifdef NDEBUG           /* required by ANSI standard */
	// no check
#else
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
			c_structure[i] = getAccConstraint().getVrnaDotBracket(i);
		}
	} else {
		int i=0;
		// first : constrained by global accessibility constraint
		for (i=0; i<start_unfold; i++) {
			c_structure[i] = getAccConstraint().getVrnaDotBracket(i);
		}
		// followed by unpaired region for interaction
		for (; i<=end_unfold; i++) {
			c_structure[i] = 'x';
		}
		// at the end : constrained by global accessibility constraint
		for (; i<len; i++) {
			c_structure[i] = getAccConstraint().getVrnaDotBracket(i);
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
	int len = (int)getSequence().size();

	// generate structure constraint
	// ('.' = 'unconstrained' and 'x' = 'unstructured/unpaired')
	char c_structure[len+1];
	c_structure[len] = '\0';
	for (int i=0; i<len; i++) {
		// copy global accessibility constraint if present
		c_structure[i] = getAccConstraint().getVrnaDotBracket(i);
	}

	// TODO add maximal BP span
	vrna_md_t curModel = vrnaHandler.getModel( getMaxLength(), getMaxLength() );

	// Vienna RNA : get final folding parameters
	paramT * foldParams = get_scaled_parameters( curModel.temperature, curModel );

	// Vienna RNA : get mfe value
	char structure[len+1];
	strncpy(structure, c_structure, len+1);
	const double min_en = fold_par(	getSequence().asString().c_str()
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
	return std::exp(-(curModel.sfact*min_en)/ vrnaHandler.getRT() /(double)len);

}

///////////////////////////////////////////////////////////////////////////////


void
AccessibilityVrna::
fillByConstraints( VrnaHandler &vrnaHandler, std::ostream * log )
{

	// get scaling factor to avoid math problems in partition function computation
	const double pfScale = getPfScale( seq, vrnaHandler );

	// TODO add maximal BP span
	vrna_md_t curModel = vrnaHandler.getModel( getMaxLength(), getMaxLength() );

	// Vienna RNA : get final partition function folding parameters
	vrna_exp_param_s* partFoldParams
		= get_boltzmann_factors( curModel.temperature
								, curModel.betaScale
								, curModel
								, pfScale);


	// compute free energy of whole structure ensemble
	E_type E_all = calc_ensemble_free_energy(-1, -1, partFoldParams);

//	const std::string prefix =" ED("+getSequence().getId()+") computation: " ;
	int count=0;
	const int seq_len = (int)getSequence().size();
	// compute ED values for _all_ regions [i,j]
	for(int i=0; i<seq_len; i++)
	{
		bool regionUnconstrained = (getAccConstraint().getVrnaDotBracket(i) == '.');
		// compute only for region lengths (j-i+1) <= maxLength
		for(int j=i; j<std::min(seq_len,i+(int)getMaxLength()); j++)
		{
			// extend knowledge about "unconstrainedness" for the region
			regionUnconstrained = regionUnconstrained && (getAccConstraint().getVrnaDotBracket(j) == '.');
			// check if unconstrained within region (i,j)
			if (regionUnconstrained) {
				// compute ED value = E(unstructured in [i,j]) - E_all
			} else {
				edValues (i,j) = (calc_ensemble_free_energy(i,j, partFoldParams) - E_all);
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
fillByRNAplfold( VrnaHandler &vrnaHandler
		, const size_t plFoldW
		, const size_t plFoldL
		, std::ostream * log )
{
#ifdef NDEBUG           /* required by ANSI standard */
	// no check
#else
	// check if structure constraint given
	if ( ! getAccConstraint().isEmpty() ) {
		throw std::runtime_error("AccessibilityVrna::fillByRNAplfold() called but structure constraint present for sequence "+getSequence().getId());
	}
	if ( plFoldW > 0 && plFoldW < 3) {
		throw std::runtime_error("AccessibilityVrna::fillByRNAplfold() : plFoldW > 0 && plFoldW < 3");
	}
#endif

	// TODO add maximal BP span
	vrna_md_t curModel = vrnaHandler.getModel( getMaxLength(), getMaxLength() );

	// copy sequence into C data structure
	char * sequence = (char *) vrna_alloc(sizeof(char) * (getSequence().size() + 1));
	for (int i=0; i<getSequence().size(); i++) {
		sequence[i] = getSequence().asString().at(i);
	}
	sequence[getSequence().size()] = '\0';

	int length = getSequence().size();

	double ** pup = NULL;
    pup       =(double **)  vrna_alloc((length+1)*sizeof(double *));
    pup[0]    =(double *)   vrna_alloc(sizeof(double)); /*I only need entry 0*/
    pup[0][0] = (int)getMaxLength(); // length of unpaired stretch in window

    vrna_plist_t * dpp = NULL; // ?? whatever..

    vrna_exp_param_t * pf_parameters = vrna_exp_params(& curModel);

    int slidingWindow = (plFoldW==0? -1 : std::min(plFoldW,getSequence().size()));

	// call folding and unpaired prob calculation
    vrna_plist_t * pl = pfl_fold_par(sequence
    		, slidingWindow // winsize
    		, (plFoldL==0? slidingWindow : (int)plFoldL) // base pair distance
    		, 0.0 // printing probability cut-off
    		, pup // the unpaired probabilities to be filled
    		, &dpp
    		, NULL // pUfp
    		, NULL // spup
    		, pf_parameters
    );

    const double RT = vrnaHandler.getRT();

//    // check if any constraint present
//    const bool noConstraintPresent = getAccConstraint().size() == 0;
    // copy data
    for (int j=1; j<=length; j++) {
//		bool regionUnconstrained = noConstraintPresent || getAccConstraint().at(j-1) == '.';
    	for (int i=std::max(1,j-(int)pup[0][0]);i<=j; i++) {
//			// extend knowledge about "unconstrainedness" for the region
//			regionUnconstrained = regionUnconstrained && (noConstraintPresent || getAccConstraint().at(i-1) == '.');
//			// check if unconstrained within region (i,j)
//			if (regionUnconstrained) {
				// compute overall unpaired probability
				double prob_unpaired = pup[j][j-i+1];
				// compute ED value = E(unstructured in [i,j]) - E_all
				edValues (i-1,j-1) = -RT*std::log(prob_unpaired);

//			} else {
//				// region covers constrained elements --> set to upper bound
//				edValues (i-1,j-1) = ED_UPPER_BOUND;
//			}

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
    free(sequence);

}

///////////////////////////////////////////////////////////////////////////////


void
AccessibilityVrna::
fillByRNAup( VrnaHandler &vrnaHandler, std::ostream * log )
{
	// TODO add maximal BP span
	vrna_md_t curModel = vrnaHandler.getModel( getMaxLength(), getMaxLength() );

	// make model parameters accessible in VRNA2-API
	vrna_md_defaults_reset( &curModel );
	fold_constrained=1;
	tetra_loop=1;
	noLonelyPairs = curModel.noLP;
	noGU = curModel.noGU;
	no_closingGU = curModel.noGUclosure;
	energy_set = curModel.energy_set;

	// TODO CHECK IF TO BE CALLED OR NOT
//	update_fold_params();
//	    vrna_params_subst();

	////////  RNAup-like (VRNA2-API) unpaired probability calculation  ///////

	char * sequence = (char *) vrna_alloc(sizeof(char) * (getSequence().size() + 1));
	char * structure = (char *) vrna_alloc(sizeof(char) * (getSequence().size() + 1));
	for (int i=0; i<getSequence().size(); i++) {
		// copy sequence
		sequence[i] = getSequence().asString().at(i);
		// copy accessibility constraint if present
		structure[i] = getAccConstraint().getVrnaDotBracket(i);
	}
	sequence[getSequence().size()] = structure[getSequence().size()] = '\0';

	// get used RT from vienna package
	const double RT = vrnaHandler.getRT();

	// copy folding constraint
	// get mfe for this sequence (for scaling reasons)
	const double min_energy = fold( sequence, structure );
	// overwrite structure again with constraint since it now holds the mfe structure
	for (int i=0; i<getSequence().size(); i++) { structure[i] = getAccConstraint().getVrnaDotBracket(i); }
//	// compute the partition function scaling value
//	const double pf_scale = std::exp(-(curModel.sfact*min_energy)/RT/getSequence().size());
	// compute partition function matrices to prepare unpaired probability calculation
	const double ensemble_energy = pf_fold(sequence, structure);
	// compute unpaired probabilities (cast-hack due to non-conform VRNA interface)
	pu_contrib* unstr_out = pf_unstru( sequence, (int)getMaxLength());
	// free folding data and Q matrices
	free_pf_arrays();
	free(structure);
	free(sequence);

	// check if any constraint present
	// compute ED values for _all_ regions [i,j]
	for (int i=(int)getSequence().size(); i>0; i--) {
		bool regionUnconstrained = getAccConstraint().getVrnaDotBracket(i-1) == '.';
		// compute only for region lengths (j-i+1) <= maxLength
		for(int j=i; j<=std::min((int)getSequence().size(),i+(int)getMaxLength())-1;j++)
		{
			// extend knowledge about "unconstrainedness" for the region
			regionUnconstrained = regionUnconstrained && (getAccConstraint().getVrnaDotBracket(j-1) == '.');
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

