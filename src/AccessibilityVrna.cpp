
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
			, const VrnaHandler & vrnaHandler
			, const size_t maxLength
			, const size_t plFoldW
			, const size_t plFoldL
			, const AccessibilityConstraint * const accConstraint
			, const bool computeES_
		)
 :
	Accessibility( seq, maxLength, accConstraint ),
	edValues( getSequence().size(), getSequence().size(), 0, getMaxLength() ),
	esValues( NULL )
{

	// check if constraint given
	// or sliding window empty
	// or larger than sequence length
	if ( (! getAccConstraint().isEmpty()) || (plFoldW==0) || (plFoldW >= getSequence().size()) ) {
		if (plFoldW > 0 && plFoldW < getSequence().size() ) {
			throw std::runtime_error("sequence '"+seq.getId()+"': accuracy constraints provided but sliding window enabled (>0), which is currently not supported");
		}
		fillByRNAup(vrnaHandler
				, (plFoldL==0? getSequence().size() : std::min(plFoldL,getSequence().size()))
				);
		// compute ES values
		if (computeES_) {
			computeES( vrnaHandler, plFoldL );
		}
		// inefficient ED value computation O(n^2)*O(n^5) for debugging
//		fillByConstraints(vrnaHandler, (plFoldW==0? getSequence().size() : std::min(plFoldW,getSequence().size())), plFoldL);
	} else {
		if (computeES_) {
			throw std::runtime_error("can not compute local structure energies (ES) for sliding window accessibility computation");
		}
		fillByRNAplfold(vrnaHandler
				, (plFoldW==0? getSequence().size() : std::min(plFoldW,getSequence().size()))
				, (plFoldL==0? getSequence().size() : std::min(plFoldL,getSequence().size()))
				);
	}

}

/////////////////////////////////////////////////////////////////////////////

AccessibilityVrna::~AccessibilityVrna()
{
	if (esValues != NULL) {
		delete esValues;
	}
}

///////////////////////////////////////////////////////////////////////////////


E_type
AccessibilityVrna::
calc_ensemble_free_energy( const int start_unfold, const int end_unfold, vrna_exp_param_s * partFoldParams )
{
#if IN_DEBUG_MODE
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
getPfScale( const RnaSequence & seq, const VrnaHandler & vrnaHandler
		, const size_t plFoldL )
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

	// add maximal BP span
	vrna_md_t curModel = vrnaHandler.getModel( plFoldL, seq.size() );

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
fillByConstraints( const VrnaHandler &vrnaHandler
		, const size_t plFoldL )
{
	VLOG(2) <<"computing accessibility via n^2 fold calls...";
	// time logging
	TIMED_FUNC_IF(timerObj, VLOG_IS_ON(9));

	if (esValues != NULL) {
		NOTIMPLEMENTED("computation of local structure energies (ES) via AccessibilityVrna::fillByConstraints() to done");
	}

	// get scaling factor to avoid math problems in partition function computation
	const double pfScale = getPfScale( seq, vrnaHandler, plFoldL );

	// add maximal BP span
	vrna_md_t curModel = vrnaHandler.getModel( plFoldL, seq.size() );

	// Vienna RNA : get final partition function folding parameters
	vrna_exp_param_s* partFoldParams
		= get_boltzmann_factors( curModel.temperature
								, curModel.betaScale
								, curModel
								, pfScale);

	// compute free energy of whole structure ensemble
	E_type E_all = calc_ensemble_free_energy(-1, -1, partFoldParams);

	int count=0;
	const int seq_len = (int)getSequence().size();
	// compute ED values for _all_ regions [i,j]
	for(int i=0; i<seq_len; i++)
	{
		bool regionUnconstrained = getAccConstraint().isUnconstrained(i);
		// compute only for region lengths (j-i+1) <= maxLength
		for(int j=i; j<std::min(seq_len,i+(int)getMaxLength()); j++)
		{
			// extend knowledge about "unconstrainedness" for the region
			regionUnconstrained = regionUnconstrained && getAccConstraint().isUnconstrained(j);
			// check if unconstrained within region (i,j)
			if (regionUnconstrained) {
				// compute ED value = E(unstructured in [i,j]) - E_all
			} else {
				edValues (i,j) = (calc_ensemble_free_energy(i,j, partFoldParams) - E_all);
				// region covers constrained elements --> set to upper bound
				edValues (i,j) = ED_UPPER_BOUND;
			}

		}
	}

}

///////////////////////////////////////////////////////////////////////////////


void
AccessibilityVrna::
fillByRNAplfold( const VrnaHandler &vrnaHandler
		, const size_t plFoldW
		, const size_t plFoldL )
{
	VLOG(2) <<"computing accessibility via plfold routines...";
	// time logging
	TIMED_FUNC_IF(timerObj, VLOG_IS_ON(9));

#if IN_DEBUG_MODE
	// check if structure constraint given
	if ( ! getAccConstraint().isEmpty() ) {
		throw std::runtime_error("AccessibilityVrna::fillByRNAplfold() called but structure constraint present for sequence "+getSequence().getId());
	}
	if (plFoldW < 3) {
		throw std::runtime_error("AccessibilityVrna::fillByRNAplfold() : plFoldW < 3");
	}
#endif

	if (esValues != NULL) {
		throw std::runtime_error("can not compute local structure energies (ES) for sliding window accessibility computation");
	}

	// add maximal BP span
	vrna_md_t curModel = vrnaHandler.getModel( plFoldL, plFoldW );

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


	// call folding and unpaired prob calculation
    vrna_plist_t * pl = pfl_fold_par(sequence
    		, plFoldW // winsize
    		, (plFoldL==0? plFoldW : std::min(plFoldW,plFoldL)) // base pair distance
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
    	for (int i=std::max(1,j-(int)pup[0][0]);i<=j; i++) {
			// compute overall unpaired probability
			double prob_unpaired = pup[j][j-i+1];
			// compute ED value = E(unstructured in [i,j]) - E_all
			edValues (i-1,j-1) = -RT*std::log(prob_unpaired);
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
fillByRNAup( const VrnaHandler &vrnaHandler
			, const size_t plFoldL )
{
	VLOG(2) <<"computing accessibility via RNAup routines...";
	// time logging
	TIMED_FUNC_IF(timerObj, VLOG_IS_ON(9));

	const int seqLength = (int)getSequence().size();

	// get model
	vrna_md_t curModel = vrnaHandler.getModel( plFoldL, seqLength );

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

	char * sequence = (char *) vrna_alloc(sizeof(char) * (seqLength + 1));
	char * structure = (char *) vrna_alloc(sizeof(char) * (seqLength + 1));
	for (int i=0; i<seqLength; i++) {
		// copy sequence
		sequence[i] = getSequence().asString().at(i);
		// copy accessibility constraint if present
		structure[i] = getAccConstraint().getVrnaDotBracket(i);
	}
	sequence[seqLength] = structure[seqLength] = '\0';

	// get used RT from vienna package
	const double RT = vrnaHandler.getRT();

	// copy folding constraint
	// get mfe for this sequence (for scaling reasons)
	const double min_energy = fold( sequence, structure );
	// overwrite structure again with constraint since it now holds the mfe structure
	for (int i=0; i<seqLength; i++) { structure[i] = getAccConstraint().getVrnaDotBracket(i); }
//	// compute the partition function scaling value
//	const double pf_scale = std::exp(-(curModel.sfact*min_energy)/RT/seqLength);
	// compute partition function matrices to prepare unpaired probability calculation
	const float ensemble_energy = pf_fold(sequence, structure);
	// compute unpaired probabilities (cast-hack due to non-conform VRNA interface)
	pu_contrib* unstr_out = pf_unstru( sequence, (int)getMaxLength());

	// free folding data and Q matrices
	free_pf_arrays();
	free(structure);
	free(sequence);

	// check if any constraint present
	// compute ED values for _all_ regions [i,j]
	for (int i=(int)seqLength; i>0; i--) {
		bool regionUnconstrained = getAccConstraint().isUnconstrained(i-1);
		// compute only for region lengths (j-i+1) <= maxLength
		for(int j=i; j<=std::min((int)seqLength,unstr_out->w);j++)
		{
			// extend knowledge about "unconstrainedness" for the region
			regionUnconstrained = regionUnconstrained && (getAccConstraint().isUnconstrained(j-1));
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

}

////////////////////////////////////////////////////////////////////////////

void
AccessibilityVrna::
computeES( const VrnaHandler & vrnaHandler, const size_t maxBpSpan )
{

	// prepare container
	if (esValues==NULL) {
		esValues = new EsMatrix(getSequence().size(), getSequence().size());
	} else {
		esValues->resize( getSequence().size(), getSequence().size() );
	}

	// sequence length
	const int seqLength = (int)getSequence().size();
	const E_type RT = vrnaHandler.getRT();
	// folding parameters
	vrna_md_t curModel = vrnaHandler.getModel( maxBpSpan, seqLength );
	// VRNA compatible data structures
	char * sequence = (char *) vrna_alloc(sizeof(char) * (seqLength + 1));
	char * structureConstraint = (char *) vrna_alloc(sizeof(char) * (seqLength + 1));
	for (int i=0; i<seqLength; i++) {
		// copy sequence
		sequence[i] = getSequence().asString().at(i);
		// copy accessibility constraint if present
		structureConstraint[i] = getAccConstraint().getVrnaDotBracket(i);
	}
	sequence[seqLength] = structureConstraint[seqLength] = '\0';
	// prepare folding data
	vrna_fold_compound_t * foldData = vrna_fold_compound( sequence, &curModel, VRNA_OPTION_PF);

	// set accessibility constraint
    unsigned int constraint_options = 0;
    constraint_options    |= VRNA_CONSTRAINT_DB
                          |  VRNA_CONSTRAINT_DB_PIPE
                          |  VRNA_CONSTRAINT_DB_DOT
                          |  VRNA_CONSTRAINT_DB_X
                          |  VRNA_CONSTRAINT_DB_ANG_BRACK
                          |  VRNA_CONSTRAINT_DB_RND_BRACK;
    vrna_constraints_add( foldData, (const char *)structureConstraint, constraint_options);

    // compute correct partition function scaling via mfe
    FLT_OR_DBL min_free_energy = vrna_mfe( foldData, NULL );
    vrna_exp_params_rescale( foldData, &min_free_energy);

	// compute partition functions
	const float ensembleE = vrna_pf( foldData, NULL );

	if (foldData->exp_matrices == NULL) {
		throw std::runtime_error("AccessibilityVrna::computeES() : partition functions after computation not available");
	}
	// copy ensemble energies of multi loop parts = ES values
	FLT_OR_DBL qm_val = 0.0;
	for (int i=0; i<seqLength; i++) {
		for (int j=i; j<seqLength; j++) {
			// get Qm value
			// indexing via iindx starts with 1 instead of 0
			qm_val = foldData->exp_matrices->qm[foldData->iindx[i+1]-j+1];
			// ES energy = -RT*log( Qm )
			// ensure Qm > 0 for computation; otherwise E_INF
			(*esValues)(i,j) =  E_equal( qm_val, 0.0)
								? E_INF
								: std::min<E_type>( E_INF,
										(E_type)(-log(qm_val)
												-(E_type)seqLength*std::log(foldData->exp_params->pf_scale)
													)*foldData->exp_params->kT/1000.0);
		}
	}
	// garbage collection
	vrna_fold_compound_free(foldData);
	free(structureConstraint);
	free(sequence);

}

/////////////////////////////////////////////////////////////////////////////

