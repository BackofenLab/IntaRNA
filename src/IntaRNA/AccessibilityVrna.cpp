
#include "IntaRNA/AccessibilityVrna.h"

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
	#include <ViennaRNA/constraints/SHAPE.h>
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

namespace IntaRNA {

/////////////////////////////////////////////////////////////////////////////

AccessibilityVrna::AccessibilityVrna(
			const RnaSequence& seq
			, const size_t maxLength
			, const AccessibilityConstraint * const accConstraint
			, const VrnaHandler & vrnaHandler
			, const size_t plFoldW
		)
 :
	Accessibility( seq, maxLength, accConstraint ),
	edValues( getSequence().size(), getSequence().size(), 0, getMaxLength() )
{

	// window-based accessibility computation
	fillByRNAplfold(vrnaHandler
			, (plFoldW==0? getSequence().size() : std::min(plFoldW,getSequence().size()))
			, getAccConstraint().getMaxBpSpan()
			);

// fillByRNAplfold computation not threadsafe
//#if INTARNA_MULITHREADING
//		#pragma omp critical(intarna_omp_callingVRNA)
//#endif
//{
//	fillByRNAplfold(..); // obsolete
//} // omp critical(intarna_omp_callingVRNA)

}

/////////////////////////////////////////////////////////////////////////////

AccessibilityVrna::~AccessibilityVrna()
{
}

///////////////////////////////////////////////////////////////////////////////


E_type
AccessibilityVrna::
calc_ensemble_free_energy( const int start_unfold, const int end_unfold, vrna_exp_param_s * partFoldParams )
{
#if INTARNA_IN_DEBUG_MODE
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
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) <<"computing accessibility via n^2 fold calls..."; }
	// time logging
	TIMED_FUNC_IF(timerObj, VLOG_IS_ON(9));

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
		const bool leftSideBlocked = getAccConstraint().isMarkedBlocked(i);
		// compute only for region lengths (j-i+1) <= maxLength
		for(int j=i; j<std::min(seq_len,i+(int)getMaxLength()); j++)
		{
			// check if ends are blocked
			if (leftSideBlocked || getAccConstraint().isMarkedBlocked(j)) {
				// region covers constrained elements --> set to upper bound
				edValues(i,j) = ED_UPPER_BOUND;
			} else {
				// compute ED value = E(unstructured in [i,j]) - E_all
				edValues(i,j) = std::max<E_type>(0.,(calc_ensemble_free_energy(i,j, partFoldParams) - E_all));
			}

		}
	}

}

///////////////////////////////////////////////////////////////////////////////

void
AccessibilityVrna::
callbackForStorage(FLT_OR_DBL   *pr,
					 int          pr_size,
					 int          j, // right 3'-end of interval (indexing starting with 1)
					 int          max,
					 unsigned int type,
					 void         *data)
{
	// check if expected data
	if (type & (VRNA_PROBS_WINDOW_UP | VRNA_ANY_LOOP)) {

		// access the storage data
		std::pair< AccessibilityVrna*, double > storageRT = *((std::pair< AccessibilityVrna*, double >*)data);
		// direct data access for computation
	    const double RT = storageRT.second;
	    EdMatrix & edValues = storageRT.first->edValues;
	    const AccessibilityConstraint & accConstr = storageRT.first->getAccConstraint();

	    // copy unpaired data for all available interval lengths
	    // but ensure interval does not contain blocked positions
	    const bool rightEndBlocked = accConstr.isMarkedBlocked(j-1);
	    for (int l = std::min(j,std::min(pr_size,std::min(max,(int)storageRT.first->getMaxLength()))); l>=1; l--) {
			// get unpaired probability
			double prob_unpaired = pr[l];
//			TODO: check for [0,1] range and correct if needed (print WARNING)
			// get left interval boundary index
			int i = j - l + 1;
			// check if interval ends are blocked positions
			// check if zero before computing its log-value
			if (rightEndBlocked || accConstr.isMarkedBlocked(i-1) || (prob_unpaired == 0.0) ) {
				// ED value = ED_UPPER_BOUND
				edValues(i-1,j-1) = ED_UPPER_BOUND;
			} else {
				// compute ED value = E(unstructured in [i,j]) - E_all
				edValues(i-1,j-1) = std::max<E_type>( 0., -RT*std::log(prob_unpaired));
			}
	    }

	} else {
#if INTARNA_IN_DEBUG_MODE
		LOG( DEBUG ) <<"AccessibilityVrna::callbackForStorage() : getting unexpected data for type " << type;
#endif
	}
}


///////////////////////////////////////////////////////////////////////////////


void
AccessibilityVrna::
fillByRNAplfold( const VrnaHandler &vrnaHandler
		, const size_t plFoldW
		, const size_t plFoldL )
{
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) <<"computing accessibility via plfold routines...";}
	// time logging
	TIMED_FUNC_IF(timerObj, VLOG_IS_ON(9));

#if INTARNA_IN_DEBUG_MODE
	if (plFoldW < 3) {
		throw std::runtime_error("AccessibilityVrna::fillByRNAplfold() : plFoldW < 3");
	}
#endif

	// add maximal BP span
	vrna_md_t curModel = vrnaHandler.getModel( plFoldL, plFoldW );

	const int length = getSequence().size();

	// copy sequence into C data structure
	char * sequence = (char *) vrna_alloc(sizeof(char) * (length + 1));
	for (int i=0; i<length; i++) {
		sequence[i] = getSequence().asString().at(i);
	}
	sequence[length] = '\0';

    // setup folding data
    vrna_fold_compound_t * fold_compound = vrna_fold_compound( sequence, &curModel, VRNA_OPTION_PF | VRNA_OPTION_WINDOW );

	// setup folding constraints
	if ( ! getAccConstraint().isEmpty() ) {
		// copy structure constraint
		char * structure = structure = (char *) vrna_alloc(sizeof(char) * (length + 1));
		for (int i=0; i<length; i++) {
		// copy accessibility constraint
		structure[i] = getAccConstraint().getVrnaDotBracket(i);
		}
		// set array end indicator
		structure[length] = '\0';

		// Adding hard constraints from pseudo dot-bracket
		unsigned int constraint_options = VRNA_CONSTRAINT_DB_DEFAULT;
		// enforce constraints
		constraint_options |= VRNA_CONSTRAINT_DB_ENFORCE_BP;

		// add constraint information to the fold compound object
		vrna_constraints_add(fold_compound, (const char *)structure, constraint_options);

		// cleanup
		free(structure);

		// check if SHAPE reactivity data available
		if (!getAccConstraint().getShapeFile().empty()) {

			// add SHAPE reactivity data
			// add SHAPE data as soft constraints
			vrna_constraints_add_SHAPE(fold_compound,
			                           getAccConstraint().getShapeFile().c_str(),
			                           getAccConstraint().getShapeMethod().c_str(), // shape_method,
			                           getAccConstraint().getShapeConversion().c_str(), // shape_conversion,
			                           0, // verbose,
			                           VRNA_OPTION_PF | VRNA_OPTION_WINDOW );


		}
	}


    // provide access to this object to be filled by the callback
    // and the normalized temperature for the Boltzmann weight computation
    std::pair< AccessibilityVrna*, double > storageRT(this, vrnaHandler.getRT());

	// call folding and unpaired prob calculation
    vrna_probs_window( fold_compound, plFoldW, VRNA_PROBS_WINDOW_UP, &callbackForStorage, (void*)(&storageRT));


    // garbage collection
    vrna_fold_compound_free(fold_compound);
    free(sequence);

}


///////////////////////////////////////////////////////////////////////////////


void
AccessibilityVrna::
fillByRNAup( const VrnaHandler &vrnaHandler
			, const size_t plFoldL )
{
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) <<"computing accessibility via RNAup routines..."; }
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
		const bool leftSideBlocked = getAccConstraint().isMarkedBlocked(i-1);
		// compute only for region lengths (j-i+1) <= maxLength
		for(int j=i; j<=std::min((int)seqLength,unstr_out->w);j++)
		{
			// check if region ends are blocked
			if (leftSideBlocked || getAccConstraint().isMarkedBlocked(j-1)) {
				// region ends blocked --> set to upper bound
				edValues(i-1,j-1) = ED_UPPER_BOUND;
			} else {
				// compute overall unpaired probability
				double prob_unpaired =
						unstr_out->H[i][j-i]+
						unstr_out->I[i][j-i]+
						unstr_out->M[i][j-i]+
						unstr_out->E[i][j-i];
				// check if zero before computing its log-value
				if ( prob_unpaired == 0.0 ) {
					// ED value = ED_UPPER_BOUND
					edValues(i-1,j-1) = ED_UPPER_BOUND;
				} else {
					// compute ED value = E(unstructured in [i,j]) - E_all
					edValues(i-1,j-1) = std::max<E_type>( 0., -RT*std::log(prob_unpaired));
				}
			}
		}
	}

	// free data
	free_pu_contrib( unstr_out );

}

/////////////////////////////////////////////////////////////////////////////

} // namespace

