
#include "InteractionEnergyVrna.h"

#include <cassert>
#include <set>

// ES computation
extern "C" {
	#include <ViennaRNA/fold_vars.h>
	#include <ViennaRNA/fold.h>
	#include <ViennaRNA/part_func.h>
	#include <ViennaRNA/structure_utils.h>
	#include <ViennaRNA/utils.h>
}


////////////////////////////////////////////////////////////////////////////

InteractionEnergyVrna::InteractionEnergyVrna(
		const Accessibility & accS1
		, const ReverseAccessibility & accS2
		, VrnaHandler &vrnaHandler
		, const size_t maxInternalLoopSize1
		, const size_t maxInternalLoopSize2
		, const bool initES
	)
 :
	InteractionEnergy(accS1, accS2, maxInternalLoopSize1, maxInternalLoopSize2)
// get final VRNA folding parameters
	, foldModel( vrnaHandler.getModel() )
	, foldParams( vrna_params( &foldModel ) )
	, RT(vrnaHandler.getRT())
	, bpCG( BP_pair[RnaSequence::getCodeForChar('C')][RnaSequence::getCodeForChar('G')] )
	, bpGC( BP_pair[RnaSequence::getCodeForChar('G')][RnaSequence::getCodeForChar('C')] )
	, esValues1(NULL)
	, esValues2(NULL)
{
	vrna_md_defaults_reset( &foldModel );

	// init ES values if needed
	if (initES) {
		// VRNA computation not completely threadsafe
#if INTARNA_MULITHREADING
		#pragma omp critical(intarna_omp_callingVRNA)
#endif
		{
			// create ES container to be filled
			esValues1 = new EsMatrix();
			esValues2 = new EsMatrix();
			// fill ES container
			computeES( accS1, *esValues1 );
			computeES( accS2, *esValues2 );
		} // omp critical(intarna_omp_callingVRNA)
	}
}

////////////////////////////////////////////////////////////////////////////

InteractionEnergyVrna::~InteractionEnergyVrna()
{
	// garbage collection
	if (foldParams != NULL) {
		free(foldParams);
		foldParams = NULL;
	}
	CLEANUP(esValues1);
	CLEANUP(esValues2);

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
		if (BP_pair[i][j] != 0) {
			basePairCodes.insert( BP_pair[i][j] );
		}
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

void
InteractionEnergyVrna::
computeES( const Accessibility & acc, InteractionEnergyVrna::EsMatrix & esValues )
{

	// prepare container
	esValues.resize( acc.getSequence().size(), acc.getSequence().size() );

	// sequence length
	const int seqLength = (int)acc.getSequence().size();
	const E_type RT = getRT();

	// VRNA compatible data structures
	char * sequence = (char *) vrna_alloc(sizeof(char) * (seqLength + 1));
	char * structureConstraint = (char *) vrna_alloc(sizeof(char) * (seqLength + 1));
	for (int i=0; i<seqLength; i++) {
		// copy sequence
		sequence[i] = acc.getSequence().asString().at(i);
		// copy accessibility constraint if present
		structureConstraint[i] = acc.getAccConstraint().getVrnaDotBracket(i);
	}
	sequence[seqLength] = structureConstraint[seqLength] = '\0';
	// prepare folding data
	vrna_md_t curModel;
	vrna_md_copy( &curModel, &foldModel );
	// set maximal base pair span
	curModel.max_bp_span = acc.getAccConstraint().getMaxBpSpan();
	if (curModel.max_bp_span == (int)acc.getSequence().size()) {
		curModel.max_bp_span = -1;
	}
	vrna_fold_compound_t * foldData = vrna_fold_compound( sequence, &foldModel, VRNA_OPTION_PF);

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
	if (foldData->exp_matrices->qm == NULL) {
		throw std::runtime_error("AccessibilityVrna::computeES() : partition functions Qm after computation not available");
	}
	// copy ensemble energies of multi loop parts = ES values
	FLT_OR_DBL qm_val = 0.0;
	const int minLoopSubseqLength = foldModel.min_loop_size + 2;
	for (int i=0; i<seqLength; i++) {
		for (int j=i; j<seqLength; j++) {
			// check if too short to enable a base pair
			if (j-i+1 < minLoopSubseqLength) {
				// make unfavorable
				esValues(i,j) = E_INF;
			} else {
				// get Qm value
				// indexing via iindx starts with 1 instead of 0
				qm_val = foldData->exp_matrices->qm[foldData->iindx[i+1]-j+1];
				if ( E_equal(qm_val, 0.) ) {
					esValues(i,j) = E_INF;
				} else {
					// ES energy = -RT*log( Qm )
					esValues(i,j) =  (E_type)( - RT*( std::log(qm_val)
													+((FLT_OR_DBL)(j-i+1))*std::log(foldData->exp_params->pf_scale)));
				}
			}
		}
	}
	// garbage collection
	vrna_fold_compound_free(foldData);
	free(structureConstraint);
	free(sequence);

}

////////////////////////////////////////////////////////////////////////////

