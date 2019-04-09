
#include "IntaRNA/InteractionEnergyVrna.h"

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

namespace IntaRNA {

////////////////////////////////////////////////////////////////////////////

InteractionEnergyVrna::InteractionEnergyVrna(
		const Accessibility & accS1
		, const ReverseAccessibility & accS2
		, VrnaHandler &vrnaHandler
		, const size_t maxInternalLoopSize1
		, const size_t maxInternalLoopSize2
		, const bool initES
		, const E_type energyAdd
		, const bool energyWithDangles
	)
 :
	InteractionEnergy(accS1, accS2, maxInternalLoopSize1, maxInternalLoopSize2, energyAdd, energyWithDangles)
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
//	23.11.2017 : should not be relevant anymore
//#if INTARNA_MULITHREADING
//		#pragma omp critical(intarna_omp_callingVRNA)
//#endif
//		{
		// create ES container to be filled
		esValues1 = new EsMatrix();
		esValues2 = new EsMatrix();
		// fill ES container
		computeES( accS1, *esValues1 );
		computeES( accS2, *esValues2 );
//		} // omp critical(intarna_omp_callingVRNA)
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
	 INTARNA_CLEANUP(esValues1);
	 INTARNA_CLEANUP(esValues2);

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
				, Evrna_2_E(E_IntLoop(	0	// unpaired region 1
						, 0	// unpaired region 2
						, *p1	// type BP (i1,i2)
						, *p2	// type BP (j2,j1)
						, 0
						, 0
						, 0
						, 0
						, foldParams))
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
				, Evrna_2_E(vrna_E_ext_stem( *p1
					  , alphabet.at(i)
					  , alphabet.at(j)
					  , foldParams
					))
				);
	}
	}
	}

	return minDangleE;
}

////////////////////////////////////////////////////////////////////////////

void
InteractionEnergyVrna::
computeES( const Accessibility & acc, InteractionEnergyVrna::EsMatrix & esToFill )
{

	// prepare container
	esToFill.resize( acc.getSequence().size(), acc.getSequence().size() );

	// sequence length
	const int seqLength = (int)acc.getSequence().size();
	const Z_type RT = getRT();

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
	if (curModel.max_bp_span >= (int)acc.getSequence().size()) {
		curModel.max_bp_span = -1;
	}
	// TODO check if VRNA_OPTION_WINDOW reasonable to speedup
	vrna_fold_compound_t * foldData = vrna_fold_compound( sequence, &foldModel, VRNA_OPTION_PF);

	// Adding hard constraints from pseudo dot-bracket
	unsigned int constraint_options = VRNA_CONSTRAINT_DB_DEFAULT;
	// enforce constraints
	constraint_options |= VRNA_CONSTRAINT_DB_ENFORCE_BP;

	vrna_constraints_add( foldData, (const char *)structureConstraint, constraint_options);

    // compute correct partition function scaling via mfe
    FLT_OR_DBL min_free_energy = vrna_mfe( foldData, NULL );
    vrna_exp_params_rescale( foldData, &min_free_energy);

	// compute partition functions
	const FLT_OR_DBL ensembleE = vrna_pf( foldData, NULL );

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
				esToFill(i,j) = E_INF;
			} else {
				// get Qm value
				// indexing via iindx starts with 1 instead of 0
				qm_val = foldData->exp_matrices->qm[foldData->iindx[i+1]-j+1];
				if ( Z_equal(Z_type(qm_val), 0.) ) {
					esToFill(i,j) = E_INF;
				} else {
					// ES energy = -RT*log( Qm )
					esToFill(i,j) =  Z_2_E( - RT* Z_type( std::log(qm_val)
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


} // namespace
