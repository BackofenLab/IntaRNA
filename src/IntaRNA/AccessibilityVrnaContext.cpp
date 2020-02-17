
#include "IntaRNA/AccessibilityVrnaContext.h"

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

AccessibilityVrnaContext::AccessibilityVrnaContext(
			const RnaSequence& seq
			, const size_t maxLength
			, const AccessibilityConstraint * const accConstraint
			, const VrnaHandler & vrnaHandler
			, const size_t plFoldW
			, const size_t maxInteriorSpan
		)
 :
	AccessibilityVrna( seq, maxLength, accConstraint, vrnaHandler, plFoldW ),
	edExteriorValues( getSequence().size(), getSequence().size(), 0, getMaxLength() ),
	maxInteriorSpan(maxInteriorSpan)
{
	// super constructor calls init()
}

/////////////////////////////////////////////////////////////////////////////

void
AccessibilityVrnaContext::
init(const VrnaHandler & vrnaHandler
		, const size_t plFoldW)
{
	// if sequence shows minimal length
	if (seq.size() > 4) {
		// window-based accessibility computation
		fillByRNAplfold(vrnaHandler
				, (plFoldW==0? getSequence().size() : std::min(plFoldW,getSequence().size()))
				, getAccConstraint().getMaxBpSpan()
				, &callbackForStorageExterior
				);
	} else {
		// init ED values for short sequences
		for (auto row = edValues.begin1(); row != edValues.end1(); row++) {
			for (auto ed = row.begin(); ed != row.end(); ed++) {
				*ed = 0;
			}
		}
		// init exterior-context ED values for short sequences
		for (auto row = edExteriorValues.begin1(); row != edExteriorValues.end1(); row++) {
			for (auto ed = row.begin(); ed != row.end(); ed++) {
				*ed = 0;
			}
		}
	}
}



/////////////////////////////////////////////////////////////////////////////

AccessibilityVrnaContext::~AccessibilityVrnaContext()
{
}

///////////////////////////////////////////////////////////////////////////////

void
AccessibilityVrnaContext::
callbackForStorageExterior(FLT_OR_DBL   *pr,
					 int          pr_size,
					 int          j, // right 3'-end of interval (indexing starting with 1)
					 int          max,
					 unsigned int type,
					 void         *data)
{
	// forward call to store full ED values
	AccessibilityVrna::callbackForStorage(pr,pr_size,j,max,type,data);

	// check if exterior data
	if (type & (VRNA_PROBS_WINDOW_UP | VRNA_EXT_LOOP)) {

		// access the storage data
		std::pair< AccessibilityVrnaContext*, FLT_OR_DBL > storageRT = *((std::pair< AccessibilityVrnaContext*, FLT_OR_DBL >*)data);
		// direct data access for computation
	    const FLT_OR_DBL RT = storageRT.second;
	    EdMatrix & edValues = storageRT.first->edExteriorValues;
	    const AccessibilityConstraint & accConstr = storageRT.first->getAccConstraint();

	    // copy unpaired data for all available interval lengths
	    // but ensure interval does not contain blocked positions
	    const bool rightEndBlocked = accConstr.isMarkedBlocked(j-1);
	    for (int l = std::min(j,std::min(pr_size,std::min(max,(int)storageRT.first->getMaxLength()))); l>=1; l--) {
			// get unpaired probability
			FLT_OR_DBL prob_unpaired = pr[l];
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
				edValues(i-1,j-1) = std::max<E_type>( 0., Z_2_E( -RT*Z_type(std::log(prob_unpaired) )));
			}
	    }

	} else {
#if INTARNA_IN_DEBUG_MODE
		LOG( DEBUG ) <<"AccessibilityVrnaContext::callbackForStorageExterior() : getting unexpected data for type " << type;
#endif
	}
}

/////////////////////////////////////////////////////////////////////////////

} // namespace

