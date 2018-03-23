
#ifndef INTARNA_ACCESSIBILITYVIENNA_H_
#define INTARNA_ACCESSIBILITYVIENNA_H_

#include "IntaRNA/Accessibility.h"
#include "IntaRNA/VrnaHandler.h"

#include <boost/numeric/ublas/banded.hpp>

#include <iostream>


extern "C" {
	#include <ViennaRNA/fold_vars.h>
	#include <ViennaRNA/params.h>
}


namespace IntaRNA {

/**
 * Computes accessibility energy values for _all_ _regions_ using the free
 * energies of structure ensembles based on partition function computations
 * via the Vienna RNA package.
 *
 * @author Martin Mann 2014
 */
class AccessibilityVrna : public Accessibility {

public:


	/**
	 * Construction
	 * @param sequence the sequence the accessibility data belongs to
	 * @param maxLength the maximal window size of accessible regions to be
	 *           considered. 0 defaults to the sequence's length.
	 * @param accConstraint if not NULL, accessibility constraint that enforces some regions
	 *        to be unstructured both in sequence and interaction
	 * @param vrnaHandler the VRNA parameter handler to be used
	 * @param plFoldW the sliding window size to be used for plFold computations
	 */
	AccessibilityVrna( const RnaSequence& sequence
			, const size_t maxLength
			, const AccessibilityConstraint * const accConstraint
			, const VrnaHandler & vrnaHandler
			, const size_t plFoldW = 0
			);

	/**
	 * destruction
	 */
	virtual ~AccessibilityVrna();

	/**
	 * Returns the accessibility energy value for the given range in the
	 * sequence, i.e. the energy difference (ED) to make the region accessible.
	 *
	 * @param from the start index of the regions (from <= to)
	 * @param to the end index of the regions (to < seq.length)
	 *
	 * @return the ED value if (j-1+1) <= maxLength or ED_UPPER_BOUND otherwise
	 *
	 * @throw std::runtime_error in case it does not hold 0 <= from <= to < seq.length
	 */
	virtual
	E_type
	getED( const size_t from, const size_t to ) const;

protected:

	//! type for the ED value matrix (upper triangular matrix banded by maxLength)
	typedef boost::numeric::ublas::banded_matrix<E_type> EdMatrix;

	//! the ED values for the given sequence
	EdMatrix edValues;

	/**
	 * Computes the free energy of the structure ensemble that is unstructured
	 * in the region [start_unfold,end_unfold] including the boundaries.
	 * If start and end are set to -1 the full structure ensemble without
	 * constraints is considered.
	 *
	 * @param start_unfold first position to be unstructured, or -1 if no
	 * structure constraint is to be set
	 * @param end_unfold last position to be unstructured, or -1 if no
	 * structure constraint is to be set
	 * @param partFoldParams the folding parameters to be used
	 *
	 * @return the energy of the structure ensemble
	 */
	E_type
	calc_ensemble_free_energy(
			const int start_unfold
			, const int end_unfold
			, vrna_exp_param_s * partFoldParams
			);

	/**
	 * Computes a scaling factor to avoid overflow in partition function
	 * computation.
	 *
	 * @param seq the sequence the parameter is for
	 * @param vrnaHandler the VRNA handler to be used
	 * @param plFoldL the maximal base pair span to be used or 0 for plFoldW
	 */
	double
	getPfScale( const RnaSequence & seq
				, const VrnaHandler & vrnaHandler
				, const size_t plFoldL );


	/**
	 * Use the old intaRNA way using n^2 constrained folding to fill ED-values
	 *
	 * @param vrnaHandler the VRNA handler to be used
	 * @param plFoldL the maximal base pair span to be used or 0 for plFoldW
	 */
	void
	fillByConstraints( const VrnaHandler &vrnaHandler
						, const size_t plFoldL );

	/**
	 * Use RNAup-like style to fill ED-values
	 *
	 * @param vrnaHandler the VRNA handler to be used
	 * @param plFoldL the maximal base pair span to be used or 0 for plFoldW
	 */
	void
	fillByRNAup( const VrnaHandler &vrnaHandler
					, const size_t plFoldL );

	/**
	 * Use RNAplfold-like style to fill ED-values
	 *
	 * @param vrnaHandler the VRNA handler to be used
	 * @param plFoldW the sliding window size to be used or 0 for full length
	 * @param plFoldL the maximal base pair span to be used or 0 for plFoldW
	 */
	void
	fillByRNAplfold( const VrnaHandler &vrnaHandler
						, const size_t plFoldW
						, const size_t plFoldL );

	/**
	 * callback function used when calling vrna_probs_window()
	 *
	 * This function will be called for each probability data set in the sliding
	 * window probability computation implementation of vrna_probs_window().
	 * The argument @a type specifies the type of probability that is passed to
	 * this function.
	 *
	 * @see vrna_probs_window()
	 *
	 * @param pr      An array of probabilities
	 * @param pr_size The length of the probability array
	 * @param j       The j-position (3'-end) of the probability intervals (indexing starting with 1)
	 * @param max     The (theoretical) maximum length of the probability array
	 * @param type    The type of probability that is passed to this function
	 * @param storageRT    Auxiliary data: should hold a pair(AccessibilityVrna*,double RT)
	 *
	 */
	static
	void
	callbackForStorage(	FLT_OR_DBL    *pr,
	                    int           pr_size,
	                    int           j,
	                    int           max,
	                    unsigned int  type,
	                    void          *storageRT);

};


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

inline
E_type
AccessibilityVrna::
getED( const size_t from, const size_t to ) const
{
	// input range check
	checkIndices(from,to);

	if ((to-from+1) <= getMaxLength()) {
		// check for constrained end positions
		if (!getAccConstraint().isAccessible(from) || !getAccConstraint().isAccessible(to)) {
			// end position blocked --> omit accessibility
			return ED_UPPER_BOUND;
		}
		// return according ED value from the precomputed matrix
		return edValues (from,to);
	} else {
		// region length exceeds maximally allowed length -> no value
		return ED_UPPER_BOUND;
	}
}

/////////////////////////////////////////////////////////////////////////////

} // namespace

#endif /* ACCESSIBILITYVIENNA_H_ */
