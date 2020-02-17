
#ifndef INTARNA_ACCESSIBILITYVIENNACONTEXT_H_
#define INTARNA_ACCESSIBILITYVIENNACONTEXT_H_

#include "IntaRNA/AccessibilityVrna.h"




namespace IntaRNA {

/**
 * Computes accessibility energy values for _all_ _regions_ using the free
 * energies of structure ensembles based on partition function computations
 * via the Vienna RNA package.
 *
 * @author Martin Mann 2014
 */
class AccessibilityVrnaContext : public AccessibilityVrna {

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
	 * @param maxInteriorSpan for sequences of lengths below the given range
	 *           full ED values are provided; otherwise only exterior-context
	 *           EDs are given
	 */
	AccessibilityVrnaContext( const RnaSequence& sequence
			, const size_t maxLength
			, const AccessibilityConstraint * const accConstraint
			, const VrnaHandler & vrnaHandler
			, const size_t plFoldW = 0
			, const size_t maxInteriorSpan = 999999
			);

	/**
	 * destruction
	 */
	virtual ~AccessibilityVrnaContext();

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

	//! the exterior-context ED values for for the given sequence
	EdMatrix edExteriorValues;

	//! maximal sequence length for which full EDs are provided; above, only
	//! exterior-context EDs are given
	const size_t maxInteriorSpan;

	virtual
	void
	init(const VrnaHandler & vrnaHandler
			, const size_t plFoldW);

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
	callbackForStorageExterior(	FLT_OR_DBL    *pr,
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
AccessibilityVrnaContext::
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
		if ( (to-from+1) <= maxInteriorSpan) {
			// return full ED value from the precomputed matrix
			return edValues (from,to);
		} else {
			// return external-only ED value from the precomputed matrix
			return edExteriorValues (from,to);
		}
	} else {
		// region length exceeds maximally allowed length -> no value
		return ED_UPPER_BOUND;
	}
}

/////////////////////////////////////////////////////////////////////////////

} // namespace

#endif /* INTARNA_ACCESSIBILITYVIENNACONTEXT_H_ */
