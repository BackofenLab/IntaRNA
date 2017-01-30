
#ifndef ACCESSIBILITYVIENNA_H_
#define ACCESSIBILITYVIENNA_H_

#include "Accessibility.h"
#include "VrnaHandler.h"

#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/triangular.hpp>

#include <iostream>


extern "C" {
	#include "ViennaRNA/fold_vars.h"
	#include "ViennaRNA/params.h"
}



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
	 * @param plFoldL the maximal base pair span to be used for plFold computations
	 * @param computeES whether or not ES values are to be computed
	 */
	AccessibilityVrna( const RnaSequence& sequence
			, const size_t maxLength
			, const AccessibilityConstraint * const accConstraint
			, const VrnaHandler & vrnaHandler
			, const size_t plFoldW = 0
			, const size_t plFoldL = 0
			, const bool computeES = false
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

	/**
	 * Provides the ensemble energy (ES) of all intramolecular substructures
	 * that can be formed within a given region of the sequence under the
	 * assumption that the region is part of an (intermolecular) multiloop,
	 * i.e. at least one base pair is formed by each substructure, given by
	 *
	 *   ES(i,j) = -RT * log( Qm[i,j] )
	 *
	 * where Qm denotes the according partition function computed by the
	 * McCaskill algorithm.
	 *
	 * If no structure can be formed within the region (Qm==0), E_INF is returned.
	 *
	 * @param i the start of the structured region
	 * @param j the end of the structured region
	 * @return the ES value for [i,j] or E_INF if no intramolecular
	 *         structure can be formed
	 */
	virtual
	E_type
	getES( const size_t i, const size_t j ) const;

protected:

	//! type for the ED value matrix (upper triangular matrix banded by maxLength)
	typedef boost::numeric::ublas::banded_matrix<E_type> EdMatrix;

	//! the ED values for the given sequence
	EdMatrix edValues;

	//! matrix to store ES values (upper triangular matrix)
	typedef boost::numeric::ublas::triangular_matrix<E_type, boost::numeric::ublas::upper> EsMatrix;

	//! the ES values for the given sequence if computed (otherwise NULL)
	EsMatrix * esValues;

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
	 * Computes the ES values and fills esValues container
	 *
	 * @param vrnaHandler the VrnaHandler to be used
	 * @param maxBpSpan the maximal distance of pairing positions
	 */
	void
	computeES( const VrnaHandler & vrnaHandler, const size_t maxBpSpan );

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
		// check for constrained positions within region
		if (!getAccConstraint().isAccessible(from,to)) {
			// position covers a blocked position --> omit accessibility
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

inline
E_type
AccessibilityVrna::
getES( const size_t from, const size_t to ) const
{
	if (esValues == NULL) {
		throw std::runtime_error("AccessibilityVrna::getES() : no ES values available (constructed without ES computation)");
	}
	// input range check
	checkIndices(from,to);
	return (*esValues)(from,to);
}

/////////////////////////////////////////////////////////////////////////////


#endif /* ACCESSIBILITYVIENNA_H_ */
