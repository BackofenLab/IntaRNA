/*
 * AccessibilityVienna.h
 *
 *  Created on: 25.06.2014
 *      Author: Mmann
 */

#ifndef ACCESSIBILITYVIENNA_H_
#define ACCESSIBILITYVIENNA_H_

#include "Accessibility.h"

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
class AccessibilityVienna : public Accessibility {

public:


	/**
	 * Construction
	 * @param sequence the sequence the accessibility data belongs to
	 * @param maxLength the maximal length of accessible regions to be
	 *           considered. 0 defaults to the sequence's length.
	 * @param accConstraint accessibility constraint that enforces some regions
	 *        to be unstructured both in sequence and interaction
	 * @param T temperature used to compute folding energies in Celsius
	 * @param log if not NULL, the stream to write log messages to
	 */
	AccessibilityVienna( const RnaSequence& sequence
			, const size_t maxLength = 0
			, const std::string & accConstraint = ""
			, const T_type T = 37.0
			, std::ostream * log = NULL
			);

	/**
	 * destruction
	 */
	virtual ~AccessibilityVienna();

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

	//! type for the ED value matrix (upper triangular matrix sufficient)
	typedef boost::numeric::ublas::triangular_matrix<E_type
				, boost::numeric::ublas::upper> EdMatrix;

	//! the ED values for the given sequence
	EdMatrix edValues;

	//! Vienna RNA package : energy model parameters
	model_detailsT modelDetails;

	//! Vienna RNA package : partition function folding parameters to be used
	//! for the energy computation
	pf_paramT* partFoldParams;


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
	 * @param T temperature used to compute folding energies in Celsius
	 *
	 * @return the energy of the structure ensemble
	 */
	E_type
	calc_ensemble_free_energy(
			const int start_unfold
			, const int end_unfold
			, const T_type T
			);

	/**
	 * Computes a scaling factor to avoid overflow in partition function
	 * computation.
	 *
	 * @param seq the sequence the parameter is for
	 * @param modelDetails the Vienna RNA energy model details to be used
	 * @param T the temperatue in Celsius to be used for energy computations
	 */
	static
	double
	getPfScale( const RnaSequence & seq
				, const model_detailsT & modelDetails
				, const T_type T );

};

#endif /* ACCESSIBILITYVIENNA_H_ */