/*
 * AccessibilityVienna.h
 *
 *  Created on: 25.06.2014
 *      Author: Mmann
 */

#ifndef ACCESSIBILITYVIENNA_H_
#define ACCESSIBILITYVIENNA_H_

#include "Accessibility.h"
#include "VrnaHandler.h"

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
	 * @param vrnaHandler the VRNA parameter handler to be used
	 * @param maxLength the maximal length of accessible regions to be
	 *           considered. 0 defaults to the sequence's length.
	 * @param accConstraint accessibility constraint that enforces some regions
	 *        to be unstructured both in sequence and interaction
	 * @param log if not NULL, the stream to write log messages to
	 */
	AccessibilityVrna( const RnaSequence& sequence
			, VrnaHandler & vrnaHandler
			, const size_t maxLength = 0
			, const std::string & accConstraint = ""
			, std::ostream * log = NULL
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

	//! type for the ED value matrix (upper triangular matrix sufficient)
	typedef boost::numeric::ublas::triangular_matrix<E_type
				, boost::numeric::ublas::upper> EdMatrix;

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
	 */
	static
	double
	getPfScale( const RnaSequence & seq
				, VrnaHandler & vrnaHandler );


	/**
	 * Use the old intaRNA way using n^2 constrained folding to fill ED-values
	 *
	 * @param vrnaHandler the VRNA handler to be used
	 * @param logStream stream to log debug output etc.
	 */
	void
	fillByConstraints( VrnaHandler &vrnaHandler, std::ostream * log );

	/**
	 * Use RNAup-like style to fill ED-values
	 *
	 * @param vrnaHandler the VRNA handler to be used
	 * @param logStream stream to log debug output etc.
	 */
	void
	fillByRNAup( VrnaHandler &vrnaHandler, std::ostream * log );

	/**
	 * Use RNAplfold-like style to fill ED-values
	 *
	 * @param vrnaHandler the VRNA handler to be used
	 * @param logStream stream to log debug output etc.
	 */
	void
	fillByRNAplfold( VrnaHandler &vrnaHandler, std::ostream * log );

};

#endif /* ACCESSIBILITYVIENNA_H_ */
