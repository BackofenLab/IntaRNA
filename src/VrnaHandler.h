/*
 * ViennaSetup.h
 *
 *  Created on: 22.09.2016
 *      Author: Mmann
 */

#ifndef VIENNASETUP_H_
#define VIENNASETUP_H_


extern "C" {
	#include <ViennaRNA/model.h>
	#include <ViennaRNA/params.h>
}


/**
 * Central hub for Vienna RNA package related stuff
 *
 * @author Martin Mann
 *
 */
class VrnaHandler {

protected:

	//! Folding temperature in Celsius
	double temperature;

	//! Maximal distance between base pair partners within one RNA
	//! or -1 if no constraint needed
	int max_bp_span;

	//! Size of the sliding window for locally folding within one RNA
	//! or -1 if no constraint needed
	int window_size;


public:

	/**
	 * Construction with global VRNA folding parameter setup
	 *
	 * @param temperature Folding temperature in Celsius
	 * @param max_bp_span Maximal distance between base pair partners within one RNA
	 * 				or -1 if no constraint needed
	 * @param window_size Size of the sliding window for locally folding within one RNA
	 * 				or -1 if no constraint needed
	 */
	VrnaHandler(
			double temperature = 37.0,
			int max_bp_span = -1,
			int window_size = -1 );

	/**
	 * destruction
	 */
	virtual ~VrnaHandler();

	/**
	 * Generates a new VRNA parameter model according to the global settings
	 * @return the model to be used for VRNA computations
	 */
	vrna_md_t
	getModel() const;


};


#endif /* VIENNASETUP_H_ */
