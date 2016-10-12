/*
 * ViennaSetup.h
 *
 *  Created on: 22.09.2016
 *      Author: Mmann
 */

#ifndef VIENNASETUP_H_
#define VIENNASETUP_H_

#include <string>

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

	//! VRNA parameter model
	vrna_md_t model;

public:

	/**
	 * Construction with global VRNA folding parameter setup
	 *
	 * @param temperature Folding temperature in Celsius
	 * @param vrnaParamFile name of a VRNA parameter file to be used for
	 *        parameter setup or NULL if defaults are to be used
	 */
	VrnaHandler(
			double temperature = 37.0,
			const std::string * const vrnaParamFile = NULL );

	/**
	 * destruction
	 */
	virtual ~VrnaHandler();

	/**
	 * Generates a new VRNA parameter model according to the global and local settings
	 *
	 * NOTE: you have to call vrna_md_defaults_reset() to broadcast
	 * the model details before using non-VRNA-API-3 functions.
	 *
	 * @param max_bp_span Maximal distance between base pair partners within one RNA
	 * 				or -1 if no constraint needed
	 * @param window_size Size of the sliding window for locally folding within one RNA
	 * 				or -1 if no constraint needed
	 * @return the model to be used for VRNA computations
	 */
	vrna_md_t
	getModel( int max_bp_span = -1, int window_size = -1 ) const;

	/**
	 * Provides RT for the current setup
	 * @return R*temperature
	 */
	double
	getRT() const;



};


#endif /* VIENNASETUP_H_ */
