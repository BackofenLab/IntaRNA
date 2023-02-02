#ifndef VIENNAHANDLER_H_
#define VIENNAHANDLER_H_

#include "IntaRNA/general.h"
#include <string>

extern "C" {
	#include <ViennaRNA/model.h>
	#include <ViennaRNA/params.h>
}

namespace IntaRNA {

/**
 * Central hub for Vienna RNA package related stuff
 *
 * @author Martin Mann
 *
 */
class VrnaHandler {

public:

	static constexpr const char* Turner99 = "Turner99";
	static constexpr const char* Turner04 = "Turner04";
	static constexpr const char* Andronescu07 = "Andronescu07";

protected:

	//! VRNA parameter model
	vrna_md_t model;

	//! the RT constant used for the current setup
	Z_type RT;

public:

	/**
	 * Construction with global VRNA folding parameter setup
	 *
	 * @param vrnaParamFile name of a VRNA parameter file to be used for
	 *        parameter setup
	 * @param temperature Folding temperature in Celsius
	 * @param noGUclosue whether or not GU base pairs are allowed at helix ends
	 * @param noLP whether or not lonely base pairs are considered
	 */
	VrnaHandler(
			Z_type temperature = 37.0
			, const std::string & vrnaParamFile = std::string(Turner04)
			, const bool noGUclosure = false
			, const bool noLP = false );

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
	 * @oaram pfScale can be used to scale pf_scale of VRNA partition function computation.
	 *                Note: only used for values >= 1.0, all other values are ignored.
	 * @return the model to be used for VRNA computations
	 */
	vrna_md_t
	getModel(
			const int max_bp_span = -1,
			const int window_size = -1,
			const double pfScale = getPfScaleDefault() ) const;

	/**
	 * Provides RT for the current setup
	 * @return R*temperature
	 */
	Z_type
	getRT() const;

	/**
	 * Provides RT for the given temperature
	 * @return R*temperature
	 */
	static
	Z_type
	getRT( const Z_type temperature );


	/**
	 * Provides the default value of the ViennaRNA package for the pf_scale factor
	 * model.sfact to be used in argument parsing etc.
	 * @return the default model.sfact value from the ViennaRNA package
	 */
	static
	double
	getPfScaleDefault();



};

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

inline
Z_type
VrnaHandler::
getRT() const
{
	return RT;
}

////////////////////////////////////////////////////////////////////////////

} // namespace

#endif /* VIENNAHANDLER_H_ */
