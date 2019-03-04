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

protected:

	//! VRNA parameter model
	vrna_md_t model;

	//! the RT constant used for the current setup
	Z_type RT;

public:

	/**
	 * Construction with global VRNA folding parameter setup
	 *
	 * @param temperature Folding temperature in Celsius
	 * @param vrnaParamFile name of a VRNA parameter file to be used for
	 *        parameter setup or NULL if defaults are to be used
	 */
	VrnaHandler(
			Z_type temperature = 37.0,
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
	Z_type
	getRT() const;

	/**
	 * Provides RT for the given temperature
	 * @return R*temperature
	 */
	static
	Z_type
	getRT( const Z_type temperature );



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
