/*
 * ViennaSetup.cpp
 *
 *  Created on: 22.09.2016
 *      Author: Mmann
 */

#include "VrnaHandler.h"
#include "general.h"

#include <iostream>

extern "C" {
//	#include <ViennaRNA/model.h>
}



////////////////////////////////////////////////////////////////////////////

VrnaHandler::
VrnaHandler( double temperature, int max_bp_span, int window_size )
	:
	temperature(temperature),
	max_bp_span(max_bp_span),
	window_size(window_size)
{


}

////////////////////////////////////////////////////////////////////////////

VrnaHandler::
~VrnaHandler() {
}

////////////////////////////////////////////////////////////////////////////




vrna_md_t
VrnaHandler::
getModel( ) const
{

	std::cerr <<"\n## VRNAHANDLER : GETMODEL\n"<<std::endl;
	// TODO overwrite locally
	vrna_md_t modelDetails;

	// get default model data
	vrna_md_set_default( &modelDetails );

	// set model details according to command line arguments

	modelDetails.temperature = this->temperature;
//	  double  temperature;                  /**<  @brief  The temperature used to scale the thermodynamic parameters */
//	  double  betaScale;                    /**<  @brief  A scaling factor for the thermodynamic temperature of the Boltzmann factors */
	modelDetails.dangles = 2;
//	  int     dangles;                      /**<  @brief  Specifies the dangle model used in any energy evaluation (0,1,2 or 3) */
//	  int     special_hp;                   /**<  @brief  Include special hairpin contributions for tri, tetra and hexaloops */
	modelDetails.noLP = 0;
//	  int     noLP;                         /**<  @brief  Only consider canonical structures, i.e. no 'lonely' base pairs */
	modelDetails.noGU = 0;
//	  int     noGU;                         /**<  @brief  Do not allow GU pairs */
	modelDetails.noGUclosure = 0;
//	  int     noGUclosure;                  /**<  @brief  Do not allow loops to be closed by GU pair */
//	  int     logML;                        /**<  @brief  Use logarithmic scaling for multi loops */
	modelDetails.circ = 0;
//	  int     circ;                         /**<  @brief  Assume RNA to be circular instead of linear */
	modelDetails.gquad = 0;
//	  int     gquad;                        /**<  @brief  Include G-quadruplexes in structure prediction */
//	  int     canonicalBPonly;              /**<  @brief  remove non-canonical bp's from constraint structures  */
//	  int     uniq_ML;                      /**<  @brief  Flag to ensure unique multibranch loop decomposition during folding */
//	  int     energy_set;                   /**<  @brief  Specifies the energy set that defines set of compatible base pairs */
//	  int     backtrack;                    /**<  @brief  Specifies whether or not secondary structures should be backtraced */
//	  char    backtrack_type;               /**<  @brief  Specifies in which matrix to backtrack */
//	  int     compute_bpp;                  /**<  @brief  Specifies whether or not backward recursions for base pair probability (bpp) computation will be performed */
//	  char    nonstandards[64];             /**<  @brief  contains allowed non standard bases */
	modelDetails.max_bp_span = this->max_bp_span;
//	  int     max_bp_span;                  /**<  @brief  maximum allowed base pair span */
//	  int     min_loop_size;                /**<  @brief  Minimum size of hairpin loops */
	modelDetails.window_size = this->window_size;
//	  int     window_size;                  /**<  @brief  Size of the sliding window for locally optimal structure predition */
//	  int     oldAliEn;                     /**<  @brief  Use old alifold energy model */
//	  int     ribo;                         /**<  @brief  Use ribosum scoring table in alifold energy model */
//	  double  cv_fact;                      /**<  @brief  Covariance scaling factor for consensus structure prediction */
//	  double  nc_fact;                      /**<  @brief  Scaling factor to weight covariance contributions of non-canonical pairs */
//	  double  sfact;                        /**<  @brief  Scaling factor for partition function scaling */
//	  int     rtype[8];                     /**<  @brief  Reverse base pair type array */
//	  short   alias[MAXALPHA+1];            /**<  @brief  alias of an integer nucleotide representation */
//	  int     pair[MAXALPHA+1][MAXALPHA+1]; /**<  @brief  Integer representation of a base pair */


	return modelDetails;


}

////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////


