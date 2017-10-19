
#ifndef INTARNA_GENERAL_H_
#define INTARNA_GENERAL_H_


//////////////  GENERAL CONFIGURE FLAGS  ////////////////

#ifdef HAVE_CONFIG_H
	#include "config.h"
#endif

#include "IntaRNA/intarna_config.h"

//! central flag whether or not debug mode is enabled
#define INTARNA_IN_DEBUG_MODE ((defined(_DEBUG)) || (!defined (NDEBUG)))

////////////////  CENTRAL LOGGING LIB  //////////////////

#include "easylogging++.h"

////////////////  GARBAGE COLLECTION  ///////////////////

#define  INTARNA_CLEANUP(pointer) if (pointer != NULL) {delete pointer; pointer=NULL;}

////////////////  ABORT NON-IMPLEMENTED  ////////////////

#include <stdexcept>

	#define INTARNA_NOT_IMPLEMENTED(message)	\
		throw std::runtime_error( \
					std::string("\nSTOP DUE TO MISSING IMPLEMENTATION : ") \
					+ message \
				);

////////////////  ABORT ON NULL POINTER  ////////////////

#include <stdexcept>

	#define INTARNA_CHECK_NOT_NULL(pointer,message) if (pointer == NULL) {	\
		throw std::runtime_error( \
					std::string("\nSTOP DUE TO NULL POINTER : ") \
					+ message \
				); }


////////////////  STRING GENERATION  ////////////////////

#include <boost/lexical_cast.hpp>

	#ifdef toString
		#error toString already defined
	#endif
	#define toString( x ) boost::lexical_cast<std::string>(x)


////////////////  NUMBER STRING LENGTH  ////////////////////

#include <cmath>

	#ifdef numStringLength
		#error numStringLength already defined
	#endif
	#define numStringLength( x ) ( (x==0) ? 1 : (1+std::floor(std::log10(x))) )


////////////////  GLOBAL TYPEDEFS  //////////////////////

#include <cmath>

namespace IntaRNA {

	//! type for energy values (energy + accessibility [ED])
	typedef float E_type;

	//! type for temperature values
	typedef E_type T_type;

} // namespace

#ifdef E_precisionEpsilon
	#error E_precisionEpsilon already defined
#endif
	//! the delta difference range to consider two energies equivalent
#define E_precisionEpsilon 1000.0*std::numeric_limits<E_type>::epsilon()

#ifdef E_equal
	#error E_equal already defined
#endif
	//! check if two energies are equal according to some epsilon
#define E_equal( e1, e2 ) ( std::abs((e1)-(e2)) < E_precisionEpsilon)

#ifdef E_isNotINF
	#error E_isNotINF already defined
#endif
	//! check if a given energy is NOT set to E_INF
#define E_isNotINF( e ) ( std::numeric_limits<E_type>::max() >= e )

#ifdef E_isINF
	#error E_isINF already defined
#endif
	//! check if a given energy is set to E_INF
#define E_isINF( e ) (  std::numeric_limits<E_type>::max() < e )


////////////////  GLOBAL CONSTANTS  /////////////////////

#include <limits>

namespace IntaRNA {

	const E_type E_INF = std::numeric_limits<E_type>::infinity();

} // namespace


////////////////  UTILITY FUNCTION  /////////////////////

#include <iostream>
#include <string>

namespace IntaRNA {

/**
 * Returns an open output stream for the given output name, i.e.
 *
 * - & std::cout : if outName == STDOUT
 * - & std::cerr : if outName == STDERR
 * - new std::fstream( outName ) : else if outName non-empty
 *
 * @param outName the name of the output to open. use STDOUT/STDERR for the
 *        respective output stream or otherwise a filename to be created.
 *
 * @return the open output stream, or NULL in error case
 *
 */
std::ostream *
newOutputStream( const std::string & outName );


/**
 * Deletes an output stream pointer returned by newOutputStream().
 * If it is a file stream, it is closed before deletion.
 * If it is pointing to std::cout or std::cerr, nothing is done.
 *
 * @param outStream the pointer to the output stream to delete
 */
void
deleteOutputStream( std::ostream * outStream );

} // namespace

#endif /* GENERAL_H_ */
