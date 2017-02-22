
#ifndef GENERAL_H_
#define GENERAL_H_


//////////////  GENERAL CONFIGURE FLAGS  ////////////////

#ifdef HAVE_CONFIG_H
	#include <config.h>
#endif

#include "intarna_config.h"

//! central compiler flag whether or not debug mode is enabled
#define IN_DEBUG_MODE ((defined(_DEBUG)) || (!defined (NDEBUG)))

////////////////  CENTRAL LOGGING LIB  //////////////////

// disable default log file creation
#ifndef ELPP_NO_DEFAULT_LOG_FILE
	#define ELPP_NO_DEFAULT_LOG_FILE 1
#endif
// disable log file argument (and parsing)
#ifndef ELPP_DISABLE_LOG_FILE_FROM_ARG
	#define ELPP_DISABLE_LOG_FILE_FROM_ARG 1
#endif
// enable debug error tracking
#if IN_DEBUG_MODE
	#define ELPP_DEBUG_ERRORS 1
#endif

#include "easylogging++.h"

////////////////  GARBAGE COLLECTION  ///////////////////

#define CLEANUP(pointer) if (pointer != NULL) {delete pointer; pointer=NULL;}

////////////////  ABORT NON-IMPLEMENTED  ////////////////

#include <stdexcept>

	#define NOTIMPLEMENTED(message)	\
		throw std::runtime_error( \
					std::string("\nSTOP DUE TO MISSING IMPLEMENTATION : ") \
					+ message \
				);

////////////////  ABORT ON NULL POINTER  ////////////////

#include <stdexcept>

	#define CHECKNOTNULL(pointer,message) if (pointer == NULL) {	\
		throw std::runtime_error( \
					std::string("\nSTOP DUE TO NULL POINTER : ") \
					+ message \
				); }


////////////////  STRING GENERATION  ////////////////////

#include <boost/lexical_cast.hpp>

	#define toString( x ) boost::lexical_cast<std::string>(x)


////////////////  NUMBER STRING LENGTH  ////////////////////

#include <cmath>

	#define numStringLength( x ) ( (x==0) ? 1 : (1+std::floor(std::log10(x))) )


////////////////  GLOBAL TYPEDEFS  //////////////////////

#include <cmath>

	//! type for energy values (energy + accessibility [ED])
	typedef float E_type;

	//! type for temperature values
	typedef E_type T_type;

	//! the delta difference range to consider two energies equivalent
#define E_precisionEpsilon 1000.0*std::numeric_limits<E_type>::epsilon()

	//! check if two energies are equal according to some epsilon
#define E_equal( e1, e2 ) ( std::abs((e1)-(e2)) < E_precisionEpsilon)

	//! check if a given energy is NOT set to E_INF
#define E_isNotINF( e ) ( std::numeric_limits<E_type>::max() >= e )

	//! check if a given energy is set to E_INF
#define E_isINF( e ) (  std::numeric_limits<E_type>::max() < e )


////////////////  GLOBAL CONSTANTS  /////////////////////

#include <limits>

	const E_type E_INF = std::numeric_limits<E_type>::infinity();


////////////////  UTILITY FUNCTION  /////////////////////

#include <iostream>
#include <string>

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



#endif /* GENERAL_H_ */
