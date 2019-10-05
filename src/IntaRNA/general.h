
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

template<class T> void INTARNA_CLEANUP( T *& pointer) { if (pointer != NULL) {delete pointer; pointer=NULL;} }

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


////////////////  GLOBAL TYPEDEFS  //////////////////////

#include <cmath>
#include <limits>

#if INTARNA_MULTIPRECISION
	#include <boost/multiprecision/float128.hpp>
#endif

namespace IntaRNA {

	//! type for energy values in kcal/mol units for in-/output only
	typedef float E_kcal_type;

	//! type for energy values (energy + accessibility [ED]) (internally)
	typedef int E_type;
	const E_type E_INF = (std::numeric_limits<E_type>::max() / 8) + 1;
	const E_type E_MAX = E_INF / 2;

	//! type for probabilities, RT and Boltzmann values
#if INTARNA_MULTIPRECISION
	typedef boost::multiprecision::float128 Z_type;
	#define Z_log boost::multiprecision::log
	#define Z_exp boost::multiprecision::exp
#else
	typedef double Z_type;
	#define Z_log std::log
	#define Z_exp std::exp
#endif
	const Z_type Z_INF = std::numeric_limits<Z_type>::infinity();

} // namespace

#ifdef IntaRNA_precisionEpsilon
	#error IntaRNA_precisionEpsilon already defined
#endif
	//! the delta difference range to consider two floating point values equivalent
#define IntaRNA_precisionEpsilon std::pow(10, -10)

#ifdef E_2_Ekcal
  #error E_2_Ekcal already defined
#endif
  //! convert internal energy type to energy value in kcal/mol units
#define E_2_Ekcal( e ) ( static_cast<E_kcal_type>(e) / 100.0 )

#ifdef Ekcal_2_E
  #error Ekcal_2_E already defined
#endif
  //! convert energy in kcal/mol units to internal energy type
#define Ekcal_2_E( e ) ( static_cast<E_type>(e * 100) )

#ifdef E_2_Z
  #error E_2_Z already defined
#endif
  //! convert E_type to Z_type
#define E_2_Z( e ) ( static_cast<Z_type>(e) / 100.0 )

#ifdef Z_2_E
  #error Z_2_E already defined
#endif
  //! convert Z_type to E_type
#define Z_2_E( e ) ( static_cast<E_type>(e * 100) )

#ifdef E_equal
	#error E_equal already defined
#endif
	//! check if two energies are equal
#define E_equal( e1, e2 ) ( e1 == e2 )

#ifdef E_isNotINF
	#error E_isNotINF already defined
#endif
	//! check if a given energy is NOT set to E_INF
#define E_isNotINF( e ) ( E_INF > e )

#ifdef E_isINF
	#error E_isINF already defined
#endif
	//! check if a given energy is set to E_INF
#define E_isINF( e ) ( E_INF <= e )



#ifdef Z_equal
	#error Z_equal already defined
#endif
	//! check if two energies are equal according to some epsilon
#if INTARNA_MULTIPRECISION
	#define Z_equal( e1, e2 ) ( boost::multiprecision::abs((e1)-(e2)) < IntaRNA_precisionEpsilon)
#else
	#define Z_equal( e1, e2 ) ( std::abs((e1)-(e2)) < IntaRNA_precisionEpsilon)
#endif
// another option from http://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
//#define Z_equal_ULP 2
//#define Z_equal( e1, e2 ) ( \
//	/* the machine epsilon has to be scaled to the magnitude of the values used */ \
//	/* and multiplied by the desired precision in ULPs (units in the last place) */ \
//	std::abs(e1-e2) < std::numeric_limits<T>::epsilon() * std::abs(e1+e2) * Z_equal_ULP \
//	/* unless the result is subnormal */ \
//	|| std::abs(e1-e2) < std::numeric_limits<T>::min() \
//)

#ifdef Z_isNotINF
	#error Z_isNotINF already defined
#endif
	//! check if a given energy is NOT set to Z_INF
#define Z_isNotINF( e ) ( std::numeric_limits<Z_type>::max() >= e )

#ifdef Z_isINF
	#error Z_isINF already defined
#endif
	//! check if a given energy is set to Z_INF
#define Z_isINF( e ) (  std::numeric_limits<Z_type>::max() < e )


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
 * If the filename ends in '.gz', gzip compression and binary output is enabled.
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
 * Deletes an output stream pointer returned by newOutputStream() and sets it to NULL.
 * If it is a file stream, it is closed before deletion.
 * If it is pointing to std::cout or std::cerr, only the pointer is reset to NULL.
 *
 * @param outStream the pointer to the output stream to delete; will be set to NULL
 */
void
deleteOutputStream( std::ostream *& outStream );

/**
 * Returns an open input stream for the given input name, i.e.
 *
 * - & std::cin : if inName == STDIN
 * - new std::fstream( inName ) : else inName non-empty
 *
 * If the filename ends in '.gz', gzip compression and binary input is enabled.
 *
 * @param inName the name of the input to open. use STDIN for the
 *        respective input stream or otherwise a filename to be created.
 *
 * @return the open input stream, or NULL in error case
 *
 */
std::istream *
newInputStream( const std::string & inName );


/**
 * Deletes an input stream pointer returned by newInputStream() and sets it to NULL.
 * If it is a file stream, it is closed before deletion.
 * If it is pointing to std::cin, only the pointer is reset to NULL.
 *
 * @param inStream the pointer to the input stream to delete; will be set to NULL
 */
void
deleteInputStream( std::istream *& inStream );

} // namespace

#endif /* GENERAL_H_ */
