/*
 * general.h
 *
 *  Created on: 24.06.2014
 *      Author: Mmann
 */

#ifndef GENERAL_H_
#define GENERAL_H_


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

#define E_equal( e1, e2 ) ( std::abs((e1)-(e2)) < 1000.0*std::numeric_limits<E_type>::epsilon())


////////////////  GLOBAL CONSTANTS  /////////////////////

#include <limits>

	const E_type E_INF = std::numeric_limits<E_type>::infinity();

#endif /* GENERAL_H_ */
