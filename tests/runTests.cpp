

#ifdef HAVE_CONFIG_H
	#include <config.h>
#endif

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

#include "easylogging++.h"

// initialize logging for binary
INITIALIZE_EASYLOGGINGPP

////////////////  SETUP TESTING LIB  /////////////////////

// make this the main test runner
#define CATCH_CONFIG_MAIN

#include "catch.hpp"



// don't do anything else, its done via the compiler..

