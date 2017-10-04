

#include <iostream>

////////////////  CENTRAL LOGGING LIB  //////////////////

// disable default log file creation
#define ELPP_NO_DEFAULT_LOG_FILE 1
// disable log file argument (and parsing)
#define ELPP_DISABLE_LOG_FILE_FROM_ARG 1
// enable debug error tracking
#define ELPP_DEBUG_ERRORS 1
// enable performance tracking
#define ELPP_FEATURE_PERFORMANCE_TRACKING 1

#include "easylogging++.h"


// initialize logging for binary
INITIALIZE_EASYLOGGINGPP


/////////////////////////////////////////////////////////////////////

void testTiming() {
	// time logging
	TIMED_FUNC(timerObj);
//	TIMED_FUNC_IF(timerObj, VLOG_IS_ON(9));

	std::cout <<"testTiming() std::out"<<std::endl;
	std::cin.get();
}


/////////////////////////////////////////////////////////////////////
/**
 * program main entry
 *
 * @param argc number of program arguments
 * @param argv array of program arguments of length argc
 */
int main(int argc, char **argv){


	// setup default logger
	el::Configurations defaultConf;
	defaultConf.setToDefault();
	// Values are always std::string
	defaultConf.setGlobally(el::ConfigurationType::Format, "# %level %msg");
	// default logger uses default configurations
	el::Loggers::reconfigureLogger("default", defaultConf);

//		// set overall logging style
//		el::Loggers::reconfigureAllLoggers(el::ConfigurationType::Format, std::string("# %level : %msg"));
//		// TODO setup log file
//		el::Loggers::reconfigureAllLoggers(el::ConfigurationType::ToFile, std::string("false"));
//		// set additional logging flags
//		el::Loggers::addFlag(el::LoggingFlag::ColoredTerminalOutput);
//		el::Loggers::addFlag(el::LoggingFlag::DisableApplicationAbortOnFatalLog);
//		el::Loggers::addFlag(el::LoggingFlag::LogDetailedCrashReason);
//		el::Loggers::addFlag(el::LoggingFlag::AllowVerboseIfModuleNotSpecified);

	// setup logging with given parameters
	START_EASYLOGGINGPP(argc, argv);

	std::cout <<"main() std::out"<<std::endl;

	LOG(INFO) <<"main() INFO";

	testTiming();

	  // all went fine
	return 0;
}

