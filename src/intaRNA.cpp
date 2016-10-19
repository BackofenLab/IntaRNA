
#include "general.h"

#include <iostream>

#include "CommandLineParsing.h"

#include "RnaSequence.h"

#include "Accessibility.h"

#include "InteractionEnergy.h"

#include "Predictor.h"

#include "OutputHandler.h"


// initialize logging for binary
INITIALIZE_EASYLOGGINGPP

/////////////////////////////////////////////////////////////////////
/**
 * program main entry
 *
 * @param argc number of program arguments
 * @param argv array of program arguments of length argc
 */
int main(int argc, char **argv) {

	try {

		// setup logging with given parameters
		START_EASYLOGGINGPP(argc, argv);
		// set overall logging style
		el::Loggers::reconfigureAllLoggers(el::ConfigurationType::Format, std::string("# %level : %msg"));
		// TODO setup log file
		el::Loggers::reconfigureAllLoggers(el::ConfigurationType::ToFile, std::string("false"));
		// set additional logging flags
		el::Loggers::addFlag(el::LoggingFlag::ColoredTerminalOutput);
		el::Loggers::addFlag(el::LoggingFlag::DisableApplicationAbortOnFatalLog);
		el::Loggers::addFlag(el::LoggingFlag::LogDetailedCrashReason);
		el::Loggers::addFlag(el::LoggingFlag::AllowVerboseIfModuleNotSpecified);

		// parse command line parameters
		CommandLineParsing parameters;
		{
			VLOG(1) <<"parsing arguments"<<"...";
			int retCode = parameters.parse( argc, argv );
			if (retCode != CommandLineParsing::ReturnCode::KEEP_GOING) {
				return retCode;
			}
		}


		// run prediction for all pairs of sequences
		// first: iterate over all query sequences
		for ( size_t queryNumber = 0; queryNumber < parameters.getQuerySequences().size(); ++queryNumber )
		{

			// get accessibility handler
			VLOG(1) <<"computing accessibility for query "<<parameters.getQuerySequences().at(queryNumber).getId()<<"...";
			Accessibility * queryAcc = parameters.getQueryAccessibility(queryNumber);
			CHECKNOTNULL(queryAcc,"query initialization failed");

			// check if we have to warn about ambiguity
			if (queryAcc->getSequence().isAmbiguous()) {
				LOG(INFO) <<"Sequence '"<<queryAcc->getSequence().getId()
						<<"' contains ambiguous nucleotide encodings. These positions are ignored for interaction computation.";
			}

			// second: iterate over all target sequences to get all pairs to predict for
			for ( size_t targetNumber = 0; targetNumber < parameters.getTargetSequences().size(); ++targetNumber )
			{

				// get target accessibility handler
				VLOG(1) <<"computing accessibility for target "<<parameters.getTargetSequences().at(targetNumber).getId()<<"...";
				Accessibility * targetAccOrig = parameters.getTargetAccessibility(targetNumber);
				CHECKNOTNULL(targetAccOrig,"target initialization failed");
				// reverse indexing of target sequence for the computation
				ReverseAccessibility * targetAcc = new ReverseAccessibility(*targetAccOrig);

				// check if we have to warn about ambiguity
				if (targetAcc->getSequence().isAmbiguous()) {
					LOG(INFO) <<"Sequence '"<<targetAcc->getSequence().getId()
							<<"' contains ambiguous nucleotide encodings. These positions are ignored for interaction computation.";
				}

				// get energy computation handler for both sequences
				InteractionEnergy* energy = parameters.getEnergyHandler( *queryAcc, *targetAcc );
				CHECKNOTNULL(energy,"energy initialization failed");

				// get output/storage handler
				OutputHandler * output = parameters.getOutputHandler( *energy );
				CHECKNOTNULL(output,"output handler initialization failed");

				// get interaction prediction handler
				Predictor * predictor = parameters.getPredictor( *energy, *output );
				CHECKNOTNULL(predictor,"predictor initialization failed");

				// run prediction
				// TODO we can also limit the prediction range eg. for streamed predictions
				VLOG(1) <<"predicting interactions...";
				predictor->predict();


				// garbage collection
				CLEANUP(predictor)
				CLEANUP(output)
				CLEANUP(energy)
				CLEANUP(targetAcc)
				CLEANUP(targetAccOrig)
			}
			// garbage collection
			CLEANUP(queryAcc)

		}



	////////////////////// exception handling ///////////////////////////
	} catch (std::exception & e) {
		LOG(ERROR) <<"Exception raised : " <<e.what() <<"\n";
		return -1;
	}

	  // all went fine
	return 0;
}

