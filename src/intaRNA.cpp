
#include "general.h"

#include <iostream>
#include <exception>

#include <boost/foreach.hpp>

#include "CommandLineParsing.h"

#include "RnaSequence.h"
#include "Accessibility.h"
#include "InteractionEnergy.h"
#include "Predictor.h"
#include "OutputHandler.h"
#include "OutputHandlerIntaRNA1detailed.h"

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
	
		// set overall logging style
		el::Loggers::reconfigureAllLoggers(el::ConfigurationType::Format, std::string("# %level : %msg"));
		// TODO setup log file
		el::Loggers::reconfigureAllLoggers(el::ConfigurationType::ToFile, std::string("false"));
		// set additional logging flags
		el::Loggers::addFlag(el::LoggingFlag::ColoredTerminalOutput);
		el::Loggers::addFlag(el::LoggingFlag::DisableApplicationAbortOnFatalLog);
		el::Loggers::addFlag(el::LoggingFlag::LogDetailedCrashReason);
		el::Loggers::addFlag(el::LoggingFlag::AllowVerboseIfModuleNotSpecified);

		// setup logging with given parameters
		START_EASYLOGGINGPP(argc, argv);

		// parse command line parameters
		CommandLineParsing parameters;
		{
			VLOG(1) <<"parsing arguments"<<"...";
			int retCode = parameters.parse( argc, argv );
			if (retCode != CommandLineParsing::ReturnCode::KEEP_GOING) {
				return retCode;
			}
		}

		size_t reportedInteractions = 0;

		// second: iterate over all target sequences to get all pairs to predict for
		for ( size_t targetNumber = 0; targetNumber < parameters.getTargetSequences().size(); ++targetNumber )
		{

			// get target accessibility handler
			VLOG(1) <<"computing accessibility for target '"<<parameters.getTargetSequences().at(targetNumber).getId()<<"'...";
			Accessibility * targetAcc = parameters.getTargetAccessibility(targetNumber);
			CHECKNOTNULL(targetAcc,"target initialization failed");

			// check if we have to warn about ambiguity
			if (targetAcc->getSequence().isAmbiguous()) {
				LOG(INFO) <<"Sequence '"<<targetAcc->getSequence().getId()
						<<"' contains ambiguous IUPAC nucleotide encodings. These positions are ignored for interaction computation and replaced by 'N'.";
			}

			// run prediction for all pairs of sequences
			// first: iterate over all query sequences
			// TODO maybe parallelize for heuristic mode (low mem per job); BUT ENSURE accessibility computation via plfold only (rnaup not threadsafe)
			for ( size_t queryNumber = 0; queryNumber < parameters.getQuerySequences().size(); ++queryNumber )
			{

				// get accessibility handler
				VLOG(1) <<"computing accessibility for query '"<<parameters.getQuerySequences().at(queryNumber).getId()<<"'...";
				Accessibility * queryAccOrig = parameters.getQueryAccessibility(queryNumber);
				CHECKNOTNULL(queryAccOrig,"query initialization failed");
				// reverse indexing of target sequence for the computation
				ReverseAccessibility * queryAcc = new ReverseAccessibility(*queryAccOrig);

				// check if we have to warn about ambiguity
				if (queryAcc->getSequence().isAmbiguous()) {
					LOG(INFO) <<"Sequence '"<<queryAcc->getSequence().getId()
							<<"' contains ambiguous nucleotide encodings. These positions are ignored for interaction computation.";
				}


				// get energy computation handler for both sequences
				InteractionEnergy* energy = parameters.getEnergyHandler( *targetAcc, *queryAcc );
				CHECKNOTNULL(energy,"energy initialization failed");

				// get output/storage handler
				OutputHandler * output = parameters.getOutputHandler( *energy );
				CHECKNOTNULL(output,"output handler initialization failed");

				// check if we have to add separator for IntaRNA v1 output
				if (reportedInteractions > 0 && dynamic_cast<OutputHandlerIntaRNA1detailed*>(output) != NULL) {
					dynamic_cast<OutputHandlerIntaRNA1detailed*>(output)->addSeparator( true );
				}

				// get interaction prediction handler
				Predictor * predictor = parameters.getPredictor( *energy, *output );
				CHECKNOTNULL(predictor,"predictor initialization failed");

				// run prediction for all range combinations
				BOOST_FOREACH(const IndexRange & tRange, parameters.getTargetRanges(targetNumber)) {
				BOOST_FOREACH(const IndexRange & qRange, parameters.getQueryRanges(queryNumber)) {

					VLOG(1) <<"predicting interactions for"
							<<" target range " <<tRange
							<<" and"
							<<" query range " <<qRange
							<<"...";

					predictor->predict(	  tRange
										, queryAcc->getReversedIndexRange(qRange)
										, parameters.getOutputConstraint()
										);

				} // target ranges
				} // query ranges

				reportedInteractions += output->reported();

				// garbage collection
				CLEANUP(predictor)
				CLEANUP(output)
				CLEANUP(energy)
				CLEANUP(queryAcc)
				CLEANUP(queryAccOrig)
			}
			// garbage collection
			CLEANUP(targetAcc)

		}



	////////////////////// exception handling ///////////////////////////
	} catch (std::exception & e) {
		LOG(DEBUG) <<"Exception raised : " <<e.what() <<"\n\n"
#if IN_DEBUG_MODE
			<<"  ==> run debugger for details\n";
		throw e;
#else
			<<"  ==> Please report to the IntaRNA development team! Thanks!\n";
		return -1;
#endif
	} catch (...) {
		std::exception_ptr eptr = std::current_exception();
		LOG(DEBUG) <<"Unknown exception raised \n\n"
#if IN_DEBUG_MODE
			<<"  ==> run debugger for details\n";
		if (eptr) {
			std::rethrow_exception(eptr);
		}
#else
			<<"  ==> Please report to the IntaRNA development team! Thanks!\n";
		return -1;
#endif

	}

	  // all went fine
	return 0;
}

