
#include "general.h"

#include <iostream>
#include <exception>

#include <omp.h>

#include <boost/foreach.hpp>

#include "CommandLineParsing.h"

#include "RnaSequence.h"
#include "Accessibility.h"
#include "InteractionEnergy.h"
#include "Predictor.h"
#include "OutputHandler.h"
#include "OutputHandlerIntaRNA1.h"

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

		// number of already reported interactions to enable IntaRNA v1 separator output
		size_t reportedInteractions = 0;

		// storage to avoid accessibility recomputation (init NULL)
		std::vector< ReverseAccessibility * > queryAcc(parameters.getQuerySequences().size(), NULL);

		// compute all query accessibilities to enable parallelization
		// do serially since not all VRNA routines are threadsafe
		for (size_t qi=0; qi<queryAcc.size(); qi++) {
			// get accessibility handler
			VLOG(1) <<"computing accessibility for query '"<<parameters.getQuerySequences().at(qi).getId()<<"'...";
			Accessibility * queryAccOrig = parameters.getQueryAccessibility(qi);
			CHECKNOTNULL(queryAccOrig,"query initialization failed");
			// reverse indexing of target sequence for the computation
			queryAcc[qi] = new ReverseAccessibility(*queryAccOrig);

			// check if we have to warn about ambiguity
			if (queryAccOrig->getSequence().isAmbiguous()) {
				LOG(INFO) <<"Sequence '"<<queryAccOrig->getSequence().getId()
						<<"' contains ambiguous nucleotide encodings. These positions are ignored for interaction computation.";
			}
		}

		// check which loop to parallelize
		const bool parallelizeTargetLoop = parameters.getTargetSequences().size() > 1;
		const bool parallelizeQueryLoop = !parallelizeTargetLoop && parameters.getQuerySequences().size() > 1;

		// run prediction for all pairs of sequences
		// first: iterate over all target sequences
		// parallelize this loop if possible; if not -> parallelize the query-loop
		# pragma omp parallel for schedule(dynamic) num_threads( parameters.getThreads() ) shared(queryAcc,reportedInteractions) if(parallelizeTargetLoop)
		for ( size_t targetNumber = 0; targetNumber < parameters.getTargetSequences().size(); ++targetNumber )
		{
			// explicit try-catch-block due to missing OMP exception forwarding
			try {
				// get target accessibility handler
				VLOG(1) <<"computing accessibility for target '"<<parameters.getTargetSequences().at(targetNumber).getId()<<"'...";

				// VRNA not completely threadsafe ...
				Accessibility * targetAcc = parameters.getTargetAccessibility(targetNumber);
				CHECKNOTNULL(targetAcc,"target initialization failed");

				// check if we have to warn about ambiguity
				if (targetAcc->getSequence().isAmbiguous()) {
					LOG(INFO) <<"Sequence '"<<targetAcc->getSequence().getId()
							<<"' contains ambiguous IUPAC nucleotide encodings. These positions are ignored for interaction computation and replaced by 'N'.";
				}

				// second: iterate over all query sequences
				// this parallelization should only be enabled if the outer target-loop is not parallelized
				# pragma omp parallel for schedule(dynamic) num_threads( parameters.getThreads() ) shared(queryAcc,reportedInteractions,targetAcc,targetNumber) if(parallelizeQueryLoop)
				for ( size_t queryNumber = 0; queryNumber < parameters.getQuerySequences().size(); ++queryNumber )
				{
					// explicit try-catch-block due to missing OMP exception forwarding
					try {
						// sanity check
						assert( queryAcc.at(queryNumber) != NULL );

						// get energy computation handler for both sequences
						InteractionEnergy* energy = parameters.getEnergyHandler( *targetAcc, *(queryAcc.at(queryNumber)) );
						CHECKNOTNULL(energy,"energy initialization failed");

						// get output/storage handler
						OutputHandler * output = parameters.getOutputHandler( *energy );
						CHECKNOTNULL(output,"output handler initialization failed");

						// check if we have to add separator for IntaRNA v1 output
						if (reportedInteractions > 0 && dynamic_cast<OutputHandlerIntaRNA1*>(output) != NULL) {
							dynamic_cast<OutputHandlerIntaRNA1*>(output)->addSeparator( true );
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
												, queryAcc.at(queryNumber)->getReversedIndexRange(qRange)
												, parameters.getOutputConstraint()
												);

						} // target ranges
						} // query ranges

						#pragma omp atomic update
						reportedInteractions += output->reported();

						// garbage collection
						CLEANUP(predictor);
						CLEANUP(output);
						CLEANUP(energy);

					////////////////////// exception handling ///////////////////////////
					} catch (std::exception & e) {
						LOG(DEBUG) <<"Exception raised for  target "
									<<parameters.getTargetSequences().at(targetNumber).getId()
									<<" and query "
									<<parameters.getQuerySequences().at(queryNumber).getId()
									<<" : " <<e.what() <<"\n\n"
				#if IN_DEBUG_MODE
							<<"  ==> run debugger for details\n";
						throw e;
				#else
							<<"  ==> Please report to the IntaRNA development team! Thanks!\n";
						return -1;
				#endif
					} catch (...) {
						std::exception_ptr eptr = std::current_exception();
						LOG(DEBUG) <<"Unknown exception raised for target "
									<<parameters.getTargetSequences().at(targetNumber).getId()
									<<" and query "
									<<parameters.getQuerySequences().at(queryNumber).getId()
									<<"\n\n"
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
				}

				// write accessibility to file if needed
				parameters.writeTargetAccessibility( *targetAcc );

				// garbage collection
				CLEANUP(targetAcc);

			////////////////////// exception handling ///////////////////////////
			} catch (std::exception & e) {
				LOG(DEBUG) <<"Exception raised for target "
							<<parameters.getTargetSequences().at(targetNumber).getId()
							<<" : " <<e.what() <<"\n\n"
		#if IN_DEBUG_MODE
					<<"  ==> run debugger for details\n";
				throw e;
		#else
					<<"  ==> Please report to the IntaRNA development team! Thanks!\n";
				return -1;
		#endif
			} catch (...) {
				std::exception_ptr eptr = std::current_exception();
				LOG(DEBUG) <<"Unknown exception raised for target "
							<<parameters.getTargetSequences().at(targetNumber).getId()
							<<"\n\n"
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
		}
		// garbage collection
		for (size_t queryNumber=0; queryNumber < queryAcc.size(); queryNumber++) {
			// this is a hack to cleanup the original accessibility object
			Accessibility* queryAccOrig = &(const_cast<Accessibility&>(queryAcc[queryNumber]->getAccessibilityOrigin()) );
			// write accessibility to file if needed
			parameters.writeQueryAccessibility( *queryAccOrig );
			CLEANUP( queryAccOrig );
			// cleanup (now broken) reverse accessibility object
			CLEANUP(queryAcc[queryNumber]);
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

