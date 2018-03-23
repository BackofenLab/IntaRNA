
#include "IntaRNA/general.h"

#include <iostream>
#include <exception>

#if INTARNA_MULITHREADING
	#include <omp.h>
#endif

#include <boost/foreach.hpp>

#include "CommandLineParsing.h"

#include "IntaRNA/RnaSequence.h"
#include "IntaRNA/Accessibility.h"
#include "IntaRNA/InteractionEnergy.h"
#include "IntaRNA/Predictor.h"
#include "IntaRNA/OutputHandler.h"
#include "IntaRNA/OutputHandlerIntaRNA1.h"
#include "IntaRNA/OutputHandlerInteractionList.h"

// initialize logging for binary
INITIALIZE_EASYLOGGINGPP

using namespace IntaRNA;

/////////////////////////////////////////////////////////////////////
/**
 * program main entry
 *
 * @param argc number of program arguments
 * @param argv array of program arguments of length argc
 */
int main(int argc, char **argv){

	try {

		// set overall logging style
		el::Loggers::reconfigureAllLoggers(el::ConfigurationType::Format, std::string("# %level : %msg"));
		// TODO setup log file
		el::Loggers::reconfigureAllLoggers(el::ConfigurationType::ToFile, std::string("false"));
		el::Loggers::reconfigureAllLoggers(el::ConfigurationType::ToStandardOutput, std::string("true"));
		// set additional logging flags
		el::Loggers::addFlag(el::LoggingFlag::DisableApplicationAbortOnFatalLog);
		el::Loggers::addFlag(el::LoggingFlag::LogDetailedCrashReason);
		el::Loggers::addFlag(el::LoggingFlag::AllowVerboseIfModuleNotSpecified);
#if INTARNA_LOG_COLORING
		el::Loggers::addFlag(el::LoggingFlag::ColoredTerminalOutput);
#endif

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

#if INTARNA_MULITHREADING
		// OMP shared variables to enable exception forwarding from within OMP parallelized for loop
		bool threadAborted = false;
		std::exception_ptr exceptionPtrDuringOmp = NULL;
		std::stringstream exceptionInfoDuringOmp;
#endif


		// number of already reported interactions to enable IntaRNA v1 separator output
		size_t reportedInteractions = 0;

		// storage to avoid accessibility recomputation (init NULL)
		std::vector< ReverseAccessibility * > queryAcc(parameters.getQuerySequences().size(), NULL);

		// compute all query accessibilities to enable parallelization
#if INTARNA_MULITHREADING
		// parallelize this loop if possible; if not -> parallelize the query-loop
		# pragma omp parallel for schedule(dynamic) num_threads( parameters.getThreads() ) shared(queryAcc,reportedInteractions,exceptionPtrDuringOmp,exceptionInfoDuringOmp)
#endif
		for (size_t qi=0; qi<queryAcc.size(); qi++) {
			// get accessibility handler
#if INTARNA_MULITHREADING
			#pragma omp flush (threadAborted)
			// explicit try-catch-block due to missing OMP exception forwarding
			if (!threadAborted) {
				try {
					// get query accessibility handler
					#pragma omp critical(intarna_omp_logOutput)
#endif
					VLOG(1) <<"computing accessibility for query '"<<parameters.getQuerySequences().at(qi).getId()<<"'...";
					Accessibility * queryAccOrig = parameters.getQueryAccessibility(qi);
					INTARNA_CHECK_NOT_NULL(queryAccOrig,"query initialization failed");
					// reverse indexing of target sequence for the computation
					queryAcc[qi] = new ReverseAccessibility(*queryAccOrig);

					// check if we have to warn about ambiguity
					if (queryAccOrig->getSequence().isAmbiguous()) {
#if INTARNA_MULITHREADING
						#pragma omp critical(intarna_omp_logOutput)
#endif
						LOG(INFO) <<"Sequence '"<<queryAccOrig->getSequence().getId()
								<<"' contains ambiguous nucleotide encodings. These positions are ignored for interaction computation.";
					}
#if INTARNA_MULITHREADING
				////////////////////// exception handling ///////////////////////////
				} catch (std::exception & e) {
					// ensure exception handling for first failed thread only
					#pragma omp critical(intarna_omp_exception)
					{
						if (!threadAborted) {
							// store exception information
							exceptionPtrDuringOmp = std::make_exception_ptr(e);
							exceptionInfoDuringOmp <<" #thread "<<omp_get_thread_num() <<" #query "<<qi <<" : "<<e.what();
							// trigger abortion of all threads
							threadAborted = true;
							#pragma omp flush (threadAborted)
						}
					} // omp critical(intarna_omp_exception)
				} catch (...) {
					// ensure exception handling for first failed thread only
					#pragma omp critical(intarna_omp_exception)
					{
						if (!threadAborted) {
							// store exception information
							exceptionPtrDuringOmp = std::current_exception();
							exceptionInfoDuringOmp <<" #thread "<<omp_get_thread_num() <<" #query "<<qi;
							// trigger abortion of all threads
							threadAborted = true;
							#pragma omp flush (threadAborted)
						}
					} // omp critical(intarna_omp_exception)
				}
			} // if not threadAborted
#endif
		}

		// check which loop to parallelize
		const bool parallelizeTargetLoop = parameters.getTargetSequences().size() > 1;
		const bool parallelizeQueryLoop = !parallelizeTargetLoop && parameters.getQuerySequences().size() > 1;


		// run prediction for all pairs of sequences
		// first: iterate over all target sequences
#if INTARNA_MULITHREADING
		// parallelize this loop if possible; if not -> parallelize the query-loop
		# pragma omp parallel for schedule(dynamic) num_threads( parameters.getThreads() ) shared(queryAcc,reportedInteractions,exceptionPtrDuringOmp,exceptionInfoDuringOmp) if(parallelizeTargetLoop)
#endif
		for ( size_t targetNumber = 0; targetNumber < parameters.getTargetSequences().size(); ++targetNumber )
		{
#if INTARNA_MULITHREADING
			#pragma omp flush (threadAborted)
			// explicit try-catch-block due to missing OMP exception forwarding
			if (!threadAborted) {
				try {
					// get target accessibility handler
					#pragma omp critical(intarna_omp_logOutput)
#endif
					{ VLOG(1) <<"computing accessibility for target '"<<parameters.getTargetSequences().at(targetNumber).getId()<<"'..."; }

					// VRNA not completely threadsafe ...
					Accessibility * targetAcc = parameters.getTargetAccessibility(targetNumber);
					INTARNA_CHECK_NOT_NULL(targetAcc,"target initialization failed");

					// check if we have to warn about ambiguity
					if (targetAcc->getSequence().isAmbiguous()) {
#if INTARNA_MULITHREADING
						#pragma omp critical(intarna_omp_logOutput)
#endif
						{ LOG(INFO) <<"Sequence '"<<targetAcc->getSequence().getId()
								<<"' contains ambiguous IUPAC nucleotide encodings. These positions are ignored for interaction computation and replaced by 'N'.";}
					}

					// second: iterate over all query sequences
#if INTARNA_MULITHREADING
					// this parallelization should only be enabled if the outer target-loop is not parallelized
					# pragma omp parallel for schedule(dynamic) num_threads( parameters.getThreads() ) shared(queryAcc,reportedInteractions,exceptionPtrDuringOmp,exceptionInfoDuringOmp,targetAcc,targetNumber) if(parallelizeQueryLoop)
#endif
					for ( size_t queryNumber = 0; queryNumber < parameters.getQuerySequences().size(); ++queryNumber )
					{
#if INTARNA_MULITHREADING
						#pragma omp flush (threadAborted)
						// explicit try-catch-block due to missing OMP exception forwarding
						if (!threadAborted) {
							try {
#endif
								// sanity check
								assert( queryAcc.at(queryNumber) != NULL );

								// get energy computation handler for both sequences
								InteractionEnergy* energy = parameters.getEnergyHandler( *targetAcc, *(queryAcc.at(queryNumber)) );
								INTARNA_CHECK_NOT_NULL(energy,"energy initialization failed");

								// get output/storage handler
								OutputHandler * output = parameters.getOutputHandler( *energy );
								INTARNA_CHECK_NOT_NULL(output,"output handler initialization failed");

								// check if we have to add separator for IntaRNA v1 output
								if (reportedInteractions > 0 && dynamic_cast<OutputHandlerIntaRNA1*>(output) != NULL) {
									dynamic_cast<OutputHandlerIntaRNA1*>(output)->addSeparator( true );
								}

								// setup collecting output handler to ensure
								// k-best output per query-target combination
								// and not per region combination if not requested
								OutputHandlerInteractionList bestInteractions(
										(parameters.reportBestPerRegion() ? std::numeric_limits<size_t>::max() : 1 )
											* parameters.getOutputConstraint().reportMax );

								// get interaction prediction handler
								Predictor * predictor = parameters.getPredictor( *energy, bestInteractions );
								INTARNA_CHECK_NOT_NULL(predictor,"predictor initialization failed");

								// run prediction for all range combinations
								// NOTE: RANGE COMBINATIONS CAN NOT BE PARALLELIZED UNLESS
								// (a) predictor.predict() is threadsafe (currently not, eg. due to offset setup)
								// (b) or for each range pair a separate predictor is created (overhead to be checked)
								BOOST_FOREACH(const IndexRange & tRange, parameters.getTargetRanges(*energy, targetNumber)) {
								BOOST_FOREACH(const IndexRange & qRange, parameters.getQueryRanges(*energy, queryNumber)) {

#if INTARNA_MULITHREADING
									#pragma omp critical(intarna_omp_logOutput)
#endif
									{ VLOG(1) <<"predicting interactions for"
											<<" target "<<targetAcc->getSequence().getId()
											<<" (range " <<(tRange+1)<<")"
											<<" and"
											<<" query "<<queryAcc.at(queryNumber)->getSequence().getId()
											<<" (range " <<(qRange+1)<<")"
											<<"..."; }

									predictor->predict(	  tRange
														, queryAcc.at(queryNumber)->getReversedIndexRange(qRange)
														, parameters.getOutputConstraint()
														);

								} // target ranges
								} // query ranges

								// update final output handler
								BOOST_FOREACH( const Interaction * inter, bestInteractions) {
									// forward all reported interactions for all regions to final output handler
									output->add(*inter);
								}

#if INTARNA_MULITHREADING
								#pragma omp atomic update
#endif
								reportedInteractions += output->reported();

								// garbage collection
								 INTARNA_CLEANUP(predictor);
								 INTARNA_CLEANUP(output);
								 INTARNA_CLEANUP(energy);

#if INTARNA_MULITHREADING
							////////////////////// exception handling ///////////////////////////
							} catch (std::exception & e) {
								// ensure exception handling for first failed thread only
								#pragma omp critical(intarna_omp_exception)
								{
									if (!threadAborted) {
										// store exception information
										exceptionPtrDuringOmp = std::make_exception_ptr(e);
										exceptionInfoDuringOmp <<" #thread "<<omp_get_thread_num() <<" #target "<<targetNumber <<" #query " <<queryNumber <<" : "<<e.what();
										// trigger abortion of all threads
										threadAborted = true;
										#pragma omp flush (threadAborted)
									}
								} // omp critical(intarna_omp_exception)
							} catch (...) {
								// ensure exception handling for first failed thread only
								#pragma omp critical(intarna_omp_exception)
								{
									if (!threadAborted) {
										// store exception information
										exceptionPtrDuringOmp = std::current_exception();
										exceptionInfoDuringOmp <<" #thread "<<omp_get_thread_num() <<" #target "<<targetNumber <<" #query " <<queryNumber;
										// trigger abortion of all threads
										threadAborted = true;
										#pragma omp flush (threadAborted)
									}
								} // omp critical(intarna_omp_exception)
							}
						} // if not threadAborted
#endif
					} // for queries

					// write accessibility to file if needed
					parameters.writeTargetAccessibility( *targetAcc );

					// garbage collection
					 INTARNA_CLEANUP(targetAcc);

#if INTARNA_MULITHREADING
				////////////////////// exception handling ///////////////////////////
				} catch (std::exception & e) {
					// ensure exception handling for first failed thread only
					#pragma omp critical(intarna_omp_exception)
					{
						if (!threadAborted) {
							// store exception information
							exceptionPtrDuringOmp = std::make_exception_ptr(e);
							exceptionInfoDuringOmp <<" #thread "<<omp_get_thread_num() <<" #target "<<targetNumber <<" : "<<e.what();
							// trigger abortion of all threads
							threadAborted = true;
							#pragma omp flush (threadAborted)
						}
					} // omp critical(intarna_omp_exception)
				} catch (...) {
					// ensure exception handling for first failed thread only
					#pragma omp critical(intarna_omp_exception)
					{
						if (!threadAborted) {
							// store exception information
							exceptionPtrDuringOmp = std::current_exception();
							exceptionInfoDuringOmp <<" #thread "<<omp_get_thread_num() <<" #target "<<targetNumber;
							// trigger abortion of all threads
							threadAborted = true;
							#pragma omp flush (threadAborted)
						}
					} // omp critical(intarna_omp_exception)
				}
			} // if not threadAborted
#endif
		} // for targets

		// garbage collection
		for (size_t queryNumber=0; queryNumber < queryAcc.size(); queryNumber++) {
			// this is a hack to cleanup the original accessibility object
			Accessibility* queryAccOrig = &(const_cast<Accessibility&>(queryAcc[queryNumber]->getAccessibilityOrigin()) );
			// write accessibility to file if needed
			parameters.writeQueryAccessibility( *queryAccOrig );
			 INTARNA_CLEANUP( queryAccOrig );
			// cleanup (now broken) reverse accessibility object
			 INTARNA_CLEANUP(queryAcc[queryNumber]);
		}

#if INTARNA_MULITHREADING
		if (threadAborted) {
			if (!exceptionInfoDuringOmp.str().empty()) {
				LOG(WARNING) <<"Exception raised for : "<<exceptionInfoDuringOmp.str();
			}
			if (exceptionPtrDuringOmp != NULL) {
				std::rethrow_exception(exceptionPtrDuringOmp);
			}
		}
#endif

	////////////////////// exception handling ///////////////////////////
	} catch (std::exception & e) {
		LOG(WARNING) <<"Exception raised : " <<e.what() <<"\n\n"
			<<"  ==> Please report to the IntaRNA development team! Thanks!\n";
		return -1;
	} catch (...) {
		std::exception_ptr eptr = std::current_exception();
		LOG(WARNING) <<"Unknown exception raised \n\n"
			<<"  ==> Please report to the IntaRNA development team! Thanks!\n";
		return -1;
	}

	  // all went fine
	return 0;
}

