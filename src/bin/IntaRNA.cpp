
#include "IntaRNA/general.h"

// initialize logging for binary
INITIALIZE_EASYLOGGINGPP

#include <iostream>
#include <exception>
#include <memory>
#include <sstream>

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
#include "IntaRNA/OutputHandlerInteractionList.h"

using namespace IntaRNA;

namespace {

struct QueryAccessibilityStorage {
	// Declaration order is intentional: reversed borrows from original and is
	// therefore destroyed first.
	std::unique_ptr<Accessibility> original;
	std::unique_ptr<ReverseAccessibility> reversed;
};

}

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
		// default log file setup
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


		// check if log file set and update all loggers before going on
		if (el::Helpers::commandLineArgs() != NULL && el::Helpers::commandLineArgs()->hasParamWithValue(el::base::consts::kDefaultLogFileParam))
		{
			// default all to file
			el::Loggers::reconfigureAllLoggers(el::ConfigurationType::ToStandardOutput, std::string("false"));
			el::Loggers::reconfigureAllLoggers(el::ConfigurationType::ToFile, std::string("true"));
			// enforec error out to standard output
			el::Loggers::reconfigureAllLoggers(el::Level::Error, el::ConfigurationType::ToStandardOutput, std::string("true"));
			el::Loggers::reconfigureAllLoggers(el::Level::Error, el::ConfigurationType::ToFile, std::string("false"));
		}

		// parse command line parameters
		CommandLineParsing parameters( CommandLineParsing::getPersonality(argc,argv) );
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

		// storage to avoid accessibility recomputation
		std::vector<QueryAccessibilityStorage> queryAcc(parameters.getQuerySequences().size());

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
					queryAcc[qi].original.reset(parameters.getQueryAccessibility(qi));
					INTARNA_CHECK_NOT_NULL(queryAcc[qi].original.get(),"query initialization failed");
					// reverse indexing of target sequence for the computation
					queryAcc[qi].reversed = std::make_unique<ReverseAccessibility>(*queryAcc[qi].original);

					// check if we have to warn about ambiguity
					if (queryAcc[qi].original->getSequence().isAmbiguous()) {
#if INTARNA_MULITHREADING
						#pragma omp critical(intarna_omp_logOutput)
#endif
						VLOG(1) <<"Sequence '"<<queryAcc[qi].original->getSequence().getId()
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
		const bool parallelizeWindowsLoop = !parallelizeTargetLoop && !parallelizeQueryLoop;


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
					std::unique_ptr<Accessibility> targetAcc(parameters.getTargetAccessibility(targetNumber));
					INTARNA_CHECK_NOT_NULL(targetAcc.get(),"target initialization failed");

					// check if we have to warn about ambiguity
					if (targetAcc->getSequence().isAmbiguous()) {
#if INTARNA_MULITHREADING
						#pragma omp critical(intarna_omp_logOutput)
#endif
						{ VLOG(1) <<"Sequence '"<<targetAcc->getSequence().getId()
								<<"' contains ambiguous IUPAC nucleotide encodings. These positions are ignored for interaction computation and are replaced by 'N'.";}
					}

					// second: iterate over all query sequences
#if INTARNA_MULITHREADING
					// this parallelization should only be enabled if the outer target-loop is not parallelized
					# pragma omp parallel for schedule(dynamic) num_threads( parameters.getThreads() ) shared(queryAcc,reportedInteractions,exceptionPtrDuringOmp,exceptionInfoDuringOmp,targetAcc,targetNumber) if(parallelizeQueryLoop)
#endif
					for ( size_t queryIdx = 0; queryIdx < parameters.getQueryNumberForTarget(targetNumber); ++queryIdx )
					{
						// get index of this query wrt. getQuerySequence() and queryAcc()
						const size_t queryNumber = parameters.getQueryIndexForTarget(queryIdx, targetNumber);
#if INTARNA_MULITHREADING
						#pragma omp flush (threadAborted)
						// explicit try-catch-block due to missing OMP exception forwarding
						if (!threadAborted) {
							try {
#endif
								// sanity check
								assert( queryAcc.at(queryNumber).reversed != nullptr );

								// get energy computation handler for both sequences
								std::unique_ptr<InteractionEnergy> energy(parameters.getEnergyHandler( *targetAcc, *queryAcc.at(queryNumber).reversed ));
								INTARNA_CHECK_NOT_NULL(energy.get(),"energy initialization failed");

								// get output/storage handler
								std::unique_ptr<OutputHandler> output(parameters.getOutputHandler( *energy ));
								INTARNA_CHECK_NOT_NULL(output.get(),"output handler initialization failed");

								// setup collecting output handler to ensure
								// k-best output per query-target combination
								// and not per region combination if not requested
								OutputHandlerInteractionList bestInteractions( parameters.getOutputConstraint(*energy),
										(parameters.reportBestPerRegion() ? std::numeric_limits<size_t>::max() : 1 )
											* parameters.getOutputConstraint(*energy).reportMax );

								// run prediction for all range combinations
								for(const IndexRange & tRange : parameters.getTargetRanges(*energy, targetNumber, *targetAcc)) {
								for(const IndexRange & qRange : parameters.getQueryRanges(*energy, queryNumber, *queryAcc.at(queryNumber).original)) {

									// get windows for both ranges
									std::vector<IndexRange> queryWindows = qRange.overlappingWindows(parameters.getWindowWidth(), parameters.getWindowOverlap());
									std::vector<IndexRange> targetWindows = tRange.overlappingWindows(parameters.getWindowWidth(), parameters.getWindowOverlap());

									// iterate over all window combinations
#if INTARNA_MULITHREADING
									// this parallelization should only be enabled if neither the outer target-loop nor the inner query-loop are parallelized
									# pragma omp parallel for schedule(dynamic) collapse(2) num_threads( parameters.getThreads() ) shared(queryAcc,reportedInteractions,exceptionPtrDuringOmp,exceptionInfoDuringOmp,targetAcc,targetNumber,queryWindows,targetWindows, bestInteractions, energy) if(parallelizeWindowsLoop)
#endif									
									for (int qNumWindow = 0; qNumWindow < queryWindows.size(); ++qNumWindow) {
									for (int tNumWindow = 0; tNumWindow < targetWindows.size(); ++tNumWindow) {									
#if INTARNA_MULITHREADING
										#pragma omp flush (threadAborted)
										// explicit try-catch-block due to missing OMP exception forwarding
										if (!threadAborted) {
											try {
#endif										
										
												IndexRange qWindow = queryWindows.at(qNumWindow);
												IndexRange tWindow = targetWindows.at(tNumWindow);
#if INTARNA_MULITHREADING
												#pragma omp critical(intarna_omp_logOutput)
#endif
												{ VLOG(1) <<"predicting interactions for"
														<<" target "<<targetAcc->getSequence().getId()
														<<" (range " <<(tWindow+1)<<")"
														<<" and"
														<<" query "<<queryAcc.at(queryNumber).reversed->getSequence().getId()
														<<" (range " <<(qWindow+1)<<")"
#if INTARNA_MULITHREADING
#if INTARNA_IN_DEBUG_MODE

														<<" in thread "<<omp_get_thread_num()
#endif
#endif
														<<" ..."; }
		
												// get interaction prediction handler
												std::unique_ptr<Predictor> predictor(parameters.getPredictor( *energy, bestInteractions ));
												INTARNA_CHECK_NOT_NULL(predictor.get(),"predictor initialization failed");
		
												// run prediction for this window combination
												predictor->predict(	  tWindow
																	, queryAcc.at(queryNumber).reversed->getReversedIndexRange(qWindow)
																	);
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
									}} // window combinations
								} // target ranges
								} // query ranges
#if INTARNA_MULITHREADING
								#pragma omp critical(intarna_omp_outputHandlerUpdate)
#endif
								{// update final output handler
									// copy partition function information if available
									output->incrementZ( bestInteractions.getZ() );
									// forward all reported interactions for all regions to final output handler
									for( const Interaction * inter : bestInteractions) {
										output->add(*inter);
									}
								}

#if INTARNA_MULITHREADING
								#pragma omp atomic update
#endif
								reportedInteractions += output->reported();

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

		// write successfully initialized query accessibilities if needed
		for (const QueryAccessibilityStorage & queryAccessibility : queryAcc) {
			if (queryAccessibility.original != nullptr && queryAccessibility.reversed != nullptr) {
				parameters.writeQueryAccessibility(*queryAccessibility.original);
			}
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
			<<"  ==> Please report (including input) to the IntaRNA development team! Thanks!\n";
		el::Loggers::flushAll();
		return -1;
	} catch (...) {
		std::exception_ptr eptr = std::current_exception();
		LOG(WARNING) <<"Unknown exception raised \n\n"
			<<"  ==> Please report (including input) to the IntaRNA development team! Thanks!\n";
		el::Loggers::flushAll();
		return -1;
	}

	  // all went fine
	el::Loggers::flushAll();
	return 0;
}

