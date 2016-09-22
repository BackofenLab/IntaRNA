/*
 * intaRNA.cpp
 *
 * Main file for IntaRNA 2.0
 *
 *  Created on: 17.06.2014
 *      Author: Martin Mann
 */


#include "config.h"
#include "general.h"

#include <iostream>

#include "CommandLineParsing.h"

#include "RnaSequence.h"

#include "Accessibility.h"
#include "AccessibilityDisabled.h"
#include "AccessibilityVienna.h"

#include "Energy.h"
#include "EnergyBasePair.h"
#include "EnergyVienna.h"

#include "Predictor.h"

#include "OutputHandler.h"
#include "OutputHandlerText.h"

/////////////////////////////////////////////////////////////////////
/**
 * program main entry
 *
 * @param argc number of program arguments
 * @param argv array of program arguments of length argc
 */
int main(int argc, char **argv) {

	try {

		// parse command line parameters
		CommandLineParsing parameters(std::cout);
		{
			int retCode = parameters.parse( argc, argv );
			if (retCode <= 0) {
				return retCode;
			}
		}


		// run prediction for all pairs of sequences
		// first: iterate over all query sequences
		for ( size_t queryNumber = 0; queryNumber < parameters.getQuerySequences().size(); ++queryNumber )
		{

			// get accessibility handler
			Accessibility * queryAcc = parameters.getQueryAccessibility(queryNumber);
			CHECKNOTNULL(queryAcc,"query initialization failed");

			// check if we have to warn about ambiguity
			if (queryAcc->getSequence().isAmbiguous()) {
				std::cout <<"# WARNING: Sequence '"<<queryAcc->getSequence().getId() <<"' contains ambiguous nucleotide encodings. These positions are ignored for interaction computation."
						<<std::endl;
			}

			// second: iterate over all target sequences to get all pairs to predict for
			for ( size_t targetNumber = 0; targetNumber < parameters.getTargetSequences().size(); ++targetNumber )
			{

				// get target accessibility handler
				Accessibility * targetAccOrig = parameters.getTargetAccessibility(targetNumber);
				CHECKNOTNULL(targetAccOrig,"target initialization failed");
				// reverse indexing of target sequence for the computation
				ReverseAccessibility * targetAcc = new ReverseAccessibility(*targetAccOrig);

				// check if we have to warn about ambiguity
				if (targetAcc->getSequence().isAmbiguous()) {
					std::cout <<"# WARNING: Sequence '"<<targetAcc->getSequence().getId() <<"' contains ambiguous nucleotide encodings. These positions are ignored for interaction computation."
							<<std::endl;
				}

				// get energy computation handler for both sequences
				Energy* energy = parameters.getEnergyHandler( *queryAcc, *targetAcc );
				CHECKNOTNULL(energy,"energy initialization failed");

				// TODO get seed computation handler
//				Seed* seed = parameters.getSeedHandler( *energy );

				// TODO get output/storage handler
				OutputHandler * output = NULL;
//				CHECKNOTNULL(output,"output handler initialization failed");

				// TEST ####################

					output = new OutputHandlerText(std::cout);
					Interaction i(queryAcc->getSequence(), targetAcc->getAccessibilityOrigin().getSequence());
					i.addInteraction(13,17);
					i.addInteraction(16,14);
					i.addInteraction(15,16);
					i.addInteraction(27,2);
					i.sort();
					output->add(i);


				// TEST END ####################


				// TODO get interaction prediction handler
				Predictor * predictor = NULL;
//				CHECKNOTNULL(predictor,"predictor initialization failed");

				// TODO run prediction

				// garbage collection
				CLEANUP(predictor)
				CLEANUP(output)
				CLEANUP(energy)
				CLEANUP(targetAcc)
				CLEANUP(targetAccOrig)
				CLEANUP(queryAcc)
			}

		}



	////////////////////// exception handling ///////////////////////////
	} catch (std::exception & e) {
		std::cerr <<"\n Exception raised : " <<e.what() <<"\n";
		return -1;
	}

	  // all went fine
	return 0;
}

