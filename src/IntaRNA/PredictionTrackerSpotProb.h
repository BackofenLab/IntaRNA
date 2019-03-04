
#ifndef INTARNA_PREDICTIONTRACKERSPOTPROB_H_
#define INTARNA_PREDICTIONTRACKERSPOTPROB_H_

#include "IntaRNA/PredictionTracker.h"
#include "IntaRNA/InteractionEnergy.h"

#include <iostream>

#include <boost/algorithm/string.hpp>

namespace IntaRNA {

/**
 * Computes probability estimations that a given set of interaction spots are
 * covered by an interaction. Furthermore it computes the probability that none
 * of the targeted spots is covered by an interaction.
 *
 * The information is written to stream on destruction.
 */
class PredictionTrackerSpotProb: public PredictionTracker
{

public:

	//! regex string encoding of a single spot string encoding
	static const std::string str_spot;
	//! regex to match a valid list of spot encodings
	static const boost::regex regexSpotString;

public:

	/**
	 * Constructs a PredictionTracker that collects probability information
	 * for a set of interaction spots by computing the Boltzmann probability
	 * of any interaction enclosing the positions.
	 *
	 * @param energy the energy function used for energy calculation
	 * @param spots the interaction spots to track for their probabilities
	 * @param outStreamName the stream name where the probability data
	 *        is to be written to. use STDOUT/STDERR for the respective stream.
	 *        Otherwise, an according file is created
	 */
	PredictionTrackerSpotProb(
				const InteractionEnergy & energy
				, const std::string & spots
				, const std::string & outStreamName
			);

	/**
	 * Constructs a PredictionTracker that collects probability information
	 * for a set of interaction spots by computing the Boltzmann probability
	 * of any interaction enclosing the positions.
	 *
	 * @param energy the energy function used for energy calculation
	 * @param spots the interaction spots to track for their probabilities
	 * @param outStream the stream where the probability data
	 *        is to be written to.
	 */
	PredictionTrackerSpotProb(
				const InteractionEnergy & energy
				, const std::string & spots
				, std::ostream & outStream
			);


	/**
	 * destruction: write the probabilities to stream.
	 */
	virtual ~PredictionTrackerSpotProb();


	/**
	 * Updates the probability information for each Predictor.updateOptima() call.
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1
	 * @param energy the overall energy of the interaction site
	 */
	virtual
	void
	updateOptimumCalled( const size_t i1, const size_t j1
						, const size_t i2, const size_t j2
						, const E_type energy
						);

protected:

	class Spot {
	public:
		//! index in seq1
		size_t idx1;
		//! index in seq2
		size_t idx2;
		//! partition function of all interactions covering the spot
		Z_type Z;

		/* default construction
		 */
		Spot()
		 : idx1(std::string::npos), idx2(std::string::npos), Z(0.0)
		{}

		/* construction from data
		 * @param idx1 index in seq1
		 * @param idx2 index in seq2
		 */
		Spot(size_t idx1, size_t idx2)
		 : idx1(idx1), idx2(idx2), Z(0.0)
		{}

		/* construction from string
		 * @param spotString string encoding of a spot
		 */
		Spot(const std::string & spotString )
		 : idx1(std::stol(spotString.substr(0,spotString.find('&'))))
			, idx2(std::stol(spotString.substr(spotString.find('&')+1)))
			, Z(0.0)
		{
#if INTARNA_IN_DEBUG_MODE
			// check regex
			if (!boost::regex_match( spotString, boost::regex(PredictionTrackerSpotProb::str_spot), boost::match_perl )) {
				throw std::runtime_error("PredictionTrackerSpotProb::Spot("+spotString+") does not match its encoding regular expression");
			}
			assert(idx1 > 0);
			assert(idx2 > 0);
#endif
			// correct user input by -1 to start indexing with 0
			idx1--;
			idx2--;
		}
	};


protected:

	//! energy handler used for predictions
	const InteractionEnergy & energy;

	//! the stream to write the probabilities to
	std::ostream * outStream;

	//! whether or not outStream is to be deleted on destruction
	const bool deleteOutStream;

	//! tracking information for each spot
	std::vector<Spot> spots;

	//! partition function of all interaction not covering a tracked spot
	Z_type noSpotZ;

	//! overall partition function of all reported interactions
	Z_type overallZ;

protected:



};

//////////////////////////////////////////////////////////////////////

} // namespace

#endif /* PREDICTIONTRACKERSPOTPROB_H_ */
