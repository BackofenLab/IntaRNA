
#ifndef INTARNA_PREDICTIONTRACKERSPOTPROBALL_H_
#define INTARNA_PREDICTIONTRACKERSPOTPROBALL_H_

#include "IntaRNA/PredictionTracker.h"
#include "IntaRNA/InteractionEnergy.h"

#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/numeric/ublas/matrix.hpp>


namespace IntaRNA {

/**
 * Collects for each intermolecular index pair the probability that this pair
 * is covered by an interaction.
 *
 * The pair data is written to stream on destruction.
 */
class PredictionTrackerSpotProbAll: public PredictionTracker
{

public:

	/**
	 * Constructs a PredictionTracker that collected for each index pair of two
	 * sequences the probability to be covered by an interaction.
	 *
	 * Note, if a filename is provided, its stream is closed on destruction of
	 * this object!
	 *
	 * @param energy the energy function used for energy calculation
	 * @param streamName the stream name where the data
	 *        is to be written to. use STDOUT/STDERR for the respective stream.
	 *        Otherwise, an according file is created (has to be non-empty)
	 * @param NA_string the output string representation if a value is not available
	 *        for profile output
	 * @param sep the column separator to be used in profile output
	 */
	PredictionTrackerSpotProbAll(
				const InteractionEnergy & energy
				, const std::string & streamName
				, const std::string & NA_string = "NA"
				, const std::string & sep = ";"
			);

	/**
	 * Constructs a PredictionTracker that collected for each index pair of two
	 * sequences the probability to be covered by an interaction.
	 *
	 * Note, the stream is NOT closed nor deleted on destruction of this object!
	 *
	 * @param energy the energy function used
	 * @param outStream the stream where the data
	 *        is to be written to (has to be non-null)
	 * @param NA_string the output string representation if a value is not available
	 *        for profile output
	 * @param sep the column separator to be used in profile output
	 */
	PredictionTrackerSpotProbAll(
				const InteractionEnergy & energy
				, std::ostream * outStream
				, const std::string & NA_string = "NA"
				, const std::string & sep = ";"
			);

	/**
	 * destruction: write the profile(s) to the according streams.
	 */
	virtual ~PredictionTrackerSpotProbAll();


	/**
	 * Updates the partition function information for each
	 * Predictor.updateOptima() call.
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

	//! energy handler used for predictions
	const InteractionEnergy & energy;

	//! whether or not the streams are to be deleted on destruction
	const bool deleteStreamsOnDestruction;

	//! the stream to write the minE-data to
	std::ostream * outStream;

	//! the output string representation if a value is not available for output
	const std::string NA_string;

	//! the output string representation of column separators
	const std::string sep;

	//! overall partition function
	Z_type overallZ;

	//! matrix type to hold the partition function for each index pair
	typedef boost::numeric::ublas::matrix<Z_type> Z2dMatrix;

	//! the index-pair-wise minimal energy values
	Z2dMatrix pairZ;



	/**
	 * Writes profile data to stream.
	 *
	 * @param out the output stream to write to
	 * @param pairZ the partition function data to write
	 * @param overallZ the overall partition function for pairZ
	 * @param energy the energy function used
	 * @param NA_string the string to be used for missing entries
	 * @param sep the column separator to be used in profile output
	 */
	static
	void
	writeData( std::ostream &out
				, const Z2dMatrix & pairZ
				, const Z_type & overallZ
				, const InteractionEnergy & energy
				, const std::string & NA_string
				, const std::string & sep
				);


};

//////////////////////////////////////////////////////////////////////

} // namespace

#endif /* INTARNA_PREDICTIONTRACKERSPOTPROBALL_H_ */
