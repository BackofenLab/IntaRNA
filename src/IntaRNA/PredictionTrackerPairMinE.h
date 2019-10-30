
#ifndef INTARNA_PREDICTIONTRACKERPAIRMINE_H_
#define INTARNA_PREDICTIONTRACKERPAIRMINE_H_

#include "IntaRNA/PredictionTracker.h"
#include "IntaRNA/InteractionEnergy.h"

#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/numeric/ublas/matrix.hpp>


namespace IntaRNA {

/**
 * Collects for each intermolecular index pair the minimal energy of any interaction
 * covering this pair (even if not forming a base pair).
 *
 * The pair data is written to stream on destruction.
 */
class PredictionTrackerPairMinE: public PredictionTracker
{

public:

	/**
	 * Constructs a PredictionTracker that collected for index pair of two
	 * sequences the minimal energy of any interaction covering this position.
	 *
	 * Note, if a filename is provided, its stream is closed on destruction of
	 * this object!
	 *
	 * @param energy the energy function used for energy calculation
	 * @param streamName the stream name where the data
	 *        is to be written to. use STDOUT/STDERR for the respective stream.
	 *        Otherwise, an according file is created (has to be non-empty)
	 * @param E_INF_string the output string representation of E_INF values in
	 *        the profile output
	 * @param sep the column separator to be used in profile output
	 */
	PredictionTrackerPairMinE(
				const InteractionEnergy & energy
				, const std::string & streamName
				, const std::string & E_INF_string = "NA"
				, const std::string & sep = ";"
			);

	/**
	 * Constructs a PredictionTracker that collected for index pair of two
	 * sequences the minimal energy of any interaction covering this position.
	 *
	 * Note, the stream is NOT closed nor deleted on destruction of this object!
	 *
	 * @param energy the energy function used
	 * @param outStream the stream where the data
	 *        is to be written to (has to be non-null)
	 * @param E_INF_string the output string representation of E_INF values in
	 *        the profile output
	 * @param sep the column separator to be used in profile output
	 */
	PredictionTrackerPairMinE(
				const InteractionEnergy & energy
				, std::ostream * outStream
				, const std::string & E_INF_string = "NA"
				, const std::string & sep = ";"
			);

	/**
	 * destruction: write the pair information to the according stream.
	 */
	virtual ~PredictionTrackerPairMinE();


	/**
	 * Updates the minE information for each Predictor.updateOptima() call.
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

	//! the output string representation of E_INF values in the profile output
	const std::string E_INF_string;

	//! the output string representation of column separators
	const std::string sep_string;

	//! matrix type to hold the mfe energies and boundaries for interaction site starts
	typedef boost::numeric::ublas::matrix<E_type> E2dMatrix;

	//! the index-pair-wise minimal energy values
	E2dMatrix pairMinE;

	/**
	 * Writes minE data to stream.
	 *
	 * @param out the output stream to write to
	 * @param pairMinE the minE data to write
	 * @param energy the energy function used
	 * @param E_INF_string the string to be used for E_INF entries
	 * @param sep_string the output string representation of column separators
	 */
	static
	void
	writeData( std::ostream &out
				, const E2dMatrix & pairMinE
				, const InteractionEnergy & energy
				, const std::string & E_INF_string
				, const std::string & sep_string );


};

//////////////////////////////////////////////////////////////////////

} // namespace

#endif /* PREDICTIONTRACKERPAIRMINE_H_ */
