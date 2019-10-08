
#ifndef INTARNA_PREDICTIONTRACKERPROFILESPOTPROB_H_
#define INTARNA_PREDICTIONTRACKERPROFILESPOTPROB_H_

#include "IntaRNA/PredictionTracker.h"
#include "IntaRNA/InteractionEnergy.h"

#include <iostream>

#include <boost/algorithm/string.hpp>

namespace IntaRNA {

/**
 * Collects for each sequence position the partition function of any interaction
 * enclosing this position.
 *
 * The profile(s) of this information (normalized by the overall partitions
 * function; ie the spot interaction probability) is written to stream on
 * destruction.
 */
class PredictionTrackerProfileSpotProb: public PredictionTracker
{

public:

	/**
	 * Constructs a PredictionTracker that collects for each positions of a
	 * sequence the partition function and thus spot probability of any
	 * interaction enclosing this position.
	 *
	 * @param energy the energy function used for energy calculation
	 * @param seq1streamName the stream name where the profile data for seq1
	 *        is to be written to. use STDOUT/STDERR for the respective stream.
	 *        Otherwise, an according file is created;
	 *        if empty, no data for seq1 is collected
	 * @param seq2streamName the stream name where the profile data for seq2
	 *        is to be written to. use STDOUT/STDERR for the respective stream.
	 *        Otherwise, an according file is created;
	 *        if empty, no data for seq2 is collected
	 * @param NA_string the output string representation if a value is not available
	 *        for profile output
	 * @param sep the column separator to be used in profile output
	 */
	PredictionTrackerProfileSpotProb(
				const InteractionEnergy & energy
				, const std::string & seq1streamName
				, const std::string & seq2streamName
				, const std::string & NA_string = "NA"
				, const std::string & sep = ";"
			);

	/**
	 * Constructs an PredictionTracker that collected for each positions of a
	 * sequence the partition function and thus spot probability of any
	 * interaction covering this position.
	 *
	 * @param energy the energy function used
	 * @param seq1stream if non-NULL, the stream where the profile data for seq1
	 *        is to be written to; if NULL, no data for seq1 is collected
	 * @param seq2stream if non-NULL, the stream where the profile data for seq2
	 *        is to be written to; if NULL, no data for seq2 is collected
	 * @param NA_string the output string representation if a value is not available
	 *        for profile output
	 * @param sep the column separator to be used in profile output
	 */
	PredictionTrackerProfileSpotProb(
				const InteractionEnergy & energy
				, std::ostream * seq1stream
				, std::ostream * seq2stream
				, const std::string & NA_string = "NA"
				, const std::string & sep = ";"
			);

	/**
	 * destruction: write the profile(s) to the according streams.
	 */
	virtual ~PredictionTrackerProfileSpotProb();


	/**
	 * Updates the profile information for each Predictor.updateOptima() call.
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

	//! if non-NULL, the stream to write the spotProb-profile for seq1 to
	std::ostream * seq1stream;

	//! if non-NULL, the stream to write the spotProb-profile for seq2 to
	std::ostream * seq2stream;

	//! the output string representation if a value is not available for profile output
	const std::string NA_string;

	//! the output string representation of a column separator
	const std::string sep;

	//! container definition for partition function profile data
	typedef std::vector<Z_type> ZProfile;

	//! the position-wise partition function values for seq1
	ZProfile seq1Z;

	//! the position-wise partition function values for seq2
	ZProfile seq2Z;

	//! the overall partition function to be used for normalization
	Z_type overallZ;

	/**
	 * Updates a given partition function profile if not empty.
	 *
	 * @param profile the partition function profile to update
	 * @param i the first index to update (inclusive)
	 * @param j the last index to update (inclusive)
	 * @param boltzmannWeight the Boltzmann weight to be used for the update
	 */
	static
	void
	updateProfile(	  ZProfile & profile
					, const size_t i
					, const size_t j
					, const Z_type boltzmannWeight);

	/**
	 * Writes profile data to stream.
	 *
	 * @param out the output stream to write to
	 * @param begin the begin iterator of the partition function data to write
	 * @param end the end iterator (exclusive) of the Z data to write
	 * @param overallZ the overall partition function to be used for
	 *           normalization
	 * @param rna the RNA the data is about
	 * @param NA_string the string to be used for missing entries
	 * @param sep the column separator to be used in profile output
	 */
	template < typename ZProfileIterator >
	static
	void
	writeProfile( std::ostream &out
				, const ZProfileIterator & begin
				, const ZProfileIterator & end
				, const Z_type overallZ
				, const RnaSequence & rna
				, const std::string & NA_string
				, const std::string & sep
				);


};

//////////////////////////////////////////////////////////////////////

template < typename ZProfileIterator >
inline
void
PredictionTrackerProfileSpotProb::
writeProfile( std::ostream &out
			, const ZProfileIterator & begin
			, const ZProfileIterator & end
			, const Z_type overallZ
			, const RnaSequence & rna
			, const std::string & NA_string
			, const std::string & sep
			)
{
	// write in CSV-like format (column data)

	const bool noZ = (overallZ==0.0);

	// print header : seq.ID ; spotProb
	out <<"idx"<<sep<<boost::replace_all_copy(rna.getId(), sep, "_")<<sep<<"spotProb" <<'\n';
	// print spot probability data
	size_t i=0;
	for (ZProfileIterator curZ = begin; curZ!=end; curZ++) {
		out
			// out index
			<<rna.getInOutIndex(i)<<sep
			// out nucleotide (index starts with 0)
			<<rna.asString().at(i)<<sep
			;
		// out infinity replacement if needed
		if ( noZ || Z_isINF( *curZ ) ) {
			out<<NA_string <<'\n';
		} else {
			out <<((*curZ)/overallZ) <<'\n';
		}
		// increase position counter
		i++;
	}
}

//////////////////////////////////////////////////////////////////////

} // namespace

#endif /* PREDICTIONTRACKERPROFILESPOTPROB_H_ */
