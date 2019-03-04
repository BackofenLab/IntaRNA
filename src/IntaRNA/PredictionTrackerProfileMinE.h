
#ifndef INTARNA_PREDICTIONTRACKERPROFILEMINE_H_
#define INTARNA_PREDICTIONTRACKERPROFILEMINE_H_

#include "IntaRNA/PredictionTracker.h"
#include "IntaRNA/InteractionEnergy.h"

#include <iostream>

#include <boost/algorithm/string.hpp>

namespace IntaRNA {

/**
 * Collects for each sequence position the minimal energy of any interaction
 * enclosing this position.
 *
 * The profile(s) of this information is written to stream on destruction.
 */
class PredictionTrackerProfileMinE: public PredictionTracker
{

public:

	/**
	 * Constructs a PredictionTracker that collects for each positions of a
	 * sequence the minimal energy of any interaction enclosing this position.
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
	 * @param E_INF_string the output string representation of E_INF values in
	 *        the profile output
	 */
	PredictionTrackerProfileMinE(
				const InteractionEnergy & energy
				, const std::string & seq1streamName
				, const std::string & seq2streamName
				, const std::string E_INF_string = "NA"
			);

	/**
	 * Constructs an PredictionTracker that collected for each positions of a
	 * sequence the minimal energy of any interaction covering this position.
	 *
	 * @param energy the energy function used
	 * @param seq1stream if non-NULL, the stream where the profile data for seq1
	 *        is to be written to; if NULL, no data for seq1 is collected
	 * @param seq2stream if non-NULL, the stream where the profile data for seq2
	 *        is to be written to; if NULL, no data for seq2 is collected
	 * @param E_INF_string the output string representation of E_INF values in
	 *        the profile output
	 */
	PredictionTrackerProfileMinE(
				const InteractionEnergy & energy
				, std::ostream * seq1stream
				, std::ostream * seq2stream
				, const std::string E_INF_string = "NA"
			);

	/**
	 * destruction: write the profile(s) to the according streams.
	 */
	virtual ~PredictionTrackerProfileMinE();


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

	//! if non-NULL, the stream to write the minE-profile for seq1 to
	std::ostream * seq1stream;

	//! if non-NULL, the stream to write the minE-profile for seq2 to
	std::ostream * seq2stream;

	//! the output string representation of E_INF values in the profile output
	const std::string E_INF_string;

	//! container definition for minE profile data
	typedef std::vector<E_type> MinEProfile;

	//! the position-wise minimal energy values for seq1
	MinEProfile seq1minE;

	//! the position-wise minimal energy values for seq2
	MinEProfile seq2minE;

	/**
	 * Updates a given minE profile if not empty.
	 *
	 * @param profile the minE profile to update
	 * @param i the first index to update (inclusive)
	 * @param j the last index to update (inclusive)
	 * @param E the energy to be used for the update
	 */
	static
	void
	updateProfile(	  MinEProfile & profile
					, const size_t i
					, const size_t j
					, const E_type E);

	/**
	 * Writes profile data to stream.
	 *
	 * @param out the output stream to write to
	 * @param begin the begin iterator of the minE data to write
	 * @param end the end iterator (exclusive) of the minE data to write
	 * @param rna the RNA the data is about
	 * @param E_INF_string the string to be used for E_INF entries
	 */
	template < typename MinEProfileIterator >
	static
	void
	writeProfile( std::ostream &out
				, const MinEProfileIterator & begin
				, const MinEProfileIterator & end
				, const RnaSequence & rna
				, const std::string & E_INF_string );


};

//////////////////////////////////////////////////////////////////////

template < typename MinEProfileIterator >
inline
void
PredictionTrackerProfileMinE::
writeProfile( std::ostream &out
			, const MinEProfileIterator & begin
			, const MinEProfileIterator & end
			, const RnaSequence & rna
			, const std::string & E_INF_string )
{
	// write in CSV-like format (column data)

	// print header : seq.ID ; minE
	out <<"idx;"<<boost::replace_all_copy(rna.getId(), ";", "_")<<";minE" <<'\n';
	// print minE data
	size_t i=1;
	for (MinEProfileIterator curE = begin; curE!=end; curE++) {
		out
			// out index
			<<i<<';'
			// out nucleotide (index starts with 0)
			<<rna.asString().at(i-1)<<';'
			;
		// out infinity replacement if needed
		if ( E_isINF( *curE ) ) {
			out<<E_INF_string <<'\n';
		} else {
			out <<E_2_Ekcal(*curE) <<'\n';
		}
		// increase position counter
		i++;
	}
}

//////////////////////////////////////////////////////////////////////

} // namespace

#endif /* PREDICTIONTRACKERPROFILEMINE_H_ */
