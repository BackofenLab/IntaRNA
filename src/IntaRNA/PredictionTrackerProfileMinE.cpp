
#include "IntaRNA/PredictionTrackerProfileMinE.h"

namespace IntaRNA {

//////////////////////////////////////////////////////////////////////

PredictionTrackerProfileMinE::
PredictionTrackerProfileMinE(
		const InteractionEnergy & energy
		, const std::string & seq1streamName
		, const std::string & seq2streamName
		, const std::string & E_INF_string
		, const std::string & sep_string
	)
 :	PredictionTracker()
	, energy(energy)
	, deleteStreamsOnDestruction(true)
	, seq1stream(NULL)
	, seq2stream(NULL)
	, E_INF_string(E_INF_string)
	, sep_string(sep_string)
	, seq1minE( seq1streamName.empty() ? 0 : energy.size1(), E_INF ) // init E_INF
	, seq2minE( seq2streamName.empty() ? 0 : energy.size2(), E_INF ) // init E_INF
{
#if INTARNA_IN_DEBUG_MODE
	// check separator
	if (sep_string.empty()) {
		throw std::runtime_error("PredictionTrackerProfileMinE() empty separator provided");
	}
#endif
	// open streams
	if (!seq1streamName.empty()) {
		seq1stream = newOutputStream( seq1streamName );
		if (seq1stream == NULL) {
			throw std::runtime_error("PredictionTrackerProfileMinE() : could not open output stream '"+seq1streamName+"' for writing");
		}
	}
	if (!seq2streamName.empty()) {
		seq2stream = newOutputStream( seq2streamName );
		if (seq2stream == NULL) {
			// cleanup seq1 stream before throwing exception
			if (seq1stream != NULL) {
				deleteOutputStream(seq1stream);
			}
			throw std::runtime_error("PredictionTrackerProfileMinE() : could not open output stream '"+seq2streamName+"' for writing");
		}
	}
}

//////////////////////////////////////////////////////////////////////

PredictionTrackerProfileMinE::
PredictionTrackerProfileMinE(
		const InteractionEnergy & energy
		, std::ostream * seq1stream
		, std::ostream * seq2stream
		, const std::string & E_INF_string
		, const std::string & sep_string
	)
 :	PredictionTracker()
	, energy(energy)
	, deleteStreamsOnDestruction(false)
	, seq1stream(seq1stream)
	, seq2stream(seq2stream)
	, E_INF_string(E_INF_string)
	, sep_string(sep_string)
	, seq1minE( seq1stream==NULL ? 0 : energy.size1(), E_INF ) // init E_INF
	, seq2minE( seq2stream==NULL ? 0 : energy.size2(), E_INF ) // init E_INF
{
#if INTARNA_IN_DEBUG_MODE
	// check separator
	if (sep_string.empty()) {
		throw std::runtime_error("PredictionTrackerProfileMinE() empty separator provided");
	}
#endif
}

//////////////////////////////////////////////////////////////////////

PredictionTrackerProfileMinE::
~PredictionTrackerProfileMinE()
{
	// write profiles to streams
	if (seq1stream != NULL) {
		writeProfile( *seq1stream
					, seq1minE.begin()
					, seq1minE.end()
					, energy.getAccessibility1().getSequence()
					, E_INF_string
					, sep_string
					);
	}
	if (seq2stream != NULL) {
		writeProfile( *seq2stream
					, seq2minE.begin()
					, seq2minE.end()
					, energy.getAccessibility2().getAccessibilityOrigin().getSequence()
					, E_INF_string
					, sep_string
					);
	}

	// clean up if file pointers were created in constructor
	if (deleteStreamsOnDestruction) {
		deleteOutputStream( seq1stream );
		deleteOutputStream( seq2stream );
	}
}

//////////////////////////////////////////////////////////////////////

void
PredictionTrackerProfileMinE::
updateOptimumCalled( const size_t i1, const size_t j1
					, const size_t i2, const size_t j2
					, const E_type curE
					)
{
	// get mapped index positions in overall sequence
	Interaction::BasePair bp_l = energy.getBasePair(
							i1 == RnaSequence::lastPos ? energy.size1()-1 : i1
							,i2 == RnaSequence::lastPos ? energy.size2()-1 : i2);
	Interaction::BasePair bp_r = energy.getBasePair(
							j1 == RnaSequence::lastPos ? energy.size1()-1 : j1
							,j2 == RnaSequence::lastPos ? energy.size2()-1 : j2);

	// update seq1 profile
	updateProfile( seq1minE
					, bp_l.first, bp_r.first
					, curE
					);

	// update seq2 profile
	updateProfile( seq2minE
					// swap base pairs due to reversed indexing
					, bp_r.second, bp_l.second
					, curE
					);
}


//////////////////////////////////////////////////////////////////////

void
PredictionTrackerProfileMinE::
updateProfile(	  MinEProfile & profile
				, const size_t i
				, const size_t j
				, const E_type E)
{
	// update if the profile is not empty (is to be filled)
	if (! profile.empty()) {
#if INTARNA_IN_DEBUG_MODE
		if (i>=profile.size() || j>=profile.size()) throw std::runtime_error("PredictionTrackerProfileMinE::updateProfile() : index range ["+toString(i)+","+toString(j)+"] exceeds sequence length "+toString(profile.size()));
		if (i>j) throw std::runtime_error("PredictionTrackerProfileMinE::updateProfile() : i "+toString(i)+" > j "+toString(j));
#endif
		// update profile data
		for (size_t k=i; k<=j; k++) {
			// check if E is smaller than current minE
			if ( E < profile.at(k) ) {
				// write new minimum
				profile[k] = E;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////


} // namespace
