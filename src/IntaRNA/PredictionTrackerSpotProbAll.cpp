
#include "IntaRNA/PredictionTrackerSpotProbAll.h"

namespace IntaRNA {

//////////////////////////////////////////////////////////////////////

PredictionTrackerSpotProbAll::
PredictionTrackerSpotProbAll(
		const InteractionEnergy & energy
		, const std::string & streamName
		, const std::string & NA_string
		, const std::string & sep
	)
 :	PredictionTracker()
	, energy(energy)
	, deleteStreamsOnDestruction(true)
	, outStream(NULL)
	, NA_string(NA_string)
	, sep(sep)
	, overallZ( Z_type(0.0) )
	, pairZ( energy.size1(), energy.size2(), (Z_type)0.0 ) // init 0
{
#if INTARNA_IN_DEBUG_MODE
	if (streamName.empty()) {
		throw std::runtime_error("PredictionTrackerSpotProbAll() : streamName empty");
	}
	// check separator
	if (sep.empty()) {
		throw std::runtime_error("PredictionTrackerSpotProbAll() empty separator provided");
	}
#endif
	// open streams
	outStream = newOutputStream( streamName );
	if (outStream == NULL) {
		throw std::runtime_error("PredictionTrackerSpotProbAll() : could not open output stream '"+streamName+"' for writing");
	}
}

//////////////////////////////////////////////////////////////////////

PredictionTrackerSpotProbAll::
PredictionTrackerSpotProbAll(
		const InteractionEnergy & energy
		, std::ostream * outStream
		, const std::string & NA_string
		, const std::string & sep
	)
 :	PredictionTracker()
	, energy(energy)
	, deleteStreamsOnDestruction(false)
	, outStream(outStream)
	, NA_string(NA_string)
	, sep(sep)
	, overallZ( Z_type(0.0) )
	, pairZ( energy.size1(), energy.size2(), (Z_type)0.0 ) // init 0
{
#if INTARNA_IN_DEBUG_MODE
	// check separator
	if (sep.empty()) {
		throw std::runtime_error("PredictionTrackerSpotProbAll() empty separator provided");
	}
#endif
}

//////////////////////////////////////////////////////////////////////

PredictionTrackerSpotProbAll::
~PredictionTrackerSpotProbAll()
{
	writeData( *outStream
				, pairZ
				, overallZ
				, energy
				, NA_string
				, sep );

	// clean up if file pointers were created in constructor
	if (deleteStreamsOnDestruction) {
		deleteOutputStream( outStream );
	}
}

//////////////////////////////////////////////////////////////////////

void
PredictionTrackerSpotProbAll::
updateOptimumCalled( const size_t i1_, const size_t j1_
					, const size_t i2_, const size_t j2_
					, const E_type curE
					)
{
	// get valid indices
	const size_t i1 =  (i1_ == RnaSequence::lastPos ? energy.size1()-1 : i1_);
	const size_t j1 =  (j1_ == RnaSequence::lastPos ? energy.size1()-1 : j1_);
	const size_t i2 =  (i2_ == RnaSequence::lastPos ? energy.size2()-1 : i2_);
	const size_t j2 =  (j2_ == RnaSequence::lastPos ? energy.size2()-1 : j2_);

#if INTARNA_IN_DEBUG_MODE
	// sanity checks
	if (i1>=energy.size1() || j1>=energy.size1()) throw std::runtime_error("PredictionTrackerSpotProbAll::updateProfile() : index-1 range ["+toString(i1)+","+toString(j1)+"] exceeds sequence-1 length "+toString(energy.size1()));
	if (i1>j1) throw std::runtime_error("PredictionTrackerSpotProbAll::updateProfile() : i1 "+toString(i1)+" > j1 "+toString(j1));
	if (i2>=energy.size2() || j2>=energy.size2()) throw std::runtime_error("PredictionTrackerSpotProbAll::updateProfile() : index-2 range ["+toString(i2)+","+toString(j2)+"] exceeds sequence-2 length "+toString(energy.size2()));
	if (i2>j2) throw std::runtime_error("PredictionTrackerSpotProbAll::updateProfile() : i2 "+toString(i2)+" > j2 "+toString(j2));
#endif

	// get Boltzmann weight contribution of this interaction
	Z_type curWeight = energy.getBoltzmannWeight( curE );

	// update overall partition function
	overallZ += curWeight;

	// update pair data
	for (size_t k1=i1; k1<=j1; k1++) {
		for (size_t k2=i2; k2<=j2; k2++) {
			// update partition function
			pairZ(k1,k2) += curWeight;
		}
	}
}


//////////////////////////////////////////////////////////////////////

void
PredictionTrackerSpotProbAll::
writeData( std::ostream &out
			, const Z2dMatrix & pairZ
			, const Z_type & overallZ
			, const InteractionEnergy & energy
			, const std::string & NA_string
			, const std::string & sep )
{
	// direct access to sequence string information
	const RnaSequence & rna1 = energy.getAccessibility1().getSequence();
	const RnaSequence & rna2 = energy.getAccessibility2().getAccessibilityOrigin().getSequence();
	const std::string & rna1Str = energy.getAccessibility1().getSequence().asString();
	const std::string & rna2Str = energy.getAccessibility2().getSequence().asString();

	// write in CSV-like format
	// header = reversed seq2
	// col[0] = seq1
	// cell[0,0] = "minE"

	// print header : spotProb; "nt_index" ... starting with index 1
	out <<"spotProb";
	for (size_t j=rna2.size(); j-- > 0; ) {
		out <<sep <<rna2Str.at(j)<<"_"<<rna2.getInOutIndex(energy.getAccessibility2().getReversedIndex(j));
	}
	out <<'\n';
	// print minE data
	for (size_t i=0; i<pairZ.size1(); i++) {
		// out nt in seq1 together with index
		out <<rna1Str.at(i)<<"_"<<rna1.getInOutIndex(i);
		for (size_t j=pairZ.size2(); j-- > 0; ) {
			// out separator
			out <<sep;
			// out infinity replacement if needed
			if ( Z_isINF( pairZ(i,j) ) ) {
				out<<NA_string;
			} else {
				// print probability = pairZ / overallZ
				out <<(pairZ(i,j)/overallZ);
			}
		}
		// line end
		out <<'\n';
	}
}

//////////////////////////////////////////////////////////////////////

} // namespace

