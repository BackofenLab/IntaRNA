
#include "PredictionTrackerPairMinE.h"

//////////////////////////////////////////////////////////////////////

PredictionTrackerPairMinE::
PredictionTrackerPairMinE(
		const InteractionEnergy & energy
		, const std::string & streamName
		, const std::string E_INF_string
	)
 :	PredictionTracker()
	, energy(energy)
	, deleteStreamsOnDestruction(true)
	, outStream(NULL)
	, E_INF_string(E_INF_string)
	, pairMinE( energy.size1(), energy.size2(), E_INF ) // init E_INF
{
#if IN_DEBUG_MODE
	if (streamName.empty()) {
		throw std::runtime_error("PredictionTrackerPairMinE() : streamName empty");
	}
#endif
	// open streams
	outStream = newOutputStream( streamName );
	if (outStream == NULL) {
		throw std::runtime_error("PredictionTrackerPairMinE() : could not open output stream '"+streamName+"' for writing");
	}
}

//////////////////////////////////////////////////////////////////////

PredictionTrackerPairMinE::
PredictionTrackerPairMinE(
		const InteractionEnergy & energy
		, std::ostream * outStream
		, const std::string E_INF_string
	)
 :	PredictionTracker()
	, energy(energy)
	, deleteStreamsOnDestruction(false)
	, outStream(outStream)
	, E_INF_string(E_INF_string)
	, pairMinE( energy.size1(), energy.size2(), E_INF ) // init E_INF
{
}

//////////////////////////////////////////////////////////////////////

PredictionTrackerPairMinE::
~PredictionTrackerPairMinE()
{
	writeData( *outStream
				, pairMinE
				, energy
				, E_INF_string );

	// clean up if file pointers were created in constructor
	if (deleteStreamsOnDestruction) {
		deleteOutputStream( outStream );
	}
}

//////////////////////////////////////////////////////////////////////

void
PredictionTrackerPairMinE::
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

#if IN_DEBUG_MODE
	// sanity checks
	if (i1>=energy.size1() || j1>=energy.size1()) throw std::runtime_error("PredictionTrackerPairMinE::updateProfile() : index-1 range ["+toString(i1)+","+toString(j1)+"] exceeds sequence-1 length "+toString(energy.size1()));
	if (i1>j1) throw std::runtime_error("PredictionTrackerPairMinE::updateProfile() : i1 "+toString(i1)+" > j1 "+toString(j1));
	if (i2>=energy.size2() || j2>=energy.size2()) throw std::runtime_error("PredictionTrackerPairMinE::updateProfile() : index-2 range ["+toString(i2)+","+toString(j2)+"] exceeds sequence-2 length "+toString(energy.size2()));
	if (i2>j2) throw std::runtime_error("PredictionTrackerPairMinE::updateProfile() : i2 "+toString(i2)+" > j2 "+toString(j2));
#endif

	// update pair data
	for (size_t k1=i1; k1<=j1; k1++) {
		for (size_t k2=i2; k2<=j2; k2++) {
			// check if E is smaller than current minE
			if ( curE < pairMinE(k1,k2) ) {
				// write new minimum
				pairMinE(k1,k2) = curE;
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////

void
PredictionTrackerPairMinE::
writeData( std::ostream &out
			, const E2dMatrix & pairMinE
			, const InteractionEnergy & energy
			, const std::string & E_INF_string )
{
	// direct access to sequence string information
	const std::string & rna1 = energy.getAccessibility1().getSequence().asString();
	const std::string & rna2 = energy.getAccessibility2().getSequence().asString();

	// write in CSV-like format
	// header = reversed seq2
	// col[0] = seq1
	// cell[0,0] = "minE"

	// print header : minE;rSeq2[0],rSeq2[1],...
	out <<"minE";
	for (auto nt = rna1.begin(); nt!=rna1.end(); nt++) {
		out <<';' <<(*nt);
	}
	out <<'\n';
	// print minE data
	for (size_t i=0; i<pairMinE.size1(); i++) {
		// out nt in seq1
		out <<rna2.at(i);
		for (size_t j=0; j<pairMinE.size2(); j++) {
			// out separator
			out <<';';
			// out infinity replacement if needed
			if ( E_isINF( pairMinE(i,j) ) ) {
				out<<E_INF_string;
			} else {
				// print energy
				out <<pairMinE(i,j);
			}
		}
		// line end
		out <<'\n';
	}
}

//////////////////////////////////////////////////////////////////////
