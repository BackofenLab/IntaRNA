
#include "IntaRNA/PredictionTrackerSpotProb.h"

#include <boost/foreach.hpp>
#include <boost/regex.hpp>

namespace IntaRNA {

//////////////////////////////////////////////////////////////////////

const std::string PredictionTrackerSpotProb::str_spot ="[123456789]\\d*&[123456789]\\d*";
const boost::regex PredictionTrackerSpotProb::regexSpotString("\\s*|"+str_spot+"(,"+str_spot+")*");

//////////////////////////////////////////////////////////////////////

PredictionTrackerSpotProb::
PredictionTrackerSpotProb(
		const InteractionEnergy & energy
		, const std::string & spotString
		, const std::string & outStreamName
	)
 :	PredictionTracker()
	, energy(energy)
	, outStream(NULL)
	, deleteOutStream(true)
	, spots()
	, noSpotZ(0.0)
	, overallZ(0.0)
{
#if INTARNA_IN_DEBUG_MODE
	// check spot encoding via regex
	if (!boost::regex_match( spotString, regexSpotString, boost::match_perl )) {
		throw std::runtime_error("PredictionTrackerSpotProb(spots="+spotString+") does not match its encoding regular expression");
	}
#endif

	// open stream
	if (!outStreamName.empty()) {
		outStream = newOutputStream( outStreamName );
		if (outStream == NULL) {
			throw std::runtime_error("PredictionTrackerSpotProb() : could not open output stream '"+outStreamName+"' for writing");
		}
	}

	// parse spot encoding
	size_t startPos = 0, splitPos = spotString.size();
	// find split positions
	while (startPos != splitPos) {
		splitPos = spotString.find(',',startPos);
		// insert interval
		spots.push_back( Spot(spotString.substr(startPos,splitPos-(splitPos==std::string::npos?0:startPos))));
		// update start of next interval encoding to parse
		startPos = splitPos + (splitPos != std::string::npos ? 1 : 0);
	}
}

//////////////////////////////////////////////////////////////////////

PredictionTrackerSpotProb::
PredictionTrackerSpotProb(
		const InteractionEnergy & energy
		, const std::string & spotString
		, std::ostream & outStream
	)
 :	PredictionTracker()
	, energy(energy)
	, outStream(&outStream)
	, deleteOutStream(false)
	, spots()
	, noSpotZ(0.0)
	, overallZ(0.0)
{
#if INTARNA_IN_DEBUG_MODE
	// check spot encoding via regex
	if (!boost::regex_match( spotString, regexSpotString, boost::match_perl )) {
		throw std::runtime_error("PredictionTrackerSpotProb(spots="+spotString+") does not match its encoding regular expression");
	}
#endif

	// parse spot encoding
	size_t startPos = 0, splitPos = spotString.size();
	// find split positions
	while (startPos != splitPos) {
		splitPos = spotString.find(',',startPos);
		// insert interval
		spots.push_back( Spot(spotString.substr(startPos,splitPos-(splitPos==std::string::npos?0:startPos))));
		// update start of next interval encoding to parse
		startPos = splitPos + (splitPos != std::string::npos ? 1 : 0);
	}
}

//////////////////////////////////////////////////////////////////////

PredictionTrackerSpotProb::
~PredictionTrackerSpotProb()
{
	// write probabilities to streams
	// probability of interactions covering no tracked spot
	(*outStream) <<"spot;probability\n";
	// handle division by zero if nothing was reported
	if (Z_equal(overallZ,0.0)) {
		overallZ = 1.0;
	}
	// probability of interactions covering no tracked spot
	(*outStream) <<"0&0;"<<(noSpotZ / overallZ)<<'\n';
	// probabilities of tracked spots
	BOOST_FOREACH( Spot & s, spots) {
		(*outStream) <<(s.idx1+1)<<'&'<<(s.idx2+1)<<';'<<(s.Z/overallZ)<<'\n';
	}
	outStream->flush();

	if (deleteOutStream) {
		// clean up if file pointers were created in constructor
		deleteOutputStream( outStream );
	}
}

//////////////////////////////////////////////////////////////////////

void
PredictionTrackerSpotProb::
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

	// get Boltzmann weight of this interaction
	const Z_type curBW = energy.getBoltzmannWeight( curE );
	// update overall Z
	overallZ += curBW;
	// update spot information
	bool noSpotCovered = true;

	// update final output handler
	BOOST_FOREACH( Spot & s, spots) {
		// check if covered
		if ( bp_l.first <= s.idx1 && s.idx1 <= bp_r.first
			&& bp_r.second <= s.idx2 && s.idx2 <= bp_l.second)
		{
			// update partition function of this spot, since it is covered
			s.Z += curBW;
			// note that a spot was covered
			noSpotCovered = false;
		}
	}

	// check if any spot was covered
	if (noSpotCovered) {
		// update noSpotZ
		noSpotZ += curBW;
	}
}


////////////////////////////////////////////////////////////////////////////


} // namespace
