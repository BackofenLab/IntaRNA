
#include "SeedHandlerExplicit.h"

#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>


namespace IntaRNA
{

/////////////////////////////////////////////////////////////////////////////

SeedHandlerExplicit::
SeedHandlerExplicit(
		const InteractionEnergy & energy
		, const SeedConstraint & seedConstraint
		)
 :	SeedHandler(energy,seedConstraint)
	, seedForLeftEnd()
{
#if INTARNA_IN_DEBUG_MODE
	if (seedConstraint.getExplicitSeeds().empty()) throw std::runtime_error("SeedHandlerExplicit() : explicit seed encoding is empty");
#endif

	// parse individual seed encodings (comma separated)
	const std::string & seed = seedConstraint.getExplicitSeeds();
	size_t s = seed.find_first_not_of(',');
	size_t e = seed.find(',',s+1);
	while (s != std::string::npos) {
		// get length of current seed encoding
		const size_t length = (e == std::string::npos) ? seed.size()-s : e-s;
		// parse encoding
		SeedData seedData( seed.substr(s,length), getInteractionEnergy() );
		// warning if seed out of sequence length bounds (and ignore)
		if ( seedData.start1 == std::string::npos || seedData.start2 == std::string::npos ) {
			LOG(WARNING) <<"explicit seed '"<<seed.substr(s,length)<<"' is out of sequence boundaries and will be ignored";
		} else
		if ( ! seedData.allBpCompl ) {
			LOG(WARNING) <<"explicit seed '"<<seed.substr(s,length)<<"' contains non-complementary base pairs and will be ignored";
		} else
		if (!seedData.isValid()) {
			LOG(WARNING) <<"explicit seed '"<<seed.substr(s,length)<<"' not valid and will be ignored";
		} else {
			// generate hash key element for storing
			Interaction::BasePair leftBP = Interaction::BasePair(seedData.start1,seedData.start2);
			// check if left boundary already known
			if (seedForLeftEnd.find(leftBP) != seedForLeftEnd.end()) {
				// warn if multiple seeds for left boundary (currently not supported, only mfe saved)
				LOG(WARNING)<<"multiple explicit seed encodings for left seed end ("
						<<(seedData.start1+1)<<","<<(getInteractionEnergy().getAccessibility2().getReversedIndex(seedData.start2)+1)
						<<") : using only seed with lowest energy";
				// check for mfe seed
				if (seedForLeftEnd.at(leftBP).energy > seedData.energy) {
					// overwrite existing data, since current data has lower energy
					seedForLeftEnd[ leftBP ] = seedData;
				}
			} else {
				// store seed data
				seedForLeftEnd[ leftBP ] = seedData;
			}
		}
		// update boundaries for next seed encoding
		s = seed.find_first_not_of(',',e);
		e = seed.find(',',s+1);
	}

	// check if any seed available after parsing
	if (seedForLeftEnd.empty()) {
		throw std::runtime_error("SeedHandlerExplicit() : no valid explicit seed encodings found");
	}
}

/////////////////////////////////////////////////////////////////////////////

SeedHandlerExplicit::
~SeedHandlerExplicit()
{
}

/////////////////////////////////////////////////////////////////////////////

size_t
SeedHandlerExplicit::
getSeedMaxBP( const std::string & seedEncoding )
{
	size_t maxBP = 0;
	if (!seedEncoding.empty()) {
#if INTARNA_IN_DEBUG_MODE
		if (!checkSeedEncoding(seedEncoding).empty()) throw std::runtime_error("SeedHandlerExplicit::getSeedMaxBP() : no valid seed encoding : "+checkSeedEncoding(seedEncoding));
#endif
		size_t curBP = 0;
		for( const char & c : seedEncoding ) {
			switch(c){
			case ',' :
			case '&' : maxBP = std::max( maxBP, curBP ); curBP = 0; break;
			case '|' : curBP++;
			}
		}
	}
	return maxBP;
}

/////////////////////////////////////////////////////////////////////////////

std::string
SeedHandlerExplicit::
checkSeedEncoding( const std::string & seed )
{
	if (seed.empty()) {
		return "no seed encoding present";
	}

	// check general encoding
	if ( ! boost::regex_match(seed,boost::regex(R"(^-?\d+[|\.]+&-?\d+[|\.]+$)"), boost::match_perl) ) {
		return "is not of regular expression format start1dotbar1&start2dotbar2, e.g. '3|||.|&6||.||'";
	}

	// check dot-bar ends
	size_t pos = seed.find_first_of('|');
	if ( seed.at(pos-1) == '.' ) {
		return "seed dot-bar encoding has to start with a base pair '|', e.g. '3|||.|&6||.||'";
	}
	pos = seed.find_last_of('|', seed.find('&'));
	if ( seed.at(pos+1) == '.' ) {
		return "seed dot-bar encoding of the first part has to end with a base pair '|', e.g. '3|||.|&6||.||'";
	}
	pos = seed.find_first_of('|', seed.find('&'));
	if ( seed.at(pos-1) == '.' ) {
		return "seed dot-bar encoding of the second part has to start with a base pair '|', e.g. '3|||.|&6||.||'";
	}
	pos = seed.find_last_of('|');
	if ( pos+1<seed.size() ) {
		return "seed dot-bar encoding has to end with a base pair '|', e.g. '3|||.|&6||.||'";
	}

	// check balancing of base pair number
	int openBps = 0;
	bool behindAmp = true;
	for (const char si : seed ) {
		switch (si) {
		case '|' : openBps += (behindAmp ? +1 : -1); break;
		case '&' : behindAmp = false; break;
		default: break;
		}
	}
	if (openBps != 0) {
		return "number of base pairs '|' differ between sequences";
	}

	// all good
	return "";
}

//////////////////////////////////////////////////////////////////////////


size_t
SeedHandlerExplicit::
fillSeed( const size_t i1min, const size_t i1max, const size_t i2min, const size_t i2max)
{
	// check how many seeds are within the interval
	size_t withinIntervals = 0;
	for (auto it = seedForLeftEnd.begin(); it != seedForLeftEnd.end(); it++) {
		const SeedData & s = it->second;
		if (	s.start1 >= i1min
				&& s.start1+s.dotBar1.size()-1 <= i1max
				&& s.start2 >= i2min
				&& s.start2+s.dotBar2.size()-1 <= i2max
				)
		{
			withinIntervals++;
		}
	}
	// return identified number
	return withinIntervals;
}

//////////////////////////////////////////////////////////////////////////


void
SeedHandlerExplicit::
traceBackSeed( Interaction & interaction
		, const size_t i1
		, const size_t i2
		) const
{
	// copy according seed data if any
	// try to access seed information
	auto seed = seedForLeftEnd.find( Interaction::BasePair(i1,i2) );
#if INTARNA_IN_DEBUG_MODE
	if (seed == seedForLeftEnd.end()) throw std::runtime_error("SeedHandlerExplicit::traceBackSeed("+toString(i1)+","+toString(i2)+") : no explicit seed available");
#endif
	// fill data structure for current seed
	const SeedData & s = seed->second;
	// check if seed can not contain at least three base pairs
	if (s.dotBar1.size() < 3 && s.dotBar2.size() < 3) {
		return;
	}
	// get positions of second base pair within seed (if any)
	size_t l1=s.dotBar1.find('|', 1),
			l2=s.dotBar2.find('|', 1);
	// ignore left-most seed base pair
	// add all remaining but the last (right-most) base pair
	while( l1 != std::string::npos && l1+1 < s.dotBar1.size() ) {
		// store seed base pair
		interaction.basePairs.push_back( getInteractionEnergy().getBasePair( i1+l1, i2+l2 ) );
		// update l1, l2
		l1 = s.dotBar1.find('|', l1+1);
		l2 = s.dotBar2.find('|', l2+1);
	}
}

//////////////////////////////////////////////////////////////////////////


E_type
SeedHandlerExplicit::
getSeedE( const size_t i1, const size_t i2 ) const
{
	// try to access seed information
	auto seed = seedForLeftEnd.find( Interaction::BasePair(i1,i2) );

	// check if found
	if ( seed != seedForLeftEnd.end() ) {
		// return according energy
		return seed->second.energy;
	}
	// not found -> mark as invalid seed position
	return E_INF;
}

//////////////////////////////////////////////////////////////////////////

inline
bool
SeedHandlerExplicit::
isSeedBound( const size_t i1, const size_t i2 ) const
{
	// search for seed entry in hash
	return seedForLeftEnd.find( Interaction::BasePair(i1,i2) ) != seedForLeftEnd.end();
}

//////////////////////////////////////////////////////////////////////////

size_t
SeedHandlerExplicit::
getSeedLength1( const size_t i1, const size_t i2 ) const
{
	// try to access seed information
	auto seed = seedForLeftEnd.find( Interaction::BasePair(i1,i2) );

	// check if found
	if ( seed != seedForLeftEnd.end() ) {
		// return according energy
		return seed->second.dotBar1.size();
	}
	// not found -> mark as invalid seed position
	return 0;
}

//////////////////////////////////////////////////////////////////////////


size_t
SeedHandlerExplicit::
getSeedLength2( const size_t i1, const size_t i2 ) const
{
	// try to access seed information
	auto seed = seedForLeftEnd.find( Interaction::BasePair(i1,i2) );

	// check if found
	if ( seed != seedForLeftEnd.end() ) {
		// return according energy
		return seed->second.dotBar2.size();
	}
	// not found -> mark as invalid seed position
	return 0;
}

//////////////////////////////////////////////////////////////////////////

SeedHandlerExplicit::
SeedData::
SeedData()
 :	start1(std::string::npos)
	, start2(std::string::npos)
	, dotBar1("")
	, dotBar2("")
	, allBpCompl(true)
	, energy(E_INF)
{
}

//////////////////////////////////////////////////////////////////////////

SeedHandlerExplicit::
SeedData::
SeedData( const std::string & seedEncoding, const InteractionEnergy & energyFunction )
 :	start1(std::string::npos)
	, start2(std::string::npos)
	, dotBar1("")
	, dotBar2("")
	, allBpCompl(true)
	, energy(E_INF)
{
#if INTARNA_IN_DEBUG_MODE
	std::string errMsg = checkSeedEncoding(seedEncoding);
	if (!errMsg.empty()) throw std::runtime_error("SeedHandlerExplicit::SeedData("+seedEncoding+") : seed encoding not valid : "+ errMsg);
#endif
	// get start 1
	size_t p1 = 0;
	size_t p2 = seedEncoding.find('|');
	// index shift to internal indexing
	try {
		start1 = energyFunction.getAccessibility1().getSequence().getIndex(boost::lexical_cast<long>(seedEncoding.substr(p1, p2-p1)));
	} catch (std::runtime_error & e) {}
	// get dotbar 1
	p1 = p2;
	p2 = seedEncoding.find_last_of('|', seedEncoding.find('&'));
	dotBar1 = seedEncoding.substr(p1, p2-p1+1);
	// get start 2
	p1 = p2+2; // skip ampersand
	p2 = seedEncoding.find('|',p1+1);
	// index shift to internal indexing
	try {
		start2 = energyFunction.getAccessibility2().getAccessibilityOrigin().getSequence().getIndex(boost::lexical_cast<long>(seedEncoding.substr(p1, p2-p1)));
	} catch (std::runtime_error & e) {}
	// get dotbar 2
	p1 = p2;
	p2 = seedEncoding.find_last_of('|');
	dotBar2 = seedEncoding.substr(p1, p2-p1+1);
	// correct start1
	if (start1 != std::string::npos && start1 >= energyFunction.getAccessibility1().getSequence().size()) {
		start1 = std::string::npos;
	}
	// correct start2
	if (start2 != std::string::npos) {
		if (start2 >= energyFunction.getAccessibility2().getSequence().size()) {
			start2 = std::string::npos;
		} else {
			// reverse index of seq2 right end of seed
			start2 = energyFunction.getAccessibility2().getReversedIndex(start2);
			// correct to left end
			if (start2+1 >= dotBar2.size()) {
				start2 = start2+1-dotBar2.size();
			} else {
				// seed exceeds sequence end
				start2 = std::string::npos;
			}
		}
	}

	// compute seed energy (excluding right-most base pair's contribution)
	if (start1 != std::string::npos && start2 != std::string::npos ) {
		// get interaction information
		Interaction inter(energyFunction.getAccessibility1().getSequence(), energyFunction.getAccessibility2().getSequence());
		size_t l1 = 0, l2 = 0;
		// get positions of second base pair within seed (if any)
		size_t r1=dotBar1.find('|', l1+1),
				r2=dotBar2.find('|', l2+1);
		// check if left bp complementary
		allBpCompl = energyFunction.areComplementary( start1+l1, start2+l2);
		// add energy contributions
		energy = 0.0;
		while( allBpCompl && r1 != std::string::npos ) {
			// check if right bp complementary
			allBpCompl = energyFunction.areComplementary( start1+r1, start2+r2);
			if (allBpCompl) {
				// add energy contribution
				energy += energyFunction.getE_interLeft( start1+l1, start1+r1, start2+l2, start2+r2 );
				// update base pair positions
				l1 = r1;
				l2 = r2;
				r1 = dotBar1.find('|', l1+1);
				r2 = dotBar2.find('|', l2+1);
			}
		}
		if( ! allBpCompl ) {
			energy = E_INF;
		}
	}
}

//////////////////////////////////////////////////////////////////////////

bool
SeedHandlerExplicit::
SeedData::
isValid() const
{
	return (start1 != std::string::npos && start2 != std::string::npos && allBpCompl);
}

//////////////////////////////////////////////////////////////////////////

bool
SeedHandlerExplicit::
updateToNextSeed( size_t & i1_out, size_t & i2_out
		, const size_t i1min, const size_t i1max
		, const size_t i2min, const size_t i2max
		) const
{
	// ensure we do have any seed
	if (seedForLeftEnd.empty()) {
		return false;
	}

	size_t i1 = i1_out, i2 = i2_out;
	// find true max value
	const size_t i1maxVal = std::min(energy.size1()-1,i1max)
				, i2maxVal = std::min(energy.size2()-1,i2max);

	// find current seed
	auto curSeedData = seedForLeftEnd.find( Interaction::BasePair(i1,i2) );
	// check if we have to provide first seed (out of bound or no seed start)
	if (curSeedData == seedForLeftEnd.end()) {
		curSeedData = seedForLeftEnd.begin();
	} else {
		// go to next seed
		curSeedData++;
	}

	// find next seed within range
	while (curSeedData != seedForLeftEnd.end()
			&& (curSeedData->first.first < i1min
				|| curSeedData->first.first > i1maxVal
				|| curSeedData->first.second < i2min
				|| curSeedData->first.second > i2maxVal
			))
	{
		curSeedData++;
	}
	// ensure we have a valid new seed
	if (curSeedData != seedForLeftEnd.end()) {
		// copy data
		i1_out = curSeedData->first.first;
		i2_out = curSeedData->first.second;
		return true;
	}

	return false;
}

//////////////////////////////////////////////////////////////////////////

} /* namespace IntaRNA */
