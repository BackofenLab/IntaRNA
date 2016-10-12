
#include "PredictorMfe.h"

#include <iostream>

////////////////////////////////////////////////////////////////////////////

PredictorMfe::PredictorMfe( const InteractionEnergy & energy, OutputHandler & output )
	: Predictor(energy,output)
	, mfeInteraction(energy.getAccessibility1().getSequence()
			,energy.getAccessibility2().getAccessibilityOrigin().getSequence())
	, i1offset(0)
	, i2offset(0)
	, minStackingEnergy( energy.getBestE_interLoop() )
	, minInitEnergy( energy.getE_init() )
	, minDangleEnergy( energy.getBestE_dangling() )
	, minEndEnergy( energy.getBestE_end() )

{

}

////////////////////////////////////////////////////////////////////////////

PredictorMfe::~PredictorMfe()
{
}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfe::
initMfe()
{
	// initialize global E minimum : should be below 0.0
	mfeInteraction.energy = 0.0;
	// ensure it holds only the boundary
	if (mfeInteraction.basePairs.size()!=2) {
		mfeInteraction.basePairs.resize(2);
	}
	// reset boundary base pairs
	mfeInteraction.basePairs[0].first = RnaSequence::lastPos;
	mfeInteraction.basePairs[0].second = RnaSequence::lastPos;
	mfeInteraction.basePairs[1].first = RnaSequence::lastPos;
	mfeInteraction.basePairs[1].second = RnaSequence::lastPos;
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe::
updateMfe( const size_t i1, const size_t j1
		, const size_t i2, const size_t j2
		, const E_type hybridE )
{
//	std::cerr <<"#DEBUG : energy( "<<i1<<"-"<<j1<<", "<<i2<<"-"<<j2<<" ) = "
//			<<hybridE
//			<<std::endl;

	// TODO check if reasonable to check only interactions with hybridE < 0
	if (hybridE > 0) {
		return;
	}

	// get final energy of current interaction
	E_type curE = energy.getE( i1,j1, i2,j2, hybridE );
//	std::cerr <<"#DEBUG : energy( "<<i1<<"-"<<j1<<", "<<i2<<"-"<<j2<<" ) = "
//			<<hybridE <<" : total = "<<curE
//			<<std::endl;

	if (curE < mfeInteraction.energy) {
//	std::cerr <<"#DEBUG : new mfe( "<<i1<<"-"<<j1<<", "<<i2<<"-"<<j2<<" ) = "
//			<<hybridE
//			<<std::endl;
		// store new global min
		mfeInteraction.energy = (curE);
		// store interaction boundaries
		// left
		mfeInteraction.basePairs[0].first = i1+i1offset;
		mfeInteraction.basePairs[0].second = energy.getAccessibility2().getReversedIndex(i2+i2offset);
		// right
		mfeInteraction.basePairs[1].first = j1+i1offset;
		mfeInteraction.basePairs[1].second = energy.getAccessibility2().getReversedIndex(j2+i2offset);
	}
}

////////////////////////////////////////////////////////////////////////////

