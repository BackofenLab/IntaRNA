
#include "IntaRNA/PredictorMfe.h"

#include <iostream>
#include <algorithm>

namespace IntaRNA {

////////////////////////////////////////////////////////////////////////////

PredictorMfe::PredictorMfe(
		const InteractionEnergy & energy
		, OutputHandler & output
		, PredictionTracker * predTracker
		)
	: Predictor(energy,output,predTracker)
	, mfeInteractions()
	, reportedInteractions()
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
initOptima( const OutputConstraint & outConstraint )
{
	// resize to the given number of interactions if overlapping reports allowed
	mfeInteractions.resize(outConstraint.reportOverlap!=OutputConstraint::ReportOverlap::OVERLAP_BOTH ? 1 : outConstraint.reportMax
			, Interaction(energy.getAccessibility1().getSequence()
					,energy.getAccessibility2().getAccessibilityOrigin().getSequence())
			);

	// init all interactions to be filled
	for (InteractionList::iterator i = mfeInteractions.begin(); i!= mfeInteractions.end(); i++) {
		// initialize global E minimum : reported interactions have to have energy below that value
		i->energy = outConstraint.maxE;
		// ensure it holds only the boundary
		if (i->basePairs.size()!=2) {
			i->basePairs.resize(2);
		}
		// reset boundary base pairs
		i->basePairs[0].first = RnaSequence::lastPos;
		i->basePairs[0].second = RnaSequence::lastPos;
		i->basePairs[1].first = RnaSequence::lastPos;
		i->basePairs[1].second = RnaSequence::lastPos;
	}

	// clear reported interaction ranges
	reportedInteractions.first.clear();
	reportedInteractions.second.clear();

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe::
updateOptima( const size_t i1, const size_t j1
		, const size_t i2, const size_t j2
		, const E_type interE
		, const bool isHybridE )
{
//	LOG(DEBUG) <<"PredictorMfe::updateOptima( "<<i1<<"-"<<j1<<", "<<i2<<"-"<<j2<<" , E = " <<interE<<" isHybridE="<<(isHybridE?"true":"false");

	// ignore invalid reports
	if (E_isINF(interE) || interE >= E_MAX) {
		return;
	}

	// check if nothing to be done
	if (mfeInteractions.size() == 0) {
		// report call if needed
		if (predTracker != NULL) {
			// get final energy of current interaction
			E_type curE = isHybridE ? energy.getE( i1,j1, i2,j2, interE ) : interE;
			if (E_isNotINF(curE)) {
				// inform about prediction
				predTracker->updateOptimumCalled( i1 + (i1==RnaSequence::lastPos ? 0 : energy.getOffset1())
												, j1 + (j1==RnaSequence::lastPos ? 0 : energy.getOffset1())
												, i2 + (i2==RnaSequence::lastPos ? 0 : energy.getOffset2())
												, j2 + (j2==RnaSequence::lastPos ? 0 : energy.getOffset2())
												, curE );
			}
		}
		return;
	}

	// get final energy of current interaction
	E_type curE = isHybridE ? energy.getE( i1,j1, i2,j2, interE ) : interE;
//	LOG(DEBUG) <<"energy( "<<i1<<"-"<<j1<<", "<<i2<<"-"<<j2<<" ) = "
//			<<interE <<" : total = "<<curE;
	// report call if needed
	if (predTracker != NULL && E_isNotINF(curE)) {
		// inform about prediction
		predTracker->updateOptimumCalled( i1 + (i1==RnaSequence::lastPos ? 0 : energy.getOffset1())
										, j1 + (j1==RnaSequence::lastPos ? 0 : energy.getOffset1())
										, i2 + (i2==RnaSequence::lastPos ? 0 : energy.getOffset2())
										, j2 + (j2==RnaSequence::lastPos ? 0 : energy.getOffset2())
										, curE );
	}


	// check if we have to care about insertion (curE <= worst E in list)
	if (curE > mfeInteractions.rbegin()->energy ) {
		return;
	}

	// create temporary interaction
	Interaction tmp( *(mfeInteractions.begin()) );
	tmp.energy = curE;
	tmp.basePairs.resize(2);
	tmp.basePairs[0] = energy.getBasePair(i1,i2);
	tmp.basePairs[1] = energy.getBasePair(j1,j2);
	// handle single bp interactions
	if (tmp.basePairs[0] == tmp.basePairs[1]) { tmp.basePairs.resize(1); }
	INTARNA_CLEANUP(tmp.seed);

	if (mfeInteractions.size() == 1) {
		if ( tmp < *(mfeInteractions.begin()) ) {
			// overwrite current minimum
			*(mfeInteractions.begin()) = tmp;
		}
	} else {

		// check for insertion position
		InteractionList::iterator insertPos = std::find_if_not( mfeInteractions.begin(), mfeInteractions.end(), [&](Interaction & i){return i < tmp;});

		if ( insertPos != mfeInteractions.end() && !( tmp == *insertPos )) {

			// check for equivalence with last element
			if (insertPos != mfeInteractions.begin()) {
				--insertPos;
				if (*insertPos == tmp) {
					// already within list; nothing else to do
					return;
				}
				++insertPos;
			}

			// replace last element with current interaction
			InteractionList::iterator lastInteraction = mfeInteractions.begin();
			std::advance(lastInteraction,mfeInteractions.size()-1);
			*lastInteraction = tmp;

			// relocate last element within list (no reallocation needed, just list-link-update)
			mfeInteractions.splice( insertPos, mfeInteractions, lastInteraction );

		}
	}
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe::
reportOptima( const OutputConstraint & outConstraint )
{
	// number of reported interactions
	size_t reported = 0;
	// get maximal report energy = mfe + deltaE + precisionEpsilon
	const E_type mfeDeltaE = (mfeInteractions.begin()->energy + outConstraint.deltaE);

	// clear reported interaction ranges
	reportedInteractions.first.clear();
	reportedInteractions.second.clear();

	// check if non-overlapping output is wanted
	if (outConstraint.reportOverlap!=OutputConstraint::ReportOverlap::OVERLAP_BOTH) {
		// check if mfe is worth reporting
		Interaction curBest = *mfeInteractions.begin();
		while( curBest.energy < outConstraint.maxE
				&& (curBest.energy < mfeDeltaE || E_equal(curBest.energy,mfeDeltaE))
				&& reported < outConstraint.reportMax )
		{
			// report current best
			// fill interaction with according base pairs
			traceBack( curBest, outConstraint );
			// report mfe interaction
			output.add( curBest, outConstraint );

			// store ranges to ensure non-overlapping of next best solution
			switch( outConstraint.reportOverlap ) {
			case OutputConstraint::ReportOverlap::OVERLAP_BOTH :
				break;
			case OutputConstraint::ReportOverlap::OVERLAP_SEQ1 :
				reportedInteractions.second.insert( IndexRange(energy.getIndex2(*curBest.basePairs.begin()),energy.getIndex2(*curBest.basePairs.rbegin())) );
				break;
			case OutputConstraint::ReportOverlap::OVERLAP_SEQ2 :
				reportedInteractions.first.insert( IndexRange(energy.getIndex1(*curBest.basePairs.begin()),energy.getIndex1(*curBest.basePairs.rbegin())) );
				break;
			case OutputConstraint::ReportOverlap::OVERLAP_NONE :
				reportedInteractions.first.insert( IndexRange(energy.getIndex1(*curBest.basePairs.begin()),energy.getIndex1(*curBest.basePairs.rbegin())) );
				reportedInteractions.second.insert( IndexRange(energy.getIndex2(*curBest.basePairs.begin()),energy.getIndex2(*curBest.basePairs.rbegin())) );
				break;
			}

			// count
			reported++;
			// get next best if not already enough found
			if (reported < outConstraint.reportMax) {
				// get next best
				getNextBest( curBest );
			}
		}

	} // non-overlapping interactions

	else // overlapping interactions are allowed
	{

		// report all (possibly overlapping) interactions with energy below 0
		assert(outConstraint.reportMax <= mfeInteractions.size());
		for (InteractionList::iterator i = mfeInteractions.begin();
				reported < outConstraint.reportMax
				&& i!= mfeInteractions.end(); i++)
		{
			// check if interaction is within allowed energy range
			if ( i->energy < outConstraint.maxE
					&& (i->energy < mfeDeltaE || E_equal(i->energy,mfeDeltaE))) {

				// fill mfe interaction with according base pairs
				traceBack( *i, outConstraint );
				// report mfe interaction
				output.add( *i, outConstraint );
				// count
				reported++;
			}
		}
	} // overlapping interactions

	// check if nothing was worth reporting but report was expected
	if (reported == 0 && outConstraint.reportMax > 0 ) {
		// TODO replace with no report and "no-interaction-message" in output destructor
		// -> no preferable interactions found !!!
		// replace mfeInteraction with no interaction
		mfeInteractions.begin()->clear();
		mfeInteractions.begin()->energy = 0.0;
		// report mfe interaction
		output.add( *(mfeInteractions.begin()), outConstraint );
	}

}

////////////////////////////////////////////////////////////////////////////


} // namespace
