
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
	, mfe4leftEnd()
	, reportedInteractions()
{

}

////////////////////////////////////////////////////////////////////////////

PredictorMfe::~PredictorMfe()
{
}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfe::
initOptima()
{
	// temporary access
	const OutputConstraint & outConstraint = output.getOutputConstraint();

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

	// clear mfe information for left boundaries
	mfe4leftEnd.clear();

	// init overallZ if needed
	if (outConstraint.needZall) {
		Zall = 0.0;
	}
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe::
updateOptima( const size_t i1, const size_t j1
		, const size_t i2, const size_t j2
		, const E_type interE
		, const bool isHybridE
		, const bool incrementZall )
{
//	LOG_IF(energy.getBasePair(i1,i2)==Interaction::BasePair(240-1,88-1) && energy.getBasePair(j1,j2)==Interaction::BasePair(268-1,59-1),DEBUG) <<" found it! E="<<interE;
	// ignore invalid reports
	if (E_isINF(interE) || interE >= E_MAX) {
		return;
	}

	// check GU ends if needed
	if (output.getOutputConstraint().noGUend && (energy.isGU(i1,i2) || energy.isGU(j1,j2)) ) {
		return;
	}

	// update Zall if needed
	if (incrementZall) {
		updateZall( i1,j1, i2,j2, interE, isHybridE );
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
		///////  handle non-overlapping suboptimal output  ////////////
		if (output.getOutputConstraint().reportOverlap != OutputConstraint::OVERLAP_BOTH)
		{
			// update mfe4leftEnd information
			updateMfe4leftEnd(i1,j1,i2,j2,tmp);
		}
	} else
	///////  handle fully overlapping suboptimal output  ////////////
	if (output.getOutputConstraint().reportOverlap == OutputConstraint::OVERLAP_BOTH) {

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

	} // >1 overlapping suboptimals
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe::
reportOptima()
{
	// store overall partition function
	if (Z_isNotINF(getZall())) {
		output.incrementZ( getZall() );
	}

	// temporary access
	const OutputConstraint & outConstraint = output.getOutputConstraint();
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
			if (outConstraint.needBPs) {
				// fill interaction with according base pairs
				traceBack( curBest );
			}
			// report mfe interaction
#if INTARNA_MULITHREADING
			#pragma omp critical(intarna_omp_predictorOutputAdd)
#endif
			{output.add( curBest );}

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

		// report all (possibly overlapping) interactions with energy below threshold
		assert(outConstraint.reportMax <= mfeInteractions.size());
		for (InteractionList::iterator i = mfeInteractions.begin();
				reported < outConstraint.reportMax
				&& i!= mfeInteractions.end(); i++)
		{
			// check if interaction is within allowed energy range
			if ( i->energy < outConstraint.maxE
					&& (i->energy < mfeDeltaE || E_equal(i->energy,mfeDeltaE))) {

				if (outConstraint.needBPs) {
					// fill mfe interaction with according base pairs
					traceBack( *i );
				}
				// report mfe interaction
#if INTARNA_MULITHREADING
				#pragma omp critical(intarna_omp_predictorOutputAdd)
#endif
				{output.add( *i );}
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
#if INTARNA_MULITHREADING
		#pragma omp critical(intarna_omp_predictorOutputAdd)
#endif
		{output.add( *(mfeInteractions.begin()) );}
	}

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe::
getNextBest( Interaction & curBest )
{

	// get original
	const E_type curBestE = curBest.energy;
	// clear seed information if present
	if (curBest.seed != NULL) {
		curBest.seed->clear();
		INTARNA_CLEANUP(curBest.seed);
	}

	// identify cell with next best non-overlapping interaction site
	// iterate (decreasingly) over all left interaction starts
	BestInteractionE * curBestCell = NULL;
	E_type curBestCellE = E_INF;
	HashIdx2E::key_type curBestCellStart;
	BestInteractionE * curCell = NULL;
	E_type curCellE = E_INF;
	for (auto curLeftBound = mfe4leftEnd.begin(); curLeftBound != mfe4leftEnd.end(); curLeftBound++) {

		// interaction is non-overlapping with any already reported interaction
		if (   reportedInteractions.first.overlaps( IndexRange(curLeftBound->first.first, curLeftBound->second.j1) )
			|| reportedInteractions.second.overlaps( IndexRange(curLeftBound->first.second, curLeftBound->second.j2) )
				)
		{
			continue;
		}
		// pointer access to cell value
		curCell = &(curLeftBound->second);
		// get overall energy of the interaction
		curCellE = curLeftBound->second.val;
		// deal with degeneracy of energy model, ie. multiple sites with bestE
		if (E_equal(curCellE,curBestCellE)) {
			// select left most to be deterministic in output
			if (curLeftBound->first.first > curBestCellStart.first
					|| (curLeftBound->first.first == curBestCellStart.first
							&& curLeftBound->first.second > curBestCellStart.second)
					|| (curLeftBound->first.first == curBestCellStart.first
							&& curLeftBound->first.second == curBestCellStart.second
							&& curLeftBound->second.j1 > curBestCell->j1)
					|| (curLeftBound->first.first == curBestCellStart.first
							&& curLeftBound->first.second == curBestCellStart.second
							&& curLeftBound->second.j1 == curBestCell->j1
							&& curLeftBound->second.j2 > (curBestCell->j2))
				)
			{
				// right of curBest
				continue;
			}
		} else
		// or energy is too low to be considered
		// or energy is higher than current best found so far
		if (curCellE < curBestE || curCellE > curBestCellE)
		{
			continue;
		}

		//// FOUND THE NEXT BETTER SOLUTION
		// overwrite current best found so far
		curBestCell = curCell;
		curBestCellE = curCellE;
		curBestCellStart = curLeftBound->first;

	} // all left boundaries of valid interactions

	// overwrite curBest
	curBest.energy = curBestCellE;
	curBest.basePairs.resize(2);
	if (E_isNotINF(curBestCellE)) {
		curBest.basePairs[0] = energy.getBasePair( curBestCellStart.first, curBestCellStart.second );
		curBest.basePairs[1] = energy.getBasePair( curBestCell->j1, curBestCell->j2 );
	}
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe::
updateMfe4leftEnd(const size_t i1, const size_t j1, const size_t i2, const size_t j2, const Interaction & tmp )
{
	auto mfe4leftEndEntry = mfe4leftEnd.find(HashIdx2E::key_type(i1,i2));
	// check if left end is unknown -> just add
	if (mfe4leftEndEntry == mfe4leftEnd.end()) {
		// store current interaction
		mfe4leftEnd[HashIdx2E::key_type(i1,i2)] = BestInteractionE(tmp.energy,j1,j2);
	} else {
	// else: update entry if needed
		// create temporary interaction for comparison
		Interaction mfe4curLeftBound(tmp);
		mfe4curLeftBound.energy = mfe4leftEndEntry->second.val;
		*(mfe4curLeftBound.basePairs.rbegin()) = energy.getBasePair(mfe4leftEndEntry->second.j1,mfe4leftEndEntry->second.j2);
		// check if current interaction is to be stored
		if ( tmp < mfe4curLeftBound ) {
			// overwrite current minimum
			mfe4leftEndEntry->second.val = tmp.energy;
			mfe4leftEndEntry->second.j1 = j1;
			mfe4leftEndEntry->second.j2 = j2;
		}
	}
}

////////////////////////////////////////////////////////////////////////////


} // namespace
