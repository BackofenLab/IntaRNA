
#include "IntaRNA/PredictorMfe.h"

#include <iostream>
#include <algorithm>
#include <typeinfo>

namespace IntaRNA {
namespace {

/** Borrow one already computed seed, while retaining all original annotations. */
class SelectedSeed : public SeedHandlerIdxOffset {
public:
	SelectedSeed(SeedHandler & original, const Interaction::BasePair & selected);
	size_t fillSeed(size_t i, size_t j, size_t k, size_t l) override;
	bool isSeedBound(size_t i, size_t j) const override;
	E_type getSeedE(size_t i, size_t j) const override;
	size_t getSeedLength1(size_t i, size_t j) const override;
	size_t getSeedLength2(size_t i, size_t j) const override;
	bool updateToNextSeed(size_t & i, size_t & j, size_t lo1, size_t hi1, size_t lo2, size_t hi2) const override;
private:
	const Interaction::BasePair selected;
};

inline SelectedSeed::SelectedSeed(SeedHandler & original, const Interaction::BasePair & selected)
	: SeedHandlerIdxOffset(&original,false), selected(selected)
{}

inline bool SelectedSeed::isSeedBound(size_t i, size_t j) const
{ return selected == Interaction::BasePair(i,j) && SeedHandlerIdxOffset::isSeedBound(i,j); }
inline E_type SelectedSeed::getSeedE(size_t i, size_t j) const
{ return isSeedBound(i,j) ? SeedHandlerIdxOffset::getSeedE(i,j) : E_INF; }
inline size_t SelectedSeed::getSeedLength1(size_t i, size_t j) const
{ return isSeedBound(i,j) ? SeedHandlerIdxOffset::getSeedLength1(i,j) : 0; }
inline size_t SelectedSeed::getSeedLength2(size_t i, size_t j) const
{ return isSeedBound(i,j) ? SeedHandlerIdxOffset::getSeedLength2(i,j) : 0; }
inline size_t SelectedSeed::fillSeed(size_t i, size_t j, size_t k, size_t l)
{ return selected.first >= i && selected.first <= j && selected.second >= k && selected.second <= l && isSeedBound(selected.first,selected.second); }
inline bool SelectedSeed::updateToNextSeed(size_t & i, size_t & j, size_t lo1, size_t hi1, size_t lo2, size_t hi2) const
{
	if (Interaction::BasePair(i,j) == selected || (selected.first < lo1 || selected.first > hi1 || selected.second < lo2 || selected.second > hi2 || !isSeedBound(selected.first,selected.second))) return false;
	i = selected.first; j = selected.second;
	return true;
}

/** Forward callbacks without transferring ownership of the caller's tracker. */
class BorrowedTracker : public PredictionTracker {
public:
	explicit BorrowedTracker(PredictionTracker & tracker);
	void updateOptimumCalled(size_t i, size_t j, size_t k, size_t l, E_type e) override;
private:
	PredictionTracker & tracker;
};
inline BorrowedTracker::BorrowedTracker(PredictionTracker & tracker) : tracker(tracker) {}
inline void BorrowedTracker::updateOptimumCalled(size_t i,size_t j,size_t k,size_t l,E_type e)
{ tracker.updateOptimumCalled(i,j,k,l,e); }

/** Discard ordinary ranking while collecting its unchanged partition total. */
class StatisticsOutput : public OutputHandler {
public:
	explicit StatisticsOutput(const OutputConstraint & constraint);
	void add(const Interaction &) override;
};
inline StatisticsOutput::StatisticsOutput(const OutputConstraint & constraint) : OutputHandler(constraint) {}
inline void StatisticsOutput::add(const Interaction &) {}

/** Tag the independently ranked candidates with their generating seed. */
class SeedOutput : public OutputHandler {
public:
	SeedOutput(const OutputConstraint & constraint, OutputHandlerSeedCandidates & destination, Interaction::BasePair seed);
	void add(const Interaction & interaction) override;
private:
	OutputHandlerSeedCandidates & destination;
	const Interaction::BasePair seed;
};
inline SeedOutput::SeedOutput(const OutputConstraint & constraint, OutputHandlerSeedCandidates & destination, Interaction::BasePair seed)
	: OutputHandler(constraint), destination(destination), seed(seed) {}
inline void SeedOutput::add(const Interaction & interaction) { destination.addSeed(seed,interaction); }

} // namespace

void PredictorMfe::initSeedCandidates()
{
	ownedSeedCandidates.reset();
	seedCandidates = nullptr;
	if (!output.getOutputConstraint().reportMaxPerSeed) return;
	seedCandidates = dynamic_cast<OutputHandlerSeedCandidates *>(&output);
	if (!seedCandidates) {
		ownedSeedCandidates = std::make_unique<OutputHandlerSeedCandidates>(output.getOutputConstraint());
		seedCandidates = ownedSeedCandidates.get();
	}
}

bool PredictorMfe::predictPerSeed(const IndexRange & r1, const IndexRange & r2,
		SeedHandlerIdxOffset & seeds,
		const std::function<Predictor *(OutputHandler &, PredictionTracker *, SeedHandler *)> & factory)
{
	const auto & c = output.getOutputConstraint();
	if (!c.reportMaxPerSeed) return false;
	// Validate before unsigned range arithmetic or seed-table access.
	energy.setOffset1(0); energy.setOffset2(0);
	const size_t n = energy.size1(), m = energy.size2();
	if (!r1.isAscending() || !r2.isAscending() || r1.from >= n || r2.from >= m)
		throw std::invalid_argument("per-seed prediction: invalid sequence range");
	initSeedCandidates();
	SeedHandler & original = seeds.getOriginalSeedHandler();
	// A custom recurrence must explicitly supply its own child factory rather
	// than silently reverting to a standard predictor when quotas are enabled.
	auto create = [&](OutputHandler & out, PredictionTracker * tracker, SeedHandler * handler) {
		std::unique_ptr<Predictor> child(factory(out,tracker,handler));
		if (typeid(*child) != typeid(*this))
			throw std::logic_error("per-seed prediction: custom predictor must override predictPerSeed with its own factory");
		return child;
	};
	// Compute statistics once with every seed and the original prediction mode.
	Zall = 0;
	if (c.needZall || predTracker) {
		const OutputConstraint ordinary(c.reportMax,c.reportOverlap,c.maxE,c.deltaE,c.bestSeedOnly,
				c.noLP,c.noGUend,c.needZall,c.needBPs,c.maxED);
		StatisticsOutput statistics(ordinary);
		auto predictor = create(statistics,
				predTracker ? new BorrowedTracker(*predTracker) : nullptr,
				new SeedHandlerIdxOffset(&original,false));
		predictor->predict(r1,r2);
		Zall = predictor->getZall();
		if (statistics.usesExactPartition()) output.setExactPartition(true);
		output.incrementZ(statistics.getZ());
	}
	// Reuse the admitted seed table; singleton adapters never refill or own it.
	const size_t last1 = std::min(r1.to,n-1), last2 = std::min(r2.to,m-1);
	if (c.reportMax && original.fillSeed(r1.from,last1,r2.from,last2)) {
		const OutputConstraint perSeed(c.reportMaxPerSeed,OutputConstraint::OVERLAP_BOTH,E_INF,E_INF,
				c.bestSeedOnly,c.noLP,c.noGUend,false,c.needBPs,c.maxED);
		size_t i = RnaSequence::lastPos, j = RnaSequence::lastPos;
		while (original.updateToNextSeed(i,j,r1.from,last1,r2.from,last2)) {
			const size_t l1 = original.getSeedLength1(i,j), l2 = original.getSeedLength2(i,j);
			if (!l1 || !l2 || l1-1 > last1-i || l2-1 > last2-j
					|| l1 > energy.getAccessibility1().getMaxLength() || l2 > energy.getAccessibility2().getMaxLength()) continue;
			SeedOutput ranked(perSeed,*seedCandidates,energy.getBasePair(i,j));
			auto predictor = create(ranked,nullptr,new SelectedSeed(original,{i,j}));
			// Legacy recurrences may visit cells beyond the output span. Exclude
			// those before they can spend a seed's quota, retaining zero behavior.
			static_cast<PredictorMfe *>(predictor.get())->perSeedPass = true;
			predictor->predict(r1,r2);
		}
	}
	if (ownedSeedCandidates) seedCandidates->report(output);
	return true;
}

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
	initSeedCandidates();
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

bool
PredictorMfe::
isValidOutputSite( const size_t i1, const size_t j1,
		const size_t i2, const size_t j2 ) const
{
	const OutputConstraint & constraint = output.getOutputConstraint();
	if (perSeedPass && (j1-i1+1 > energy.getAccessibility1().getMaxLength()
			|| j2-i2+1 > energy.getAccessibility2().getMaxLength())) return false;
	if (constraint.noGUend && (energy.isGU(i1,i2) || energy.isGU(j1,j2))) {
		return false;
	}
	return constraint.maxED >= Accessibility::ED_UPPER_BOUND
			|| (energy.getED1(i1,j1) <= constraint.maxED
					&& energy.getED2(i2,j2) <= constraint.maxED);
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

	if (!isValidOutputSite(i1, j1, i2, j2)) {
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
	if (output.getOutputConstraint().reportMax <= 1
		&& curE > mfeInteractions.rbegin()->energy )
	{
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
	// Ensemble candidates were ranked by anchor before shared site pruning.
	if (seedCandidates) {
		if (ownedSeedCandidates) seedCandidates->report(output);
		return;
	}

	// A zero ranked-output count still publishes ensemble statistics.
	if (mfeInteractions.empty()) return;
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
			// Selectors return two boundaries, which coincide for a single pair.
			// Normalize before either traceback or boundary-only output validation.
			if (curBest.basePairs.size() == 2 && curBest.basePairs.front() == curBest.basePairs.back()) {
				curBest.basePairs.resize(1);
			}
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
