
#ifndef INTARNA_PREDICTORMFE_H_
#define INTARNA_PREDICTORMFE_H_


#include "IntaRNA/Predictor.h"

#include "IntaRNA/IndexRangeList.h"
#include "IntaRNA/OutputHandlerSeedCandidates.h"
#include "IntaRNA/SeedHandlerIdxOffset.h"

#include <functional>
#include <list>
#include <memory>
#include <utility>

#include <boost/unordered_map.hpp>

namespace IntaRNA {

/**
 * Generic Predictor interface for MFE interaction computation to avoid
 * code redundancy
 *
 * @author Martin Mann
 *
 */
class PredictorMfe : public Predictor {

protected:
	/** Predict independently with each admitted seed, preserving this mode's
	 * recurrences and qualifying traceback. Statistics/trackers use one ordinary
	 * prediction; seed-specific passes do not contribute to them.
	 * @param r1 target region in original coordinates
	 * @param r2 reversed-query region in original coordinates
	 * @param seeds owned seed wrapper of this predictor
	 * @param factory constructs this predictor type, owning the supplied tracker
	 *            and seed adapter; children have per-seed selection disabled
	 * @return true if positive-quota prediction was handled completely
	 * @throws std::logic_error if a custom predictor subclass does not override
	 *            this hook to supply a factory preserving its custom recurrences
	 */
	virtual bool predictPerSeed(const IndexRange & r1, const IndexRange & r2,
			SeedHandlerIdxOffset & seeds,
			const std::function<Predictor *(OutputHandler &, PredictionTracker *, SeedHandler *)> & factory);

	/** Initialize local selection, or use a caller's shared regional collector. */
	void initSeedCandidates();
	//! Active collector; null on the unchanged default path.
	OutputHandlerSeedCandidates * seedCandidates = nullptr;
	//! Local collector for API callers that do not provide a regional collector.
	std::unique_ptr<OutputHandlerSeedCandidates> ownedSeedCandidates;
	//! Enforce complete-site spans before ranking a seed-constrained child pass.
	bool perSeedPass = false;

	/**
	 * Describes the currently best interaction found for a left interaction
	 * boundary i1,i2
	 */
	template < class ValueType >
	class BestInteraction {
	public:

		/**
		 * Init object
		 * @param val the value to be stored (default = INF [or MAX if not supported])
		 * @param j1 first index
		 * @param j1 second index
		 */
		BestInteraction( const ValueType val=(std::numeric_limits<ValueType>::has_infinity ? std::numeric_limits<ValueType>::infinity() : std::numeric_limits<ValueType>::max())
				, const size_t j1=RnaSequence::lastPos, const size_t j2=RnaSequence::lastPos )
			: val(val), j1(j1), j2(j2)
		{}

	public:
		//! value to be stored for the interaction, e.g. energy
		ValueType val;
		//! right end of the interaction in seq1
		size_t j1;
		//! right end of the interaction in seq2
		size_t j2;
	};

	//! BestInteraction that stores an energy value
	typedef BestInteraction<E_type> BestInteractionE;
	//! BestInteraction that stores a partition function value
	typedef BestInteraction<Z_type> BestInteractionZ;



public:

	/**
	 * Construction call the super constructor
	 *
	 * @param energy the interaction energy handler
	 * @param output the output handler to report optimal interactions to
	 * @param predTracker the prediction tracker to be used or NULL if no
	 *         tracking is to be done; if non-NULL, the tracker gets deleted
	 *         on this->destruction.
	 */
	PredictorMfe( const InteractionEnergy & energy
				, OutputHandler & output
				, PredictionTracker * predTracker );

	virtual ~PredictorMfe() override;

protected:


	//! access to the interaction energy handler of the super class
	using Predictor::energy;

	//! access to the output handler of the super class
	using Predictor::output;

	//! access to the prediction tracker of the super class
	using Predictor::predTracker;

	//! list of interactions
	typedef std::list<Interaction> InteractionList;

	//! mfe interaction boundaries
	InteractionList mfeInteractions;

	//! hash to map index pairs to BestInteractionE entries
	typedef boost::unordered_map< Interaction::BasePair, BestInteractionE, Interaction::BasePair::Hash > HashIdx2E;

	//! if non-overlapping output is required, this data structure is filled
	//! to find non-overlapping interactions
	HashIdx2E mfe4leftEnd;

	//! index ranges of reported interactions to identify non-overlapping
	//! interactions (first = seq1, second = seq2)
	//! NOTE: the indices for seq2 are reversed
	std::pair< IndexRangeList, IndexRangeList > reportedInteractions;

	/**
	 * Initializes the global energy minimum storage
	 */
	virtual
	void
	initOptima() override;

	/**
	 * Checks terminal GU and accessibility constraints for a complete site.
	 * Recursion cells may be useful extensions without being valid output sites.
	 * @param i1 inclusive start in sequence 1, relative to the energy offset
	 * @param j1 inclusive end in sequence 1, relative to the energy offset
	 * @param i2 inclusive start in reversed sequence 2, relative to its offset
	 * @param j2 inclusive end in reversed sequence 2, relative to its offset
	 * @return whether both terminal pairs and ED penalties satisfy output limits
	 */
	bool
	isValidOutputSite( const size_t i1, const size_t j1,
			const size_t i2, const size_t j2 ) const;

	/**
	 * updates the global optimum to be the mfe interaction if needed
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1
	 * @param energy the energy of the interaction
	 * @param isHybridE whether or not the given energy is only the
	 *        hybridization energy (init+loops) or the total interaction energy
	 * @param incrementZall whether or not Zall is to be incremented (if needed)
	 */
	virtual
	void
	updateOptima( const size_t i1, const size_t j1
				, const size_t i2, const size_t j2
				, const E_type energy
				, const bool isHybridE
				, const bool incrementZall ) override;


	/**
	 * Fills a given interaction (boundaries given) with the according
	 * hybridizing base pairs.
	 * Note, the
	 * @param interaction IN/OUT the interaction to fill
	 */
	virtual
	void
	traceBack( Interaction & interaction ) = 0;


	/**
	 * Identifies the next best interaction with an energy equal to or higher
	 * than the given interaction. The new interaction will not overlap any
	 * index range stored in reportedInteractions.
	 *
	 * @param curBest IN/OUT the current best interaction to be replaced with one
	 *        of equal or higher energy not overlapping with any reported
	 *        interaction so far; an interaction with energy E_INF is set, if
	 *        there is no better interaction left
	 */
	virtual
	void
	getNextBest( Interaction & curBest );


	/**
	 * Updates the mfe4leftEnd data container.
	 *
	 * Overwrite, if no update is to be done.
	 *
	 * @param i1 interaction start in seq1
	 * @param j1 interaction end in seq1
	 * @param i2 interaction start in seq2
	 * @param i2 interaction end in seq2
	 * @param curInteraction the interaction information to be used for update
	 */
	virtual
	void
	updateMfe4leftEnd(const size_t i1, const size_t j1
					, const size_t i2, const size_t j2
					, const Interaction & curInteraction );

	/**
	 * Calls for the stored mfe and suboptimal solutions traceBack(i)
	 * and pushes the according interactions to the output handler.
	 * For non-overlapping interaction enumeration, getNextBest() is called
	 * iteratively.
	 */
	virtual
	void
	reportOptima() override;

};

} // namespace

#endif /* INTARNA_PREDICTORMFE_H_ */
