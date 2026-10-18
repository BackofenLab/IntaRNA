
#ifndef INTARNA_PREDICTORMFEENS_H_
#define INTARNA_PREDICTORMFEENS_H_


#include "IntaRNA/PredictorMfe.h"

#include "IntaRNA/IndexRangeList.h"

#include <list>
#include <unordered_map>
#include <utility>
#include <vector>
#include <boost/functional/hash.hpp>

namespace IntaRNA {

/**
 * Generic Predictor interface for ensemble-based MFE interaction computation
 *
 * @author Martin Raden
 *
 */
class PredictorMfeEns : public PredictorMfe {


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
	PredictorMfeEns( const InteractionEnergy & energy
				, OutputHandler & output
				, PredictionTracker * predTracker );

	virtual ~PredictorMfeEns();

protected:

	//! data container to encode a site with respective partition function
	struct ZPartition {
		size_t i1;
		size_t j1;
		size_t i2;
		size_t j2;
		Z_type partZ;
	};

	//! access to the interaction energy handler of the super class
	using PredictorMfe::energy;

	//! access to the output handler of the super class
	using PredictorMfe::output;

	//! access to the prediction tracker of the super class
	using PredictorMfe::predTracker;

	//! map storing the partition of Zall for all considered interaction sites
	std::unordered_map<Interaction::Boundary, Z_type, Interaction::Boundary::Hash> Z_partition;

	/** Checked exact-path update with the forward/outside boundary coefficient.
	 * Preserves virtual updateZ dispatch. Observing/delegating overrides work;
	 * reweighting must override the objective coefficient instead of partFunct.
	 */
	void updateExactCompleteZ(size_t i1, size_t j1, size_t i2, size_t j2,
			Z_type hybrid, Z_type coefficient);

	//! scoped coefficient and hybrid expected by the checked complete update
	const std::pair<Z_type,Z_type> * exactContribution = nullptr;

	//! whether updateZ() currently receives one complete boundary partition
	bool updateZisComplete;


	//! a bucket supplies j1, so its key needs only the other three coordinates
	struct FrontierBoundary {
		size_t i1, i2, j2;
		bool operator==(const FrontierBoundary &) const = default;
		struct Hash {
			size_t operator()(const FrontierBoundary & boundary) const;
		};
	};

	//! partial sums grouped by target right boundary; used only without trackers
	std::vector<std::unordered_map<FrontierBoundary,Z_type,FrontierBoundary::Hash>> Z_frontier;
	//! first right-boundary bucket that has not yet been finalized
	size_t Z_frontierNext = 0;


	/**
	 * Initializes the hybridization partition functions.
	 * Will be called by predict().
	 */
	virtual
	void
	initZ();

	/**
	 * Updates the local hybridization partition functions as well as Zall.
	 *
	 * Note: avoid other Zall updates via updateOptima(.., incrementZall=false)!
	 *
	 * Note: if called multiple time for the same boundaries then the
	 * reported partition functions have to represent disjoint interaction sets!
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1
	 * @param partFunct the partition function of the interaction
	 * @param isHybridZ whether or not the given Z is only covering
	 *        hybridization energy terms (init+loops) or the total
	 *        interaction energy
	 */
	virtual
	void
	updateZ( const size_t i1, const size_t j1
				, const size_t i2, const size_t j2
				, const Z_type partFunct
				, const bool isHybridZ );

	/**
	 * Consumes a complete partition for one interaction boundary. Exact 2D
	 * prediction finalizes each boundary once and therefore does not need to
	 * retain it until reportOptima(). The virtual updateZ() hook is called so
	 * subclasses can observe or alter the update; overrides that delegate to the
	 * base implementation inherit streaming. Tracker-enabled predictions retain
	 * the partition so their established callback timing and ordering stay intact.
	 */
	void
	updateCompleteZ( const size_t i1, const size_t j1
				, const size_t i2, const size_t j2
				, const Z_type partFunct
				, const bool isHybridZ );

	/** Enable partial-partition buckets by target right boundary. Tracker runs
	 * keep the original map to preserve callback ordering. Call advanceZFrontier
	 * only when no remaining seed can start below the supplied target index.
	 * @param size1 number of target positions in the current prediction region
	 */
	void beginZFrontier(size_t size1);

	/** Consume completed boundary sums without adding their mass to Zall again.
	 * @param minimumSeedStart smallest target start of any remaining seed
	 */
	void advanceZFrontier(size_t minimumSeedStart);

	/**
	 * Calls for the stored Z_partition information updateOptima() before
	 * calling reportOptima() from its super class.
	 */
	virtual
	void
	reportOptima();

	/**
	 * Reports interaction boundaries only (no base pair tracking)
	 * @param interaction the interaction to be traced
	 */
	virtual
	void
	traceBack( Interaction & interaction );

private:

	/**
	 * Applies output filters, converts the partition representation and adds
	 * the accepted contribution to Zall. Returns the partition without ED.
	 */
	bool
	addPartitionContribution( const size_t i1, const size_t j1
				, const size_t i2, const size_t j2
				, const Z_type partFunct
				, const bool isHybridZ
				, Z_type & partZ_noED );

	/**
	 * Calls updateOptima() for each buffered entry of Z_partition. This hook is
	 * still called by reportOptima(), but complete tracker-free partitions were
	 * already consumed by updateZ() and are therefore not present in the map.
	 */
	virtual
	void
	updateOptimaUsingZ();

};

inline
size_t
PredictorMfeEns::FrontierBoundary::Hash::operator()(const FrontierBoundary & boundary) const
{
	size_t key=0;
	boost::hash_combine(key,boundary.i1);
	boost::hash_combine(key,boundary.i2);
	boost::hash_combine(key,boundary.j2);
	return key;
}

} // namespace

#endif /* INTARNA_PREDICTORMFEENS_H_ */
