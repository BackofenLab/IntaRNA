#ifndef INTARNA_OUTPUTHANDLERSEEDCANDIDATES_H_
#define INTARNA_OUTPUTHANDLERSEEDCANDIDATES_H_

#include "IntaRNA/OutputHandler.h"
#include <unordered_map>
#include <vector>

namespace IntaRNA {

/** Collect each seed's best distinct sites across prediction calls/windows.
 * Seed keys and interactions use original sequence coordinates. Final selection
 * deduplicates sites, retaining their lowest eligible energy and its structure,
 * then applies energy, overlap and count restrictions without refilling quotas.
 * Partition statistics are accumulated independently through OutputHandler.
 */
class OutputHandlerSeedCandidates : public OutputHandler {
public:
	/** @param constraint final restrictions, including a positive per-seed quota */
	explicit OutputHandlerSeedCandidates(const OutputConstraint & constraint);
	/** Ordinary ranked reports have no seed identity and are ignored.
	 * @param interaction ordinary report from the statistics-only prediction
	 */
	void add(const Interaction & interaction) override;
	/** Retain a candidate before any global output pruning.
	 * @param seed generating seed start in original coordinates
	 * @param interaction candidate, including its qualifying conformation/energy
	 */
	void addSeed(const Interaction::BasePair & seed, const Interaction & interaction);
	/** Publish the selected union, leaving partition forwarding to the caller.
	 * @param destination output receiving at most reportMax distinct sites
	 */
	void report(OutputHandler & destination) const;
	/** Reset quotas for an independent region combination; retain statistics. */
	void clear();
private:
	//! Only the best K distinct sites per seed are stored, with deep-copy ownership.
	std::unordered_map<Interaction::BasePair, std::vector<Interaction>, Interaction::BasePair::Hash> candidates;
};

} // namespace IntaRNA
#endif
