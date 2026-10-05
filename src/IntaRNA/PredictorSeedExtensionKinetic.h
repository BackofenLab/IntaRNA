#ifndef INTARNA_PREDICTORSEEDEXTENSIONKINETIC_H_
#define INTARNA_PREDICTORSEEDEXTENSIONKINETIC_H_

#include "IntaRNA/PredictorMfe.h"
#include "IntaRNA/SeedHandlerIdxOffset.h"

#include <array>
#include <cstdint>
#include <map>
#include <vector>

namespace IntaRNA {

/**
 * Deterministic, strictly downhill extension of every feasible seed.
 *
 * At each step both ends compete using the complete interaction-energy
 * difference, including accessibility, terminal penalties and both weighted
 * dangling ends. Score A uses this difference directly; B divides by
 * 1+s1+s2; C (C1) divides by 1+2*max(s1,s2). These are heuristic move rankings,
 * not physical rates or a calibrated folding-time model. Only negative
 * energy differences are accepted, including for the normalized scores.
 *
 * Extensions always add one stacked pair or two stacked pairs, including
 * across a loop. Seeds and their energies are trusted as supplied by the seed
 * handler, even when an explicit seed contains lonely pairs. Ties prefer left,
 * smaller s1+s2, smaller s1, then the single-pair move. Every valid
 * visited prefix is eligible for normal MFE/suboptimal reporting; traceback
 * reproduces the actual greedy path. Equilibrium partition-function output
 * is unsupported.
 *
 * Candidate enumeration uses the active energy model and separate loop/span
 * limits for both RNAs. It deliberately has no loop-only energy pruning:
 * such bounds omit favorable mandatory stacks and changes to the opposite
 * dangling end, and accessibility differences need not be monotone.
 */
class PredictorSeedExtensionKinetic : public PredictorMfe {
public:

	/**
	 * Constructs a predictor, taking ownership of tracker and seed handler.
	 * @param energy energy model, which must outlive this predictor
	 * @param output output handler, which must outlive this predictor
	 * @param predTracker owned tracker, or NULL
	 * @param seedHandler owned, non-NULL seed handler
	 * @param score move ranking: A, B, or C (the C1 formula)
	 * @throws std::invalid_argument for a NULL seed handler, unknown score or
	 *         output requiring an equilibrium partition function
	 */
	PredictorSeedExtensionKinetic(const InteractionEnergy & energy,
			OutputHandler & output, PredictionTracker * predTracker,
			SeedHandler * seedHandler, const char score = 'A');

	/** Frees the owned seed handler and prediction tracker. */
	virtual ~PredictorSeedExtensionKinetic();

	/**
	 * Extends all feasible seeds within inclusive, zero-based sequence ranges.
	 * Sequence 2 uses the reversed indexing of the energy model. Repeated
	 * calls reset all trajectories, cached output and index offsets.
	 * @param r1 permitted range in sequence 1
	 * @param r2 permitted range in reversed sequence 2
	 * @throws std::invalid_argument for an invalid or empty input range
	 */
	void predict(const IndexRange & r1 = IndexRange(0, RnaSequence::lastPos),
			const IndexRange & r2 = IndexRange(0, RnaSequence::lastPos)) override;

protected:
	/**
	 * Restores the exact stored greedy path and annotates its contained seeds.
	 * @param interaction interaction boundaries to expand
	 */
	void traceBack(Interaction & interaction) override;

	/**
	 * Finds the best cached prefix disjoint from already reported intervals.
	 * @param interaction current report, replaced by the next report or E_INF
	 */
	void getNextBest(Interaction & interaction) override;

	//! Inclusive boundaries (i1,j1,i2,j2), using local energy indices.
	using Boundary = std::array<size_t, 4>;
	//! Best actual path for each visited, reportable set of boundaries.
	using InteractionCache = std::map<Boundary, Interaction>;

	/** A feasible move; delta is widened before subtraction. */
	struct Candidate {
		Boundary bounds = {};
		size_t close1 = 0;
		size_t close2 = 0;
		size_t s1 = 0;
		size_t s2 = 0;
		bool left = true;
		bool macro = false;
		bool topologyKnown = false;
		bool topologyAllowed = false;
		bool active = false;
		bool localKnown = false;
		E_type local = E_INF;
		E_type hybrid = E_INF;
		E_type total = E_INF;
		std::int64_t delta = 0;
	};

	/**
	 * Optional filter before pairing and loop-energy lookups. The default
	 * enumerates every move; subclasses may implement heuristic pruning.
	 * @param candidate geometrically valid extension
	 * @param bounds current boundaries
	 * @return whether to skip this two-pair move and all moves with both gaps
	 *         at least as large, in the current state (requires monotone ED)
	 */
	virtual bool prune(const Candidate & candidate, const Boundary & bounds) const;

private:
	/** Geometry, shared pair checks and local energies for one unchanged end. */
	struct SideCandidates {
		std::vector<Candidate> moves;
		std::vector<signed char> complementary;
		size_t columns = 0;
	};

	//! Owned seed handler with offsets matching this->energy.
	SeedHandlerIdxOffset seedHandler;
	//! The selected deterministic move-ranking formula.
	const char score;
	//! Paths retained independently of the number of requested reports.
	InteractionCache interactions;
	/** @return local, inclusive boundaries of a nonempty interaction */
	Boundary getBoundary(const Interaction & interaction) const;

	/**
	 * Runs a complete greedy trajectory and retains its reportable prefixes.
	 * @param interaction initial seed, modified to its final state
	 * @param hybrid initial seed hybridization energy including initiation
	 * @param last1 last permitted local index in sequence 1
	 * @param last2 last permitted local index in reversed sequence 2
	 */
	void extendSeed(Interaction & interaction, E_type hybrid,
			size_t last1, size_t last2);

	/**
	 * Builds the rectangular gap table plus the single-stack move for one end.
	 * @param side table to replace
	 * @param bounds current boundaries
	 * @param left whether this is the left end
	 * @param last1 last permitted local index in sequence 1
	 * @param last2 last permitted local index in reversed sequence 2
	 */
	void buildCandidates(SideCandidates & side, const Boundary & bounds,
			bool left, size_t last1, size_t last2) const;

	/**
	 * Resolves shared pair checks before energies, then refreshes complete
	 * energies and selects the best downhill candidate during that traversal.
	 * @param side candidate and pairing cache for one end
	 * @param bounds current boundaries (the opposite end may have changed)
	 * @param hybrid current hybridization energy including initiation
	 * @param total current complete interaction energy
	 * @return best move within side, or NULL if none is downhill
	 */
	const Candidate * updateCandidates(SideCandidates & side, const Boundary & bounds,
			E_type hybrid, E_type total) const;

	/** @return whether a move wins by exact score and deterministic ties */
	bool isBetter(const Candidate & candidate, const Candidate & best) const;

	/**
	 * Reduces a valid prefix by boundaries, total energy and full-path ties.
	 * @param interaction actual path to consider for reporting
	 */
	void retain(const Interaction & interaction);
};

} // namespace IntaRNA

#endif /* INTARNA_PREDICTORSEEDEXTENSIONKINETIC_H_ */
