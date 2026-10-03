#ifndef INTARNA_PREDICTORSEEDEXTENSIONKINETIC_H_
#define INTARNA_PREDICTORSEEDEXTENSIONKINETIC_H_

#include "IntaRNA/PredictorMfe.h"
#include "IntaRNA/SeedHandlerIdxOffset.h"

#include <array>
#include <cstdint>
#include <map>

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
 * With noLP, every starting seed must already contain no lonely pair and a
 * nonstacking extension adds its closing pair and the following stack
 * atomically. Ties prefer left, smaller s1+s2, then smaller s1. Every valid
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

private:
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
		E_type hybrid = E_INF;
		E_type total = E_INF;
		std::int64_t delta = 0;
	};

	//! Owned seed handler with offsets matching this->energy.
	SeedHandlerIdxOffset seedHandler;
	//! The selected deterministic move-ranking formula.
	const char score;
	//! Paths retained independently of the number of requested reports.
	InteractionCache interactions;
	//! Valid starting seeds with recomputed energies, keyed by original indices.
	std::map<Interaction::BasePair, Interaction::Seed> validSeeds;

	/** @return local, inclusive boundaries of a nonempty interaction */
	Boundary getBoundary(const Interaction & interaction) const;

	/**
	 * Checks the complete seed path, including noLP and loop GU constraints.
	 * @param interaction seed path, in original sequence coordinates
	 * @param hybrid receives the traced path's loop energies plus initiation
	 * @return whether every pair and adjacent loop is feasible
	 */
	bool isValidSeed(const Interaction & interaction, E_type & hybrid) const;

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
	 * Checks one geometrically bounded move and evaluates its complete energy.
	 * @param candidate move indices/side/gaps; energies are filled on success
	 * @param bounds current interaction boundaries
	 * @param hybrid current hybridization energy, including initiation
	 * @param total current full interaction energy
	 * @return whether the entire move is feasible and strictly downhill
	 */
	bool evaluate(Candidate & candidate, const Boundary & bounds,
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
