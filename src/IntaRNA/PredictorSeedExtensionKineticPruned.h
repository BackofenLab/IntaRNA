#ifndef INTARNA_PREDICTORSEEDEXTENSIONKINETICPRUNED_H_
#define INTARNA_PREDICTORSEEDEXTENSIONKINETICPRUNED_H_

#include "IntaRNA/PredictorSeedExtensionKinetic.h"

namespace IntaRNA {

/**
 * Experimental kinetic extension with root-pair-specific loop/stack bounds.
 *
 * Before complementarity checks, a suffix-minimum table and the current ED
 * increment can reject a rectangle of larger two-pair moves. Local estimates
 * use the active ViennaRNA parameters (or the base-pair energy). Pruning
 * assumes monotone ED and ignores changes in terminal/dangling contributions:
 * it is intentionally heuristic and can change the path compared with K.
 * Every surviving move is still checked with its complete energy. Single
 * stacks are always evaluated. Unknown energy subclasses and loop limits
 * above 30 fall back to exhaustive K enumeration.
 */
class PredictorSeedExtensionKineticPruned : public PredictorSeedExtensionKinetic {
public:
	/**
	 * Constructs the experimental predictor and precomputes local bounds.
	 * @param energy energy model that must outlive this predictor
	 * @param output output handler that must outlive this predictor
	 * @param predTracker owned prediction tracker, or NULL
	 * @param seedHandler owned, non-NULL seed handler
	 * @param score move ranking A, B or C
	 */
	PredictorSeedExtensionKineticPruned(const InteractionEnergy & energy,
			OutputHandler & output, PredictionTracker * predTracker,
			SeedHandler * seedHandler, char score = 'A');

	/**
	 * Resets the ED memo and runs kinetic prediction within the given ranges.
	 * @param r1 permitted inclusive range in sequence 1
	 * @param r2 permitted inclusive range in reversed sequence 2
	 */
	void predict(const IndexRange & r1 = IndexRange(0, RnaSequence::lastPos),
			const IndexRange & r2 = IndexRange(0, RnaSequence::lastPos)) override;

protected:
	/**
	 * Tests the local suffix bound plus ED increase, omitting terminal and
	 * dangling changes. A true result rejects all componentwise larger gaps.
	 * @param candidate geometrically valid two-pair extension
	 * @param bounds current interaction boundaries
	 * @return whether to prune the candidate and its larger-gap rectangle
	 */
	bool prune(const Candidate & candidate, const Boundary & bounds) const override;

private:
	//! Left/right tables for the six oriented canonical base-pair types.
	std::array<std::vector<E_type>, 12> lowerBounds;
	//! Rectangular gap-table dimensions; zero disables the heuristic.
	size_t rows = 0, columns = 0;
	//! Current-state ED is shared by all candidates at both ends.
	mutable Boundary edBounds = {};
	mutable bool edKnown = false;
	mutable std::int64_t currentED = 0;
};

} // namespace IntaRNA

#endif /* INTARNA_PREDICTORSEEDEXTENSIONKINETICPRUNED_H_ */
