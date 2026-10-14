#ifndef INTARNA_SEEDEDPARTITIONFUNCTION_H_
#define INTARNA_SEEDEDPARTITIONFUNCTION_H_
#include "IntaRNA/Interaction.h"
#include "IntaRNA/Matrix.h"
#include "IntaRNA/StackSeedDomain.h"
#include <functional>

namespace IntaRNA {
/** Nonnegative first-seed extension inside/outside computation. A path belongs
 * to its first completed admitted seed; seed-free left prefixes are joined to
 * unrestricted right extensions. Co-ending seeds use the shortest as owner.
 * Native energy ownership is supplied by the
 * adapter: initiation once, one transition per successive pair, one boundary
 * factor. All callbacks must be deterministic; zero denotes forbidden weight.
 */
class SeededPartitionFunction {
public:
	//! Local zero-based pair, second strand in reversed coordinates.
	using Pair = Interaction::BasePair;
	/** Active rectangle, inclusive spans, loop limits and structural restriction. */
	struct Domain { size_t n,m,span1,span2,loop1,loop2; bool noLP; };
	/** Energy adapter. Complete boundary callbacks are made once per positive
	 * seeded hybrid partition and receive the SAME coefficient used outside.
	 * No callback may alter the forward weights after observing them.
	 */
	struct Weights {
		std::function<bool(Pair)> valid;
		std::function<Z_type(Pair)> initiation;
		std::function<Z_type(Pair,Pair)> transition, boundary;
		std::function<void(Pair,Pair,Z_type,Z_type)> complete;
		/** Optional frozen heuristic selector, called once per canonical anchor.
		 * Arguments are its start and end; the result is a retained right endpoint.
		 * Keep paths with no left extension, no right extension, or that endpoint.
		 * Empty means all right endpoints (exact). Selection cannot depend on pair
		 * marking; out-of-domain choices contribute no two-sided extension.
		 */
		std::function<Pair(Pair,Pair)> selectRight;
	};
	/** Raw unnormalized result; pair masses remain empty in partition-only mode. */
	struct Result { Z_type z=0; Matrix<Z_type> mass; };
	/** Compute the union of seeded paths, optionally differentiating every pair.
	 * @param domain active dimensions and structural constraints
	 * @param seeds active admitted seed occurrences
	 * @param weights weight/admission contract (no integer energy reconstruction)
	 * @param pairs request the outside pass and raw local pair masses
	 * @return successful raw result, including exact zero for an empty ensemble
	 * @throws std::range_error on numerical failure; no partial result is returned
	 */
	static Result compute(const Domain & domain, const StackSeedDomain & seeds,
			const Weights & weights, bool pairs);
};
}
#endif
