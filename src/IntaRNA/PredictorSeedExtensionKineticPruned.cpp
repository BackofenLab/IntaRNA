#include "IntaRNA/PredictorSeedExtensionKineticPruned.h"

#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/InteractionEnergyVrna.h"

#include <algorithm>
#include <typeinfo>
#include <limits>

namespace IntaRNA {

PredictorSeedExtensionKineticPruned::PredictorSeedExtensionKineticPruned(
		const InteractionEnergy & model, OutputHandler & output,
		PredictionTracker * predTracker, SeedHandler * seedHandler, const char score)
	: PredictorSeedExtensionKinetic(model, output, predTracker, seedHandler, score)
{
	// Parameter-based bounds cannot describe arbitrary overrides of getE or
	// getE_interLeft. Do not silently apply them to derived energy models.
	const bool vrna = typeid(model) == typeid(InteractionEnergyVrna);
	const bool basePair = typeid(model) == typeid(InteractionEnergyBasePair);
	if ((!vrna && !basePair) || model.getMaxInternalLoopSize1() > 30
			|| model.getMaxInternalLoopSize2() > 30) {
		return;
	}
	rows = std::min(model.getMaxInternalLoopSize1(), model.size1() > 2 ? model.size1()-2 : 0)+1;
	columns = std::min(model.getMaxInternalLoopSize2(), model.size2() > 2 ? model.size2()-2 : 0)+1;
	// E_IntLoop takes a mutable pointer in supported ViennaRNA versions.
	// Work on a copy and retain no pointer to the model's parameter storage.
	vrna_param_t params;
	if (vrna) params = static_cast<const InteractionEnergyVrna &>(model).getVrnaParams();
	constexpr int reverse[] = {0,2,1,4,3,6,5};
	for (size_t side = 0; side < 2; ++side) {
		for (int root = 1; root <= 6; ++root) {
			auto & table = lowerBounds[side*6+root-1];
			table.assign(rows*columns, E_INF);
			for (size_t s1 = 0; s1 < rows; ++s1) {
				for (size_t s2 = 0; s2 < columns; ++s2) {
					E_type & bound = table[s1*columns+s2];
					if (basePair) {
						const std::int64_t two = 2*std::int64_t(model.getE_init());
						bound = static_cast<E_type>(std::clamp(two, std::int64_t(std::numeric_limits<E_type>::min()), std::int64_t(E_INF)));
						continue;
					}
					for (int close = 1; close <= 6; ++close) {
						E_type stack = E_INF;
						for (int outer = 1; outer <= 6; ++outer) {
							stack = std::min(stack, static_cast<E_type>(side == 0
									? params.stack[outer][reverse[close]] : params.stack[close][reverse[outer]]));
						}
						// Relax sequence consistency of the four mismatch bases.
						// This can only lower the local loop/stack minimum.
						for (int a = 1; a <= 4; ++a) for (int b = 1; b <= 4; ++b)
						for (int c = 1; c <= 4; ++c) for (int d = 1; d <= 4; ++d) {
							const E_type loop = E_IntLoop(s1, s2,
									side == 0 ? close : root, reverse[side == 0 ? root : close], a,b,c,d, &params);
							if (E_isNotINF(loop) && E_isNotINF(stack)) {
								const std::int64_t local = std::int64_t(loop)+stack;
								bound = std::min(bound, static_cast<E_type>(std::clamp(local,
										std::int64_t(std::numeric_limits<E_type>::min()), std::int64_t(E_INF))));
							}
						}
					}
				}
			}
			// Componentwise suffix minima make the bound nondecreasing along
			// either gap. With monotone ED this permits rectangle pruning.
			for (size_t i = rows; i-- > 0;) for (size_t j = columns; j-- > 0;) {
				E_type & value = table[i*columns+j];
				if (i+1 < rows) value = std::min(value, table[(i+1)*columns+j]);
				if (j+1 < columns) value = std::min(value, table[i*columns+j+1]);
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////

void
PredictorSeedExtensionKineticPruned::predict(const IndexRange & r1, const IndexRange & r2)
{
	edKnown = false;
	PredictorSeedExtensionKinetic::predict(r1, r2);
}

//////////////////////////////////////////////////////////////////////////

bool
PredictorSeedExtensionKineticPruned::prune(const Candidate & candidate, const Boundary & bounds) const
{
	if (rows == 0) return false;
	// Use global coordinates for the ED memo, including repeated predict()
	// calls with different range offsets but identical local boundaries.
	const Boundary global{bounds[0]+energy.getOffset1(), bounds[1]+energy.getOffset1(),
			bounds[2]+energy.getOffset2(), bounds[3]+energy.getOffset2()};
	if (!edKnown || edBounds != global) {
		edKnown = true; edBounds = global;
		currentED = std::int64_t(energy.getED1(bounds[0],bounds[1])) + energy.getED2(bounds[2],bounds[3]);
	}
	const int root = BP_pair[energy.getAccessibility1().getSequence().asCodes().at(global[candidate.left ? 0 : 1])]
			[energy.getAccessibility2().getSequence().asCodes().at(global[candidate.left ? 2 : 3])];
	if (root < 1 || root > 6 || candidate.s1 >= rows || candidate.s2 >= columns) return false;
	const E_type bound = lowerBounds[(candidate.left ? 0 : 6)+root-1][candidate.s1*columns+candidate.s2];
	const Boundary & b = candidate.bounds;
	const E_type ed1 = energy.getED1(b[0],b[1]), ed2 = energy.getED2(b[2],b[3]);
	return E_isINF(bound) || ed1 >= Accessibility::ED_UPPER_BOUND || ed2 >= Accessibility::ED_UPPER_BOUND
			|| std::int64_t(bound)+ed1+ed2-currentED >= 0;
}

} // namespace IntaRNA
