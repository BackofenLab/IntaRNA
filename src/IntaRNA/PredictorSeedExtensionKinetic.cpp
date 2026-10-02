#include "IntaRNA/PredictorSeedExtensionKinetic.h"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <utility>

#include <boost/multiprecision/cpp_int.hpp>

namespace IntaRNA {
namespace {

SeedHandler * checkedSeedHandler(SeedHandler * handler)
{
	if (handler == NULL) {
		throw std::invalid_argument("PredictorSeedExtensionKinetic requires a seed handler");
	}
	return handler;
}

// Do not add infinity sentinels or overflow the internal integer energy type.
E_type addEnergy(const E_type first, const E_type second)
{
	if (E_isINF(first) || E_isINF(second)) {
		return E_INF;
	}
	const std::int64_t sum = std::int64_t(first) + std::int64_t(second);
	return sum >= E_INF || sum < std::numeric_limits<E_type>::min()
			? E_INF : static_cast<E_type>(sum);
}

} // namespace

//////////////////////////////////////////////////////////////////////////

PredictorSeedExtensionKinetic::PredictorSeedExtensionKinetic(
		const InteractionEnergy & energy, OutputHandler & output,
		PredictionTracker * predTracker, SeedHandler * seedHandlerInstance,
		const char score)
	: PredictorMfe(energy, output, predTracker)
	, seedHandler(checkedSeedHandler(seedHandlerInstance))
	, score(score)
	, interactions()
	, validSeeds()
{
	if (score != 'A' && score != 'B' && score != 'C') {
		throw std::invalid_argument("PredictorSeedExtensionKinetic score must be A, B or C");
	}
	if (output.getOutputConstraint().needZall) {
		throw std::invalid_argument("PredictorSeedExtensionKinetic does not compute an equilibrium partition function");
	}
}

//////////////////////////////////////////////////////////////////////////

PredictorSeedExtensionKinetic::~PredictorSeedExtensionKinetic()
{
}

//////////////////////////////////////////////////////////////////////////

void
PredictorSeedExtensionKinetic::predict(const IndexRange & r1, const IndexRange & r2)
{
	const size_t size1 = energy.getAccessibility1().getSequence().size();
	const size_t size2 = energy.getAccessibility2().getSequence().size();
	if (!r1.isAscending() || !r2.isAscending()
			|| r1.from >= size1 || r2.from >= size2) {
		throw std::invalid_argument("PredictorSeedExtensionKinetic::predict(): invalid sequence range");
	}

	energy.setOffset1(r1.from);
	energy.setOffset2(r2.from);
	seedHandler.setOffset1(r1.from);
	seedHandler.setOffset2(r2.from);
	const size_t last1 = std::min(r1.to, size1 - 1) - r1.from;
	const size_t last2 = std::min(r2.to, size2 - 1) - r2.from;
	interactions.clear();
	validSeeds.clear();
	initOptima();

	if (seedHandler.fillSeed(0, last1, 0, last2) != 0) {
		size_t i1 = RnaSequence::lastPos, i2 = RnaSequence::lastPos;
		while (seedHandler.updateToNextSeed(i1, i2, 0, last1, 0, last2)) {
			const size_t length1 = seedHandler.getSeedLength1(i1, i2);
			const size_t length2 = seedHandler.getSeedLength2(i1, i2);
			if (length1 == 0 || length2 == 0
					|| length1 - 1 > last1 - i1 || length2 - 1 > last2 - i2
					|| length1 > energy.getAccessibility1().getMaxLength()
					|| length2 > energy.getAccessibility2().getMaxLength()) {
				continue;
			}
			const size_t j1 = i1 + length1 - 1;
			const size_t j2 = i2 + length2 - 1;
			if (E_isINF(seedHandler.getSeedE(i1, i2))) {
				continue;
			}

			Interaction interaction(energy.getAccessibility1().getSequence(),
					energy.getAccessibility2().getAccessibilityOrigin().getSequence());
			interaction.basePairs.push_back(energy.getBasePair(i1, i2));
			seedHandler.traceBackSeed(interaction, i1, i2);
			if (i1 != j1 || i2 != j2) {
				interaction.basePairs.push_back(energy.getBasePair(j1, j2));
			}
			interaction.sort();
			E_type hybrid = E_INF;
			if (!isValidSeed(interaction, hybrid)
					|| getBoundary(interaction) != Boundary{i1, j1, i2, j2}) {
				continue;
			}
			interaction.energy = energy.getE(i1, j1, i2, j2, hybrid);
			if (E_isINF(interaction.energy)) {
				continue;
			}
			interaction.setSeedRange(interaction.basePairs.front(),
					interaction.basePairs.back(), interaction.energy);
			validSeeds.emplace(interaction.basePairs.front(),
					Interaction::Seed(interaction.basePairs.front(),
							interaction.basePairs.back(), interaction.energy));
			extendSeed(interaction, hybrid, last1, last2);
		}
	}

	// Reduce identical boundaries before updating optima: an earlier, inferior
	// path must never be traced back using a later replacement's base pairs.
	for (const auto & entry : interactions) {
		const Boundary & b = entry.first;
		updateOptima(b[0], b[1], b[2], b[3], entry.second.energy, false, false);
	}
	// The generic reporter assumes a nonempty optimum list. Zero reports can
	// still be useful for prediction trackers and must not dereference it.
	if (output.getOutputConstraint().reportMax != 0) {
		reportOptima();
	}
}

//////////////////////////////////////////////////////////////////////////

PredictorSeedExtensionKinetic::Boundary
PredictorSeedExtensionKinetic::getBoundary(const Interaction & interaction) const
{
	return Boundary{energy.getIndex1(interaction.basePairs.front()),
			energy.getIndex1(interaction.basePairs.back()),
			energy.getIndex2(interaction.basePairs.front()),
			energy.getIndex2(interaction.basePairs.back())};
}

//////////////////////////////////////////////////////////////////////////

bool
PredictorSeedExtensionKinetic::isValidSeed(const Interaction & interaction, E_type & hybrid) const
{
	if (interaction.basePairs.empty() || !interaction.isValid()) {
		return false;
	}
	const auto & pairs = interaction.basePairs;
	const auto & constraint = output.getOutputConstraint();
	// Reconstruct rather than trust a cached explicit-seed energy: the
	// trajectory's energy must correspond to precisely the traced structure.
	hybrid = energy.getE_init();
	if (E_isINF(hybrid)) {
		return false;
	}
	for (size_t p = 0; p < pairs.size(); ++p) {
		const size_t i1 = energy.getIndex1(pairs[p]);
		const size_t i2 = energy.getIndex2(pairs[p]);
		if (i1 >= energy.size1() || i2 >= energy.size2()
				|| !energy.areComplementary(i1, i2)) {
			return false;
		}
		const bool stackedLeft = p > 0
				&& pairs[p].first - pairs[p-1].first == 1
				&& pairs[p-1].second - pairs[p].second == 1;
		const bool stackedRight = p + 1 < pairs.size()
				&& pairs[p+1].first - pairs[p].first == 1
				&& pairs[p].second - pairs[p+1].second == 1;
		if (constraint.noLP && !stackedLeft && !stackedRight) {
			return false;
		}
		if (p == 0) {
			continue;
		}
		const size_t previous1 = energy.getIndex1(pairs[p-1]);
		const size_t previous2 = energy.getIndex2(pairs[p-1]);
		if (!stackedLeft && constraint.noGUend
				&& (energy.isGU(previous1, previous2) || energy.isGU(i1, i2))) {
			return false;
		}
		hybrid = addEnergy(hybrid, energy.getE_interLeft(previous1, i1, previous2, i2));
		if (E_isINF(hybrid)) {
			return false;
		}
	}
	return true;
}

//////////////////////////////////////////////////////////////////////////

void
PredictorSeedExtensionKinetic::extendSeed(Interaction & interaction,
		E_type hybrid, const size_t last1, const size_t last2)
{
	const size_t maxLength1 = energy.getAccessibility1().getMaxLength();
	const size_t maxLength2 = energy.getAccessibility2().getMaxLength();
	while (true) {
		retain(interaction);
		const Boundary bounds = getBoundary(interaction);
		const size_t remaining1 = maxLength1 - (bounds[1] - bounds[0] + 1);
		const size_t remaining2 = maxLength2 - (bounds[3] - bounds[2] + 1);
		bool found = false;
		Candidate best;
		for (unsigned int side = 0; side < 2; ++side) {
			const bool left = side == 0;
			const size_t space1 = std::min(remaining1, left ? bounds[0] : last1 - bounds[1]);
			const size_t space2 = std::min(remaining2, left ? bounds[2] : last2 - bounds[3]);
			if (space1 == 0 || space2 == 0) {
				continue;
			}
			// A GU boundary can still grow by stacking, but cannot start a
			// nonstacking loop when either active constraint forbids it.
			const bool stackOnly = (output.getOutputConstraint().noGUend
					|| !energy.isInternalLoopGUallowed())
					&& energy.isGU(bounds[left ? 0 : 1], bounds[left ? 2 : 3]);
			const size_t maxGap1 = stackOnly ? 0
					: std::min(energy.getMaxInternalLoopSize1(), space1 - 1);
			const size_t maxGap2 = stackOnly ? 0
					: std::min(energy.getMaxInternalLoopSize2(), space2 - 1);
			for (size_t s1 = 0; s1 <= maxGap1; ++s1) {
				for (size_t s2 = 0; s2 <= maxGap2; ++s2) {
					Candidate candidate;
					candidate.left = left;
					candidate.s1 = s1;
					candidate.s2 = s2;
					candidate.macro = output.getOutputConstraint().noLP && (s1 != 0 || s2 != 0);
					const size_t addedPairs = candidate.macro ? 2 : 1;
					if (space1 < addedPairs || space2 < addedPairs
							|| s1 > space1 - addedPairs || s2 > space2 - addedPairs) {
						continue;
					}
					candidate.bounds = bounds;
					if (left) {
						candidate.bounds[0] -= s1 + addedPairs;
						candidate.bounds[2] -= s2 + addedPairs;
						candidate.close1 = bounds[0] - s1 - 1;
						candidate.close2 = bounds[2] - s2 - 1;
					} else {
						candidate.bounds[1] += s1 + addedPairs;
						candidate.bounds[3] += s2 + addedPairs;
						candidate.close1 = bounds[1] + s1 + 1;
						candidate.close2 = bounds[3] + s2 + 1;
					}
					if (evaluate(candidate, bounds, hybrid, interaction.energy)
							&& (!found || isBetter(candidate, best))) {
						best = candidate;
						found = true;
					}
				}
			}
		}
		if (!found) {
			break;
		}
		const Interaction::BasePair close = energy.getBasePair(best.close1, best.close2);
		if (best.left) {
			interaction.basePairs.insert(interaction.basePairs.begin(), close);
			if (best.macro) {
				interaction.basePairs.insert(interaction.basePairs.begin(),
						energy.getBasePair(best.bounds[0], best.bounds[2]));
			}
		} else {
			interaction.basePairs.push_back(close);
			if (best.macro) {
				interaction.basePairs.push_back(energy.getBasePair(best.bounds[1], best.bounds[3]));
			}
		}
		hybrid = best.hybrid;
		interaction.energy = best.total;
	}
}

//////////////////////////////////////////////////////////////////////////

bool
PredictorSeedExtensionKinetic::evaluate(Candidate & candidate,
		const Boundary & bounds, const E_type hybrid, const E_type total) const
{
	const size_t old1 = bounds[candidate.left ? 0 : 1];
	const size_t old2 = bounds[candidate.left ? 2 : 3];
	const size_t outer1 = candidate.bounds[candidate.left ? 0 : 1];
	const size_t outer2 = candidate.bounds[candidate.left ? 2 : 3];
	if (!energy.areComplementary(candidate.close1, candidate.close2)
			|| (candidate.macro && !energy.areComplementary(outer1, outer2))) {
		return false;
	}
	if (output.getOutputConstraint().noGUend && (candidate.s1 != 0 || candidate.s2 != 0)
			&& (energy.isGU(old1, old2) || energy.isGU(candidate.close1, candidate.close2))) {
		return false;
	}
	const E_type loop = candidate.left
			? energy.getE_interLeft(candidate.close1, old1, candidate.close2, old2)
			: energy.getE_interLeft(old1, candidate.close1, old2, candidate.close2);
	candidate.hybrid = addEnergy(hybrid, loop);
	if (candidate.macro && E_isNotINF(candidate.hybrid)) {
		const E_type stack = candidate.left
				? energy.getE_interLeft(outer1, candidate.close1, outer2, candidate.close2)
				: energy.getE_interLeft(candidate.close1, outer1, candidate.close2, outer2);
		candidate.hybrid = addEnergy(candidate.hybrid, stack);
	}
	if (E_isINF(candidate.hybrid)) {
		return false;
	}
	const Boundary & b = candidate.bounds;
	candidate.total = energy.getE(b[0], b[1], b[2], b[3], candidate.hybrid);
	if (E_isINF(candidate.total)) {
		return false;
	}
	candidate.delta = std::int64_t(candidate.total) - std::int64_t(total);
	return candidate.delta < 0;
}

//////////////////////////////////////////////////////////////////////////

bool
PredictorSeedExtensionKinetic::isBetter(const Candidate & candidate, const Candidate & best) const
{
	// A 32-bit energy difference times a 65-bit gap denominator fits in
	// 128 bits, including for public-API loop limits beyond the CLI limits.
	using Wide = boost::multiprecision::int128_t;
	const auto denominator = [this](const Candidate & c) -> Wide {
		if (score == 'B') {
			return Wide(1) + Wide(c.s1) + Wide(c.s2);
		}
		if (score == 'C') {
			return Wide(1) + 2 * Wide(std::max(c.s1, c.s2));
		}
		return Wide(1);
	};
	const Wide lhs = Wide(candidate.delta) * denominator(best);
	const Wide rhs = Wide(best.delta) * denominator(candidate);
	if (lhs != rhs) {
		return lhs < rhs;
	}
	if (candidate.left != best.left) {
		return candidate.left;
	}
	const Wide size = Wide(candidate.s1) + Wide(candidate.s2);
	const Wide bestSize = Wide(best.s1) + Wide(best.s2);
	return size != bestSize ? size < bestSize : candidate.s1 < best.s1;
}

//////////////////////////////////////////////////////////////////////////

void
PredictorSeedExtensionKinetic::retain(const Interaction & interaction)
{
	const Boundary b = getBoundary(interaction);
	const auto & constraint = output.getOutputConstraint();
	if (interaction.energy >= E_MAX
			|| (constraint.noGUend && (energy.isGU(b[0], b[2]) || energy.isGU(b[1], b[3])))
			|| energy.getED1(b[0], b[1]) > constraint.maxED
			|| energy.getED2(b[2], b[3]) > constraint.maxED) {
		return;
	}
	auto existing = interactions.find(b);
	if (existing == interactions.end()) {
		interactions.emplace(b, interaction);
	} else if (interaction.energy < existing->second.energy
			|| (interaction.energy == existing->second.energy
					&& interaction.basePairs < existing->second.basePairs)) {
		existing->second = interaction;
	}
}

//////////////////////////////////////////////////////////////////////////

void
PredictorSeedExtensionKinetic::traceBack(Interaction & interaction)
{
	if (interaction.basePairs.empty()) {
		return;
	}
	const auto path = interactions.find(getBoundary(interaction));
	if (path == interactions.end() || path->second.energy != interaction.energy) {
		throw std::runtime_error("PredictorSeedExtensionKinetic::traceBack(): no matching greedy path");
	}
	interaction = path->second;
	seedHandler.addSeeds(interaction);
	// The generic annotator can recognize explicit seeds that were rejected
	// as starting states (e.g. a lonely seed end stacked only by extension).
	// Retain only validated starts and use their reconstructed energies.
	if (interaction.seed != NULL) {
		Interaction::SeedSet validAnnotations;
		for (const Interaction::Seed & seed : *interaction.seed) {
			const auto valid = validSeeds.find(seed.bp_i);
			if (valid != validSeeds.end() && valid->second.bp_j == seed.bp_j) {
				validAnnotations.insert(valid->second);
			}
		}
		*interaction.seed = std::move(validAnnotations);
	}
}

//////////////////////////////////////////////////////////////////////////

void
PredictorSeedExtensionKinetic::getNextBest(Interaction & interaction)
{
	const E_type previousEnergy = interaction.energy;
	const Interaction * best = NULL;
	for (const auto & entry : interactions) {
		const Boundary & b = entry.first;
		const Interaction & candidate = entry.second;
		if (candidate.energy < previousEnergy
				|| reportedInteractions.first.overlaps(IndexRange(b[0], b[1]))
				|| reportedInteractions.second.overlaps(IndexRange(b[2], b[3]))) {
			continue;
		}
		if (best == NULL || candidate < *best) {
			best = &candidate;
		}
	}
	if (best == NULL) {
		interaction.clear();
		interaction.energy = E_INF;
		return;
	}
	interaction = *best;
	const Interaction::BasePair right = interaction.basePairs.back();
	interaction.basePairs.resize(interaction.basePairs.size() == 1 ? 1 : 2);
	interaction.basePairs.back() = right;
	INTARNA_CLEANUP(interaction.seed);
}

} // namespace IntaRNA
