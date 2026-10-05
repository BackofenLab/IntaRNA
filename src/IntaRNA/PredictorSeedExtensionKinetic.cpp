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
			const E_type hybrid = addEnergy(seedHandler.getSeedE(i1, i2), energy.getE_init());
			interaction.energy = energy.getE(i1, j1, i2, j2, hybrid);
			if (E_isINF(interaction.energy)) {
				continue;
			}
			interaction.setSeedRange(interaction.basePairs.front(),
					interaction.basePairs.back(), interaction.energy);
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
PredictorSeedExtensionKinetic::prune(const Candidate &, const Boundary &) const
{
	return false;
}

//////////////////////////////////////////////////////////////////////////

void
PredictorSeedExtensionKinetic::extendSeed(Interaction & interaction,
		E_type hybrid, const size_t last1, const size_t last2)
{
	std::array<SideCandidates, 2> sides;
	Boundary bounds = getBoundary(interaction);
	buildCandidates(sides[0], bounds, true, last1, last2);
	buildCandidates(sides[1], bounds, false, last1, last2);
	while (true) {
		retain(interaction);
		const Candidate * left = updateCandidates(sides[0], bounds, hybrid, interaction.energy);
		const Candidate * right = updateCandidates(sides[1], bounds, hybrid, interaction.energy);
		if (left == NULL && right == NULL) {
			break;
		}
		const Candidate best = left != NULL && (right == NULL || isBetter(*left, *right)) ? *left : *right;
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
		bounds = best.bounds;
		// The opposite end keeps its geometry, pair checks and loop energies.
		// Its full energy must still be refreshed (ED and BOTH dangles change).
		buildCandidates(sides[best.left ? 0 : 1], bounds, best.left, last1, last2);
	}
}

//////////////////////////////////////////////////////////////////////////

void
PredictorSeedExtensionKinetic::buildCandidates(SideCandidates & side,
		const Boundary & bounds, const bool left, const size_t last1, const size_t last2) const
{
	side.moves.clear();
	const size_t space1 = std::min(energy.getAccessibility1().getMaxLength() - (bounds[1]-bounds[0]+1),
			left ? bounds[0] : last1-bounds[1]);
	const size_t space2 = std::min(energy.getAccessibility2().getMaxLength() - (bounds[3]-bounds[2]+1),
			left ? bounds[2] : last2-bounds[3]);
	if (space1 == 0 || space2 == 0) {
		return;
	}
	const bool stackOnly = (output.getOutputConstraint().noGUend || !energy.isInternalLoopGUallowed())
			&& energy.isGU(bounds[left ? 0 : 1], bounds[left ? 2 : 3]);
	const size_t maxGap1 = space1 < 2 || stackOnly ? 0 : std::min(energy.getMaxInternalLoopSize1(), space1-2);
	const size_t maxGap2 = space2 < 2 || stackOnly ? 0 : std::min(energy.getMaxInternalLoopSize2(), space2-2);
	side.columns = maxGap2+2;
	side.complementary.assign((maxGap1+2)*side.columns, -1);
	const auto append = [&](size_t s1, size_t s2, bool macro) {
		Candidate c;
		c.left = left; c.s1 = s1; c.s2 = s2; c.macro = macro;
		c.bounds = bounds;
		const size_t pairs = macro ? 2 : 1;
		if (left) {
			c.close1 = bounds[0]-s1-1; c.close2 = bounds[2]-s2-1;
			c.bounds[0] -= s1+pairs; c.bounds[2] -= s2+pairs;
		} else {
			c.close1 = bounds[1]+s1+1; c.close2 = bounds[3]+s2+1;
			c.bounds[1] += s1+pairs; c.bounds[3] += s2+pairs;
		}
		side.moves.push_back(c);
	};
	append(0, 0, false);
	if (space1 >= 2 && space2 >= 2) {
		for (size_t s1 = 0; s1 <= maxGap1; ++s1) {
			for (size_t s2 = 0; s2 <= maxGap2; ++s2) {
				append(s1, s2, true);
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////

const PredictorSeedExtensionKinetic::Candidate *
PredictorSeedExtensionKinetic::updateCandidates(SideCandidates & side,
		const Boundary & bounds, const E_type hybrid, const E_type total) const
{
	// Phase one: each position pair is tested at most once per unchanged end,
	// even when it is the closing pair of one move and outer pair of another.
	size_t stopGap2 = std::numeric_limits<size_t>::max();
	for (Candidate & c : side.moves) {
		c.bounds[c.left ? 1 : 0] = bounds[c.left ? 1 : 0];
		c.bounds[c.left ? 3 : 2] = bounds[c.left ? 3 : 2];
		c.active = c.bounds[1]-c.bounds[0]+1 <= energy.getAccessibility1().getMaxLength()
				&& c.bounds[3]-c.bounds[2]+1 <= energy.getAccessibility2().getMaxLength();
		if (c.active && c.macro) {
			if (c.s2 >= stopGap2) {
				c.active = false;
			} else if (prune(c, bounds)) {
				// A suffix bound rejects this rectangle of larger gaps without
				// further ED, complementarity or loop-energy lookups.
				stopGap2 = c.s2;
				c.active = false;
			}
		}
		if (!c.active || c.topologyKnown) {
			continue;
		}
		const auto complementary = [&](size_t s1, size_t s2) {
			signed char & cached = side.complementary[s1*side.columns+s2];
			if (cached < 0) {
				cached = energy.areComplementary(c.left ? bounds[0]-s1-1 : bounds[1]+s1+1,
						c.left ? bounds[2]-s2-1 : bounds[3]+s2+1);
			}
			return cached != 0;
		};
		c.topologyKnown = true;
		c.topologyAllowed = complementary(c.s1, c.s2)
				&& (!c.macro || complementary(c.s1+1, c.s2+1));
		if (c.topologyAllowed && output.getOutputConstraint().noGUend && (c.s1 != 0 || c.s2 != 0)) {
			c.topologyAllowed = !energy.isGU(c.close1, c.close2);
		}
	}
	const Candidate * best = NULL;
	for (Candidate & c : side.moves) {
		if (!c.active || !c.topologyAllowed) {
			continue;
		}
		if (!c.localKnown) {
			c.localKnown = true;
			c.local = c.left ? energy.getE_interLeft(c.close1, bounds[0], c.close2, bounds[2])
					: energy.getE_interLeft(bounds[1], c.close1, bounds[3], c.close2);
			if (c.macro && E_isNotINF(c.local)) {
				c.local = addEnergy(c.local, c.left
						? energy.getE_interLeft(c.bounds[0], c.close1, c.bounds[2], c.close2)
						: energy.getE_interLeft(c.close1, c.bounds[1], c.close2, c.bounds[3]));
			}
		}
		c.hybrid = addEnergy(hybrid, c.local);
		if (E_isINF(c.hybrid)) {
			continue;
		}
		c.total = energy.getE(c.bounds[0], c.bounds[1], c.bounds[2], c.bounds[3], c.hybrid);
		if (E_isINF(c.total)) {
			continue;
		}
		c.delta = std::int64_t(c.total)-std::int64_t(total);
		if (c.delta < 0 && (best == NULL || isBetter(c, *best))) {
			best = &c;
		}
	}
	return best;
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
	if (size != bestSize) return size < bestSize;
	if (candidate.s1 != best.s1) return candidate.s1 < best.s1;
	// Identical shape/score: retain the shorter move first.
	return !candidate.macro && best.macro;
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
