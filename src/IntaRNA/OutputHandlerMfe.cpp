#include "IntaRNA/OutputHandlerMfe.h"

#include "IntaRNA/OutputHandlerInteractionList.h"
#include "IntaRNA/PredictorMfe2d.h"
#include "IntaRNA/PredictorMfe2dSeed.h"
#include "IntaRNA/SeedHandlerExplicit.h"
#include "IntaRNA/SeedHandlerIdxOffset.h"

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace IntaRNA {
namespace {

const OutputConstraint & requireConstraint(OutputHandler * output)
{
	if (output == NULL) throw std::invalid_argument("OutputHandlerMfe: output is NULL");
	return output->getOutputConstraint();
}

/** Expose fixed-boundary traceback without running a complete prediction. */
class FixedBoundaryMfe : public PredictorMfe2d {
public:
	using PredictorMfe2d::PredictorMfe2d;
	bool trace(Interaction & interaction);
};

bool
FixedBoundaryMfe::trace(Interaction & interaction)
{
	interaction.sort();
	const auto first = interaction.basePairs.front();
	const auto last = interaction.basePairs.back();
	energy.setOffset1(first.first);
	energy.setOffset2(energy.getAccessibility2().getReversedIndex(first.second));
	const size_t n = last.first-first.first+1;
	const size_t m = first.second-last.second+1;
	hybridE_pq.resize(n,m);

	// A singleton has no interior recurrence but still has a complete energy.
	if (n == 1 && m == 1) {
		interaction.basePairs.resize(1);
		interaction.energy = energy.getE(0,0,0,0,energy.getE_init());
		return E_isNotINF(interaction.energy);
	}

	fillHybridE(n-1,m-1,0,0,false);
	if (E_isINF(hybridE_pq(0,0))) return false;
	PredictorMfe2d::traceBack(interaction);
	interaction.energy = energy.getE(0,n-1,0,m-1,hybridE_pq(0,0));
	return E_isNotINF(interaction.energy);
}

/** Fixed-boundary MFE traceback constrained to contain an admitted seed. */
class FixedBoundarySeededMfe : public PredictorMfe2dSeed {
public:
	FixedBoundarySeededMfe(const InteractionEnergy & energy, OutputHandler & output,
			SeedHandler * seedHandler);
	bool trace(Interaction & interaction);
};

FixedBoundarySeededMfe::FixedBoundarySeededMfe(const InteractionEnergy & energy,
		OutputHandler & output, SeedHandler * seedHandler)
	: PredictorMfe2dSeed(energy,output,NULL,seedHandler)
{}

bool
FixedBoundarySeededMfe::trace(Interaction & interaction)
{
	interaction.sort();
	const auto first = interaction.basePairs.front();
	const auto last = interaction.basePairs.back();
	const size_t offset2 = energy.getAccessibility2().getReversedIndex(first.second);
	energy.setOffset1(first.first);
	energy.setOffset2(offset2);
	seedHandler.setOffset1(first.first);
	seedHandler.setOffset2(offset2);
	const size_t n = last.first-first.first+1;
	const size_t m = first.second-last.second+1;
	hybridE_pq.resize(n,m);
	hybridE_pq_seed.resize(n,m);
	if (seedHandler.fillSeed(0,n-1,0,m-1) == 0) return false;

	// Explicit singleton seeds need no interior recurrence.
	if (n == 1 && m == 1) {
		if (!seedHandler.isSeedBound(0,0)) return false;
		const E_type hybrid = seedHandler.getSeedE(0,0)+energy.getE_init();
		interaction.basePairs.resize(1);
		interaction.energy = energy.getE(0,0,0,0,hybrid);
		interaction.setSeedRange(interaction.basePairs.front(),interaction.basePairs.front(),interaction.energy);
		return E_isNotINF(interaction.energy);
	}

	fillHybridE(n-1,m-1,0,0,false);
	if (E_isINF(hybridE_pq_seed(0,0))) return false;
	interaction.energy = energy.getE(0,n-1,0,m-1,hybridE_pq_seed(0,0));
	PredictorMfe2dSeed::traceBack(interaction);
	return E_isNotINF(interaction.energy);
}

// Explicit seeds can contain singletons. The standard seeded MFE recurrence
// represents a singleton inside a larger interaction through an adjacent pair.
class ExplicitMfeSeedConstraint : public SeedConstraint {
public:
	explicit ExplicitMfeSeedConstraint(const SeedConstraint & original);
};

ExplicitMfeSeedConstraint::ExplicitMfeSeedConstraint(const SeedConstraint & original)
	: SeedConstraint(original)
{
	if (!getExplicitSeeds().empty()) bp = 2;
}

} // namespace

OutputHandlerMfe::OutputHandlerMfe(OutputHandler * output,
		const InteractionEnergy & energy, SeedHandler * seedHandler)
	: OutputHandler(requireConstraint(output))
	, output(output)
	, energy(energy)
	, seedHandler(seedHandler)
{}

OutputHandlerMfe::~OutputHandlerMfe() = default;

void
OutputHandlerMfe::traceMfe(Interaction & interaction)
{
	const OutputConstraint & constraint = getOutputConstraint();
	OutputHandlerInteractionList unused(constraint,1);
	if (!seedHandler) {
		FixedBoundaryMfe trace(energy,unused,NULL);
		if (!trace.trace(interaction))
			throw std::runtime_error("ensemble output: no MFE structure in selected site");
		return;
	}

	const auto & explicitSeeds = seedHandler->getConstraint().getExplicitSeeds();
	if (explicitSeeds.empty()) {
		// The predictor owns its seed argument, so interpose a non-owning adapter.
		FixedBoundarySeededMfe trace(energy,unused,
				new SeedHandlerIdxOffset(seedHandler.get(),false));
		if (!trace.trace(interaction))
			throw std::runtime_error("ensemble output: no MFE structure in selected seeded site");
		return;
	}

	// Trace every explicit alternative that occurs within the selected site and
	// retain the lowest-energy representative. Mixed seed lengths and repeated
	// starts therefore remain distinguishable.
	ExplicitMfeSeedConstraint seedConstraint(seedHandler->getConstraint());
	std::vector<std::string> variants;
	const auto first = interaction.basePairs.front();
	const auto last = interaction.basePairs.back();
	if (first == last) {
		FixedBoundarySeededMfe trace(energy,unused,
				new SeedHandlerExplicit(energy,seedConstraint));
		if (!trace.trace(interaction))
			throw std::runtime_error("ensemble output: no MFE structure in selected explicit-seed site");
		return;
	}
	const size_t lo2 = energy.getAccessibility2().getReversedIndex(first.second);
	const size_t hi2 = energy.getAccessibility2().getReversedIndex(last.second);
	std::istringstream encodings(explicitSeeds);
	for (std::string encoding; std::getline(encodings,encoding,','); ) {
		seedConstraint.getExplicitSeeds() = encoding;
		SeedHandlerExplicit seed(energy,seedConstraint);
		size_t a = RnaSequence::lastPos, b = RnaSequence::lastPos;
		if (!seed.updateToNextSeed(a,b,first.first,last.first,lo2,hi2)) continue;
		if (seed.getSeedLength1(a,b) > 1) {
			variants.push_back(encoding);
			continue;
		}
		// Expand a singleton into each feasible adjacent two-pair seed. With
		// noLP the neighbor must stack; otherwise any allowed interior loop works.
		for (int direction : {-1,1}) {
			const size_t max1 = constraint.noLP ? 1 : energy.getMaxInternalLoopSize1()+1;
			const size_t max2 = constraint.noLP ? 1 : energy.getMaxInternalLoopSize2()+1;
			for (size_t d1=1; d1<=max1; ++d1) for (size_t d2=1; d2<=max2; ++d2) {
				if ((direction<0 && (d1>a-first.first || d2>b-lo2))
						|| (direction>0 && (d1>last.first-a || d2>hi2-b))) continue;
				const size_t c = direction<0 ? a-d1 : a+d1;
				const size_t d = direction<0 ? b-d2 : b+d2;
				if (!energy.areComplementary(c,d)) continue;
				variants.push_back(std::to_string(interaction.s1->getInOutIndex(std::min(a,c)))
						+"|"+std::string(d1-1,'.')+"|&"
						+std::to_string(interaction.s2->getInOutIndex(
								energy.getAccessibility2().getReversedIndex(std::max(b,d))))
						+"|"+std::string(d2-1,'.')+"|");
			}
		}
	}

	E_type best = E_INF;
	Interaction representative(interaction);
	for (const auto & encoding : variants) {
		seedConstraint.getExplicitSeeds() = encoding;
		FixedBoundarySeededMfe trace(energy,unused,
				new SeedHandlerExplicit(energy,seedConstraint));
		Interaction candidate(interaction);
		if (trace.trace(candidate) && candidate.energy < best) {
			best = candidate.energy;
			representative = candidate;
		}
	}
	if (E_isINF(best))
		throw std::runtime_error("ensemble output: no MFE structure in selected explicit-seed site");
	interaction = representative;
}

void
OutputHandlerMfe::add(const Interaction & interaction)
{
	if (interaction.isEmpty() || interaction.hasEnsembleEnergy()) {
		output->add(interaction);
		return;
	}
	Interaction representative(interaction);
	representative.ensembleEnergy = interaction.energy;
	// Ensemble predictors can already have traced seed pairs for output. Keep
	// only the selected site boundaries so the fixed-boundary MFE traceback
	// starts from the same canonical input for every predictor variant.
	representative.sort();
	const Interaction::BasePair first = representative.basePairs.front();
	const Interaction::BasePair last = representative.basePairs.back();
	representative.basePairs.clear();
	representative.basePairs.push_back(first);
	if (last != first) representative.basePairs.push_back(last);
	INTARNA_CLEANUP(representative.seed);
	traceMfe(representative);
	output->add(representative);
}

size_t
OutputHandlerMfe::reported() const
{
	return output->reported();
}

void
OutputHandlerMfe::incrementZ(const Z_type subZ)
{
	output->incrementZ(subZ);
}

Z_type
OutputHandlerMfe::getZ() const
{
	return output->getZ();
}

void
OutputHandlerMfe::setExactPartition(bool exact)
{
	OutputHandler::setExactPartition(exact);
	output->setExactPartition(exact);
}

} // namespace IntaRNA
