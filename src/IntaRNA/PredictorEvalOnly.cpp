#include "IntaRNA/PredictorEvalOnly.h"

#include <charconv>
#include <set>
#include <stdexcept>
#include <string_view>

namespace IntaRNA {

namespace {

// Parse one strand without performing arithmetic on unvalidated external indices.
std::vector<size_t> parseStrand( const std::string_view strand, const RnaSequence & sequence )
{
	const size_t structureStart = strand.find_first_of(".|");
	if (structureStart == std::string_view::npos || structureStart == 0) {
		throw std::invalid_argument("--rri: expected start1dotbar1&start2dotbar2, e.g. 1|||&1|||");
	}
	long start = 0;
	const auto parsed = std::from_chars(strand.data(), strand.data()+structureStart, start);
	if (parsed.ec != std::errc() || parsed.ptr != strand.data()+structureStart
			|| start < sequence.getInOutIndex(0)
			|| start > sequence.getInOutIndex(sequence.size()-1)
			|| (start == 0 && sequence.getInOutIndex(0) < 0)) {
		throw std::invalid_argument("--rri: invalid or out-of-range start index");
	}
	const size_t offset = sequence.getIndex(start);
	const auto structure = strand.substr(structureStart);
	if (structure.size() > sequence.size()-offset || structure.find_first_not_of(".|") != std::string_view::npos) {
		throw std::invalid_argument("--rri: dot-bar strand exceeds sequence length or contains invalid symbols");
	}
	std::vector<size_t> paired;
	for (size_t p=0; p<structure.size(); ++p) {
		if (structure[p] == '|') paired.push_back(offset+p);
	}
	return paired;
}

} // namespace

std::vector<Interaction>
PredictorEvalOnly::parseInteractions( const std::string & encoding,
		const RnaSequence & target, const RnaSequence & query )
{
	std::vector<Interaction> result;
	size_t start = 0;
	do {
		const size_t end = encoding.find(':', start);
		const auto entry = std::string_view(encoding).substr(start,
				end == std::string::npos ? end : end-start);
		const size_t separator = entry.find('&');
		if (separator == std::string_view::npos || entry.find('&', separator+1) != std::string_view::npos) {
			throw std::invalid_argument("--rri: each interaction requires exactly one '&' separator");
		}
		const auto paired1 = parseStrand(entry.substr(0, separator), target);
		const auto paired2 = parseStrand(entry.substr(separator+1), query);
		if (paired1.empty() || paired1.size() != paired2.size()) {
			throw std::invalid_argument("--rri: strands must contain the same nonzero number of pairing bars");
		}
		Interaction interaction(target, query);
		for (size_t p=0; p<paired1.size(); ++p) {
			const size_t q = paired2[paired2.size()-1-p];
			if (!RnaSequence::areComplementary(target, query, paired1[p], q)) {
				throw std::invalid_argument("--rri: non-complementary base pair at target "
						+toString(target.getInOutIndex(paired1[p]))+", query "+toString(query.getInOutIndex(q)));
			}
			interaction.basePairs.emplace_back(paired1[p], q);
		}
		result.push_back(interaction);
		if (end == std::string::npos) break;
		start = end+1;
	} while (true);
	return result;
}

PredictorEvalOnly::PredictorEvalOnly( const InteractionEnergy & energy,
		OutputHandler & output, PredictionTracker * predTracker,
		const std::vector<Interaction> & input )
	: Predictor(energy, output, predTracker)
{
	const auto & target = energy.getAccessibility1().getSequence();
	const auto & query = energy.getAccessibility2().getAccessibilityOrigin().getSequence();
	if (input.empty()) throw std::invalid_argument("PredictorEvalOnly: no interactions provided");
	std::set<Interaction::PairingVec> seen;
	for (const auto & interaction : input) {
		if (!interaction.s1 || !interaction.s2 || !interaction.isValid()
				|| interaction.s1->asString() != target.asString()
				|| interaction.s2->asString() != query.asString()) {
			throw std::invalid_argument("PredictorEvalOnly: invalid interaction or incompatible sequences");
		}
		for (const auto & bp : interaction.basePairs) {
			if (bp.first >= target.size() || bp.second >= query.size()
					|| !RnaSequence::areComplementary(target, query, bp.first, bp.second)) {
				throw std::invalid_argument("PredictorEvalOnly: out-of-range or non-complementary base pair");
			}
		}
		if (seen.insert(interaction.basePairs).second) {
			interactions.emplace_back(target, query);
			interactions.back().basePairs = interaction.basePairs;
		}
	}
}

void
PredictorEvalOnly::predict( const IndexRange &, const IndexRange & )
{
	initOptima();
	// Finish evaluation before reporting anything, so an unevaluable structure
	// cannot leave a partially reported list or partition function.
	for (auto & interaction : interactions) {
		E_type hybridE = energy.getE_init();
		for (size_t p=1; p<interaction.basePairs.size(); ++p) {
			const auto & left = interaction.basePairs[p-1];
			const auto & right = interaction.basePairs[p];
			const E_type loopE = energy.getE_interLeft(energy.getIndex1(left), energy.getIndex1(right),
					energy.getIndex2(left), energy.getIndex2(right));
			if (E_isINF(loopE)) {
				throw std::runtime_error("PredictorEvalOnly: interaction loop cannot be evaluated by the selected energy model");
			}
			hybridE += loopE;
		}
		const auto & left = interaction.basePairs.front();
		const auto & right = interaction.basePairs.back();
		interaction.energy = energy.getE(energy.getIndex1(left), energy.getIndex1(right),
				energy.getIndex2(left), energy.getIndex2(right), hybridE);
		if (E_isINF(interaction.energy)) {
			throw std::runtime_error("PredictorEvalOnly: interaction cannot be evaluated with the selected accessibility model; check accessibility windows, constraints and input coverage");
		}
	}
	for (const auto & interaction : interactions) {
		const auto & left = interaction.basePairs.front();
		const auto & right = interaction.basePairs.back();
		updateOptima(energy.getIndex1(left), energy.getIndex1(right),
				energy.getIndex2(left), energy.getIndex2(right), interaction.energy, false, true);
	}
	reportOptima();
}

void
PredictorEvalOnly::initOptima()
{
	Zall = 0;
}

void
PredictorEvalOnly::updateOptima( const size_t i1, const size_t j1,
		const size_t i2, const size_t j2, const E_type interactionE,
		const bool isHybridE, const bool incrementZ )
{
	const E_type totalE = isHybridE ? energy.getE(i1, j1, i2, j2, interactionE) : interactionE;
	if (incrementZ) incrementZall(energy.getBoltzmannWeight(totalE));
	if (predTracker != NULL) predTracker->updateOptimumCalled(i1, j1, i2, j2, totalE);
}

void
PredictorEvalOnly::reportOptima()
{
	output.incrementZ(Zall);
	for (const auto & interaction : interactions) output.add(interaction);
}

} // namespace IntaRNA
