#include "IntaRNA/OutputHandlerSeedCandidates.h"
#include <algorithm>
#include <exception>
#include <map>

namespace IntaRNA {
namespace {
// Compare complete candidates deterministically, including degenerate traces.
bool lessCandidate(const Interaction & a, const Interaction & b)
{
	if (a < b) return true;
	if (b < a) return false;
	return a.basePairs < b.basePairs;
}

bool sameSite(const Interaction & a, const Interaction & b)
{
	return a.basePairs.front() == b.basePairs.front() && a.basePairs.back() == b.basePairs.back();
}
} // namespace

OutputHandlerSeedCandidates::OutputHandlerSeedCandidates(const OutputConstraint & constraint)
	: OutputHandler(constraint)
{
	if (!constraint.reportMaxPerSeed)
		throw std::invalid_argument("per-seed selection requires a positive quota");
}

void OutputHandlerSeedCandidates::add(const Interaction &)
{}

void OutputHandlerSeedCandidates::addSeed(const Interaction::BasePair & seed, const Interaction & interaction)
{
	if (interaction.isEmpty() || E_isINF(interaction.energy) || !outConstraint.reportMax) return;
	// Window workers share this collector. Forward exceptions outside the block.
	std::exception_ptr failure;
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_seedCandidates)
#endif
	{
		try {
			auto & list = candidates[seed];
			const auto existing = std::find_if(list.begin(), list.end(), [&](const Interaction & i) { return sameSite(i,interaction); });
			const bool improves = existing != list.end() ? lessCandidate(interaction,*existing)
					: list.size() < outConstraint.reportMaxPerSeed || lessCandidate(interaction,list.back());
			if (improves) {
				if (existing != list.end()) list.erase(existing);
				const auto position = std::lower_bound(list.begin(),list.end(),interaction,lessCandidate);
				list.insert(position,interaction);
				if (list.size() > outConstraint.reportMaxPerSeed) list.pop_back();
			}
		} catch (...) { failure = std::current_exception(); }
	}
	if (failure) std::rethrow_exception(failure);
}

void OutputHandlerSeedCandidates::report(OutputHandler & destination) const
{
	// Merge equal sites before output filtering, preserving the winning trace.
	using Site = std::pair<Interaction::BasePair,Interaction::BasePair>;
	std::map<Site,Interaction> sites;
	std::vector<Interaction::BasePair> seeds;
	for (const auto & entry : candidates) seeds.push_back(entry.first);
	std::sort(seeds.begin(),seeds.end());
	for (const auto & seed : seeds) for (const auto & interaction : candidates.at(seed)) {
		const Site site(interaction.basePairs.front(),interaction.basePairs.back());
		auto [it,inserted] = sites.try_emplace(site,interaction);
		if (!inserted && lessCandidate(interaction,it->second)) it->second = interaction;
	}
	std::vector<Interaction> ranked;
	for (const auto & site : sites) ranked.push_back(site.second);
	std::sort(ranked.begin(),ranked.end(),lessCandidate);
	if (ranked.empty()) return;

	// Apply final restrictions once to the selected union, without quota refills.
	const auto & c = outConstraint;
	const auto cutoff = std::int64_t(ranked.front().energy) + c.deltaE;
	std::vector<const Interaction *> reported;
	for (const auto & interaction : ranked) {
		if (reported.size() >= c.reportMax || interaction.energy >= c.maxE || interaction.energy > cutoff) break;
		bool overlaps = false;
		for (const auto * previous : reported) {
			const bool target = interaction.basePairs.front().first <= previous->basePairs.back().first
					&& previous->basePairs.front().first <= interaction.basePairs.back().first;
			const bool query = interaction.basePairs.back().second <= previous->basePairs.front().second
					&& previous->basePairs.back().second <= interaction.basePairs.front().second;
			overlaps = overlaps || (target && (c.reportOverlap == OutputConstraint::OVERLAP_NONE || c.reportOverlap == OutputConstraint::OVERLAP_SEQ2))
					|| (query && (c.reportOverlap == OutputConstraint::OVERLAP_NONE || c.reportOverlap == OutputConstraint::OVERLAP_SEQ1));
		}
		if (overlaps) continue;
		destination.add(interaction);
		reported.push_back(&interaction);
	}
}

void OutputHandlerSeedCandidates::clear()
{
	candidates.clear();
}
} // namespace IntaRNA
