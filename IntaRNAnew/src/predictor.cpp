#include "intarnanew/predictor.hpp"

#include "intarnanew/helix_blocks.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <ranges>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <utility>

namespace intarnanew {
namespace {

[[nodiscard]] auto logAdd(const double left, const double right) noexcept -> double {
    if (!std::isfinite(left) && left < 0.0) return right;
    if (!std::isfinite(right) && right < 0.0) return left;
    const auto high = std::max(left, right);
    const auto low = std::min(left, right);
    return high + std::log1p(std::exp(low - high));
}

struct ExplicitSeed {
    std::vector<BasePair> pairs;
};

struct PathState {
    std::vector<BasePair> path;
    Energy hybridEnergy{};
    double logWeight{};
    bool seedMatched{};
    bool lastHasLeftStack{};
    std::vector<SeedMatch> seeds;
};

struct StateKey {
    Index startTarget{};
    Index startQuery{};
    bool seedMatched{};
    bool lastHasLeftStack{};
    Index seedProgress{};
    std::vector<BasePair> suffix;

    friend auto operator==(const StateKey&, const StateKey&) -> bool = default;
};

struct StateKeyHash {
    [[nodiscard]] auto operator()(const StateKey& key) const noexcept -> std::size_t {
        auto value = std::hash<Index>{}(key.startTarget);
        const auto mix = [&](const std::size_t next) {
            value ^= next + 0x9e3779b97f4a7c15ULL + (value << 6U) + (value >> 2U);
        };
        mix(std::hash<Index>{}(key.startQuery));
        mix(std::hash<bool>{}(key.seedMatched));
        mix(std::hash<bool>{}(key.lastHasLeftStack));
        mix(std::hash<Index>{}(key.seedProgress));
        for (const auto& pair : key.suffix) {
            mix(std::hash<Index>{}(pair.target));
            mix(std::hash<Index>{}(pair.query));
        }
        return value;
    }
};

struct SiteKey {
    Index targetBegin{};
    Index targetEnd{};
    Index queryBegin{};
    Index queryEnd{};

    friend auto operator<=>(const SiteKey&, const SiteKey&) = default;
};

struct SiteAccumulator {
    Interaction representative;
    double logWeight{-infinity};
};

[[nodiscard]] auto isStack(const BasePair left, const BasePair right) noexcept -> bool {
    return right.target == left.target + 1U && left.query == right.query + 1U;
}

[[nodiscard]] auto pathSignature(const std::span<const BasePair> path) -> std::string {
    std::string result;
    result.reserve(path.size() * 16U);
    for (const auto pair : path) {
        result += std::to_string(pair.target);
        result.push_back('&');
        result += std::to_string(pair.query);
        result.push_back(',');
    }
    return result;
}

[[nodiscard]] auto decimalLexicalCompare(Index left, Index right) noexcept -> int {
    if (left == right) return 0;
    Index leftDivisor{1U};
    Index rightDivisor{1U};
    while (leftDivisor <= left / 10U) leftDivisor *= 10U;
    while (rightDivisor <= right / 10U) rightDivisor *= 10U;
    while (leftDivisor != 0U && rightDivisor != 0U) {
        const auto leftDigit = left / leftDivisor;
        const auto rightDigit = right / rightDivisor;
        if (leftDigit != rightDigit) return leftDigit < rightDigit ? -1 : 1;
        left %= leftDivisor;
        right %= rightDivisor;
        leftDivisor /= 10U;
        rightDivisor /= 10U;
    }
    return leftDivisor == 0U ? -1 : 1;
}

[[nodiscard]] auto pathSignatureLess(
    const std::span<const BasePair> left,
    const std::span<const BasePair> right) noexcept -> bool {
    const auto common = std::min(left.size(), right.size());
    for (Index index{}; index < common; ++index) {
        if (const auto order = decimalLexicalCompare(left[index].target, right[index].target);
            order != 0) return order < 0;
        if (const auto order = decimalLexicalCompare(left[index].query, right[index].query);
            order != 0) return order < 0;
    }
    return left.size() < right.size();
}

[[nodiscard]] auto effectiveLength(const SideConfig& side) noexcept -> Index {
    Index result = side.interactionLengthMax;
    if (side.accessibility == AccessibilityKind::compute && side.accessibilityWindow != 0U) {
        result = result == 0U ? side.accessibilityWindow : std::min(result, side.accessibilityWindow);
    }
    return result;
}

[[nodiscard]] auto parseSigned(const std::string_view text) -> std::expected<long long, std::string> {
    long long value{};
    const auto [position, error] = std::from_chars(text.data(), text.data() + text.size(), value);
    if (error != std::errc{} || position != text.data() + text.size()) {
        return std::unexpected("invalid sequence coordinate '" + std::string(text) + "'");
    }
    return value;
}

[[nodiscard]] auto splitTopLevel(const std::string_view text, const char delimiter) -> std::vector<std::string_view> {
    std::vector<std::string_view> result;
    std::size_t start{};
    for (std::size_t index = 0; index <= text.size(); ++index) {
        if (index == text.size() || text[index] == delimiter) {
            result.push_back(text.substr(start, index - start));
            start = index + 1U;
        }
    }
    return result;
}

[[nodiscard]] auto explicitSeeds(
    const Sequence& target,
    const Sequence& query,
    const std::string& specification) -> std::expected<std::vector<ExplicitSeed>, std::string> {
    std::vector<ExplicitSeed> result;
    if (specification.empty()) return result;
    for (const auto encoded : splitTopLevel(specification, ',')) {
        const auto ampersand = encoded.find('&');
        if (ampersand == std::string_view::npos) {
            return std::unexpected("explicit seed must contain target&query encodings");
        }
        const auto decodeSide = [](const std::string_view side)
            -> std::expected<std::pair<long long, std::string_view>, std::string> {
            const auto patternStart = side.find_first_of("|.");
            if (patternStart == std::string_view::npos || patternStart == 0U) {
                return std::unexpected("explicit seed side needs a start coordinate and |/. pattern");
            }
            auto coordinate = parseSigned(side.substr(0U, patternStart));
            if (!coordinate) return std::unexpected(coordinate.error());
            const auto pattern = side.substr(patternStart);
            if (!std::ranges::all_of(pattern, [](const char symbol) { return symbol == '|' || symbol == '.'; })) {
                return std::unexpected("explicit seed pattern contains an invalid symbol");
            }
            return std::pair{*coordinate, pattern};
        };
        auto targetSide = decodeSide(encoded.substr(0U, ampersand));
        auto querySide = decodeSide(encoded.substr(ampersand + 1U));
        if (!targetSide) return std::unexpected(targetSide.error());
        if (!querySide) return std::unexpected(querySide.error());

        std::vector<Index> targetBars;
        std::vector<Index> queryBars;
        for (Index offset = 0; offset < targetSide->second.size(); ++offset) {
            if (targetSide->second[offset] == '|') {
                auto position = target.internalIndex(targetSide->first + static_cast<long long>(offset));
                if (!position) return std::unexpected(position.error());
                targetBars.push_back(*position);
            }
        }
        for (Index offset = 0; offset < querySide->second.size(); ++offset) {
            if (querySide->second[offset] == '|') {
                auto position = query.internalIndex(querySide->first + static_cast<long long>(offset));
                if (!position) return std::unexpected(position.error());
                queryBars.push_back(*position);
            }
        }
        if (targetBars.empty() || targetBars.size() != queryBars.size()) {
            return std::unexpected("explicit target/query seed patterns need equal nonzero pair counts");
        }
        ExplicitSeed seed;
        seed.pairs.reserve(targetBars.size());
        for (Index index = 0; index < targetBars.size(); ++index) {
            seed.pairs.push_back({targetBars[index], queryBars[queryBars.size() - index - 1U]});
        }
        result.push_back(std::move(seed));
    }
    return result;
}

[[nodiscard]] auto explicitSeedMatches(
    const std::vector<BasePair>& path,
    const std::vector<ExplicitSeed>& seeds) -> std::vector<SeedMatch> {
    std::vector<SeedMatch> result;
    for (const auto& seed : seeds) {
        if (seed.pairs.size() > path.size()) continue;
        for (Index start = 0; start + seed.pairs.size() <= path.size(); ++start) {
            if (std::equal(seed.pairs.begin(), seed.pairs.end(), path.begin() + static_cast<std::ptrdiff_t>(start))) {
                result.push_back({start, start + seed.pairs.size() - 1U, 0.0});
            }
        }
    }
    return result;
}

void mergeSeedMatches(
    std::vector<SeedMatch>& retained,
    std::vector<SeedMatch> candidates) {
    for (auto& candidate : candidates) {
        const auto duplicate = std::ranges::find_if(retained, [&](const SeedMatch& existing) {
            return existing.firstPair == candidate.firstPair &&
                   existing.lastPair == candidate.lastPair;
        });
        if (duplicate == retained.end()) {
            retained.push_back(std::move(candidate));
        } else if (seedMatchLess(candidate, *duplicate)) {
            *duplicate = std::move(candidate);
        }
    }
    std::ranges::sort(retained, seedMatchLess);
}

[[nodiscard]] auto guEndViolation(
    const Sequence& target,
    const Sequence& query,
    const std::vector<BasePair>& path) noexcept -> bool {
    if (path.empty()) return false;
    if (isGuPair(target[path.front().target], query[path.front().query]) ||
        isGuPair(target[path.back().target], query[path.back().query])) return true;
    for (Index index = 1U; index < path.size(); ++index) {
        if (!isStack(path[index - 1U], path[index]) &&
            (isGuPair(target[path[index - 1U].target], query[path[index - 1U].query]) ||
             isGuPair(target[path[index].target], query[path[index].query]))) return true;
    }
    return false;
}

[[nodiscard]] auto stateKey(
    const PathState& state,
    const Config& config,
    const std::span<const ExplicitSeed> anchoredSeeds) -> StateKey {
    StateKey key{
        state.path.front().target,
        state.path.front().query,
        state.seedMatched,
        state.lastHasLeftStack,
        0U,
        {},
    };
    if (config.model == InteractionModel::helixBlocks) {
        // Helix admissibility is path dependent. Keep paths distinct until the
        // documented block decomposition has been checked.
        key.suffix = state.path;
    } else if (!anchoredSeeds.empty() && !state.seedMatched) {
        // Before an explicit anchor starts, future transitions only depend on
        // the current pair, so one MFE path is sufficient. Once an anchor
        // prefix has started, retain that prefix as the small automaton state.
        // This avoids carrying arbitrary path suffixes while ensuring that a
        // weaker path capable of completing the requested anchor is not
        // merged into an incompatible one.
        Index retained{};
        for (const auto& seed : anchoredSeeds) {
            const auto maximum = std::min(state.path.size(), seed.pairs.size() - 1U);
            for (Index length = maximum; length > retained; --length) {
                if (std::equal(
                        state.path.end() - static_cast<std::ptrdiff_t>(length),
                        state.path.end(), seed.pairs.begin())) {
                    retained = length;
                    break;
                }
            }
        }
        if (retained != 0U) {
            key.suffix.assign(
                state.path.end() - static_cast<std::ptrdiff_t>(retained), state.path.end());
        }
    } else if (config.seed.required && !state.seedMatched) {
        if (config.seed.maxUnpaired == 0U &&
            config.seed.queryMaxUnpaired <= 0 &&
            config.seed.targetMaxUnpaired <= 0) {
            // For the common canonical seed, progress is only the trailing
            // stacked-run length. The current cell already identifies its
            // final pair, so hashing/copying the full suffix is redundant.
            Index progress = state.path.empty() ? 0U : 1U;
            for (Index index = state.path.size(); index > 1U &&
                 progress + 1U < config.seed.basePairs; --index) {
                if (!isStack(state.path[index - 2U], state.path[index - 1U])) break;
                ++progress;
            }
            key.seedProgress = progress;
        } else {
            const auto retained = std::min(state.path.size(), config.seed.basePairs > 0U
                ? config.seed.basePairs - 1U
                : 0U);
            key.suffix.assign(
                state.path.end() - static_cast<std::ptrdiff_t>(retained), state.path.end());
        }
    }
    return key;
}

void mergeState(
    std::unordered_map<StateKey, PathState, StateKeyHash>& states,
    PathState candidate,
    const Config& config,
    const std::span<const ExplicitSeed> anchoredSeeds) {
    auto key = stateKey(candidate, config, anchoredSeeds);
    const auto found = states.find(key);
    if (found == states.end()) {
        states.emplace(std::move(key), std::move(candidate));
        return;
    }
    found->second.logWeight = logAdd(found->second.logWeight, candidate.logWeight);
    if (candidate.hybridEnergy < found->second.hybridEnergy - 1e-12 ||
        (std::abs(candidate.hybridEnergy - found->second.hybridEnergy) <= 1e-12 &&
         pathSignatureLess(candidate.path, found->second.path))) {
        const auto combined = found->second.logWeight;
        found->second = std::move(candidate);
        found->second.logWeight = combined;
    }
}

[[nodiscard]] auto overlapsSelected(
    const Interaction& candidate,
    const std::vector<Interaction>& selected,
    const OverlapPolicy policy) -> bool {
    if (policy == OverlapPolicy::both) return false;
    return std::ranges::any_of(selected, [&](const Interaction& existing) {
        const bool targetOverlap = candidate.targetRange().overlaps(existing.targetRange());
        const bool queryOverlap = candidate.queryRange().overlaps(existing.queryRange());
        if (policy == OverlapPolicy::neither) return targetOverlap || queryOverlap;
        if (policy == OverlapPolicy::target) return queryOverlap;
        return targetOverlap;
    });
}

[[nodiscard]] auto interactionLess(const Interaction& left, const Interaction& right) -> bool {
    const auto energyDifference = left.energy.total() - right.energy.total();
    if (std::abs(energyDifference) > 1e-9) return energyDifference < 0.0;
    const auto lt = left.targetRange();
    const auto rt = right.targetRange();
    const auto lq = left.queryRange();
    const auto rq = right.queryRange();
    return std::tie(lt.begin, lt.end, lq.begin, lq.end) < std::tie(rt.begin, rt.end, rq.begin, rq.end);
}

// Legacy-compatible suboptimal traversal is asymmetric after the global
// optimum was chosen: if target positions must not overlap, equal-energy
// follow-up candidates are visited from the high query coordinate downwards.
// This mirrors rerunning the dynamic program on the remaining target domains,
// while keeping the first optimum and all other policies in canonical order.
[[nodiscard]] auto interactionFollowupLess(
    const Interaction& left,
    const Interaction& right) -> bool {
    const auto energyDifference = left.energy.total() - right.energy.total();
    if (std::abs(energyDifference) > 1e-9) return energyDifference < 0.0;
    const auto lt = left.targetRange();
    const auto rt = right.targetRange();
    if (const auto targetOrder = std::tie(lt.begin, lt.end) <=>
                                 std::tie(rt.begin, rt.end);
        targetOrder != 0) {
        return targetOrder < 0;
    }
    const auto lq = left.queryRange();
    const auto rq = right.queryRange();
    return std::tie(rq.begin, rq.end) < std::tie(lq.begin, lq.end);
}

[[nodiscard]] auto interactionSiteSignature(const Interaction& interaction) -> std::string {
    const auto target = interaction.targetRange();
    const auto query = interaction.queryRange();
    return interaction.targetId + "\n" + interaction.queryId + "\n" +
           std::to_string(target.begin) + ":" + std::to_string(target.end) + "&" +
           std::to_string(query.begin) + ":" + std::to_string(query.end);
}

[[nodiscard]] auto ensembleLess(const Interaction& left, const Interaction& right) -> bool {
    if (std::abs(left.ensembleFreeEnergy - right.ensembleFreeEnergy) > 1e-9) {
        return left.ensembleFreeEnergy < right.ensembleFreeEnergy;
    }
    return interactionLess(left, right);
}

[[nodiscard]] auto ensembleFollowupLess(
    const Interaction& left,
    const Interaction& right) -> bool {
    if (std::abs(left.ensembleFreeEnergy - right.ensembleFreeEnergy) > 1e-9) {
        return left.ensembleFreeEnergy < right.ensembleFreeEnergy;
    }
    return interactionFollowupLess(left, right);
}

void orderForOutput(
    std::vector<const Interaction*>& interactions,
    const Config& config) {
    const auto primaryLess = [&](const Interaction* left, const Interaction* right) {
        return config.model == InteractionModel::ensemble
            ? ensembleLess(*left, *right) : interactionLess(*left, *right);
    };
    std::ranges::sort(interactions, primaryLess);
    if (interactions.size() < 2U ||
        (config.output.overlap != OverlapPolicy::neither &&
         config.output.overlap != OverlapPolicy::query)) {
        return;
    }
    const auto followupLess = [&](const Interaction* left, const Interaction* right) {
        return config.model == InteractionModel::ensemble
            ? ensembleFollowupLess(*left, *right)
            : interactionFollowupLess(*left, *right);
    };
    std::ranges::sort(interactions.begin() + 1, interactions.end(), followupLess);
}

void initializeMonomerEnsembles(
    PredictionResult& result,
    const AccessibilityProvider& targetAccessibility,
    const AccessibilityProvider& queryAccessibility) noexcept {
    if (const auto logPartition = targetAccessibility.ensembleLogPartition()) {
        result.targetLogPartition = *logPartition;
    }
    if (const auto freeEnergy = targetAccessibility.ensembleFreeEnergy()) {
        result.targetEnsembleFreeEnergy = *freeEnergy;
    }
    if (const auto logPartition = queryAccessibility.ensembleLogPartition()) {
        result.queryLogPartition = *logPartition;
    }
    if (const auto freeEnergy = queryAccessibility.ensembleFreeEnergy()) {
        result.queryEnsembleFreeEnergy = *freeEnergy;
    }
}

// For selected-only exact base-pair/model-S output, a site's MFE is just the
// negative maximum pair count plus one constant. A scalar max-plus recurrence
// and bounded top-k therefore preserve all observable results without the
// generic engine's paths, hash states, Boltzmann sums, or all-site materialization.
struct CompactSite {
    Index targetBegin{};
    Index targetEnd{};
    Index queryBegin{};
    Index queryEnd{};
    std::uint16_t pairCount{};
};

[[nodiscard]] auto compactSiteLess(
    const CompactSite& left,
    const CompactSite& right,
    const Energy additiveEnergy) noexcept -> bool {
    const auto leftEnergy = -static_cast<Energy>(left.pairCount) + additiveEnergy;
    const auto rightEnergy = -static_cast<Energy>(right.pairCount) + additiveEnergy;
    const auto energyDifference = leftEnergy - rightEnergy;
    if (std::abs(energyDifference) > 1e-9) return energyDifference < 0.0;
    return std::tie(left.targetBegin, left.targetEnd, left.queryBegin, left.queryEnd) <
           std::tie(right.targetBegin, right.targetEnd, right.queryBegin, right.queryEnd);
}

[[nodiscard]] auto compactSpan(const SideConfig& side, const Index domainSize) noexcept -> Index {
    const auto limit = effectiveLength(side);
    return limit == 0U ? domainSize : std::min(limit, domainSize);
}

[[nodiscard]] auto compactCanPair(const char target, const char query) noexcept -> bool {
    switch (target) {
        case 'A': return query == 'U';
        case 'C': return query == 'G';
        case 'G': return query == 'C' || query == 'U';
        case 'U': return query == 'A' || query == 'G';
        default: return false;
    }
}

[[nodiscard]] auto supportsCompactExactBasePair(
    const Config& config,
    const AccessibilityProvider& targetAccessibility,
    const AccessibilityProvider& queryAccessibility,
    const Interval targetDomain,
    const Interval queryDomain) noexcept -> bool {
    const auto* targetDisabled =
        dynamic_cast<const DisabledAccessibility*>(&targetAccessibility);
    const auto* queryDisabled =
        dynamic_cast<const DisabledAccessibility*>(&queryAccessibility);
    const auto targetSpan = compactSpan(config.target, targetDomain.size());
    const auto querySpan = compactSpan(config.query, queryDomain.size());
    return !config.predictionRequirements.retainAllSites &&
           !config.predictionRequirements.computeInteractionPartition &&
           !config.predictionRequirements.traceback &&
           config.mode == PredictionMode::exact &&
           config.model == InteractionModel::singleSite &&
           config.energy == EnergyKind::basePair &&
           !config.seed.required &&
           config.seed.explicitSeeds.empty() &&
           config.target.accessibility == AccessibilityKind::disabled &&
           config.query.accessibility == AccessibilityKind::disabled &&
           config.target.accessibilityConstraint.empty() &&
           config.query.accessibilityConstraint.empty() &&
           config.target.regions.empty() &&
           config.query.regions.empty() &&
           config.target.regionLengthMax == 0U &&
           config.query.regionLengthMax == 0U &&
           config.windowWidth == 0U &&
           !config.output.noLonelyPairs &&
           !config.output.noGuAtEnds &&
           config.output.minUnpairedProbability <= 0.0 &&
           config.output.overlap == OverlapPolicy::both &&
           targetDisabled != nullptr && queryDisabled != nullptr &&
           targetDisabled->unconstrained(targetDomain) &&
           queryDisabled->unconstrained(queryDomain) &&
           targetSpan != 0U && querySpan != 0U &&
           targetSpan <= std::numeric_limits<std::uint16_t>::max() &&
           querySpan <= std::numeric_limits<std::uint16_t>::max();
}

[[nodiscard]] auto compactExactBasePairPrediction(
    const Config& config,
    const Sequence& target,
    const Sequence& query,
    const Interval targetDomain,
    const Interval queryDomain,
    PredictionResult result) -> PredictionResult {
    const Index resultLimit = config.output.number;
    if (resultLimit == 0U) return result;

    const Index targetSpan = compactSpan(config.target, targetDomain.size());
    const Index querySpan = compactSpan(config.query, queryDomain.size());
    if (targetSpan == 0U || querySpan == 0U) return result;
    const Index targetTransition = config.target.interactionLoopMax >= targetSpan - 1U
        ? targetSpan : config.target.interactionLoopMax + 1U;
    const Index queryTransition = config.query.interactionLoopMax >= querySpan - 1U
        ? querySpan : config.query.interactionLoopMax + 1U;
    const Index scoreRows = targetTransition;
    if (scoreRows > std::numeric_limits<Index>::max() / querySpan) {
        throw std::length_error("compact predictor matrix size overflows addressable memory");
    }

    std::vector<std::uint16_t> scores(scoreRows * querySpan);
    std::vector<std::uint16_t> vertical(querySpan);
    std::vector<Index> deque(querySpan);
    std::vector<CompactSite> retained;
    retained.reserve(std::min<Index>(resultLimit, 1024U));

    const auto totalEnergy = [&](const std::uint16_t pairCount) noexcept -> Energy {
        return -static_cast<Energy>(pairCount) + config.additiveEnergy;
    };
    const auto siteLess = [&](const CompactSite& left, const CompactSite& right) noexcept {
        return compactSiteLess(left, right, config.additiveEnergy);
    };
    const auto retain = [&](const CompactSite candidate) {
        if (!std::isfinite(totalEnergy(candidate.pairCount)) ||
            totalEnergy(candidate.pairCount) >= config.output.maxEnergy - 1e-12) {
            return;
        }
        if (retained.size() == resultLimit &&
            !siteLess(candidate, retained.back())) {
            return;
        }
        const auto position = std::ranges::lower_bound(
            retained, candidate, siteLess);
        retained.insert(position, candidate);
        if (retained.size() > resultLimit) retained.pop_back();
    };

    const Index queryLength = queryDomain.size();
    for (Index targetBegin = targetDomain.begin;
         targetBegin <= targetDomain.end;
         ++targetBegin) {
        const Index maximumTargetOffset =
            std::min(targetSpan - 1U, targetDomain.end - targetBegin);
        for (Index queryBeginReverse = 0U;
             queryBeginReverse < queryLength;
             ++queryBeginReverse) {
            const Index queryEnd = queryDomain.end - queryBeginReverse;
            if (!compactCanPair(target.str()[targetBegin], query.str()[queryEnd])) continue;
            const Index maximumQueryOffset =
                std::min(querySpan - 1U, queryLength - queryBeginReverse - 1U);

            std::ranges::fill(scores, std::uint16_t{});
            scores.front() = 1U;
            retain({targetBegin, targetBegin, queryEnd, queryEnd, 1U});

            for (Index targetOffset = 1U;
                 targetOffset <= maximumTargetOffset;
                 ++targetOffset) {
                const Index firstPreviousTarget =
                    targetOffset > targetTransition
                        ? targetOffset - targetTransition
                        : 0U;
                for (Index queryOffset = 0U;
                     queryOffset <= maximumQueryOffset;
                     ++queryOffset) {
                    std::uint16_t best{};
                    for (Index previousTarget = firstPreviousTarget;
                         previousTarget < targetOffset;
                         ++previousTarget) {
                        best = std::max(
                            best, scores[(previousTarget % scoreRows) * querySpan + queryOffset]);
                    }
                    vertical[queryOffset] = best;
                }

                const Index currentRow = targetOffset % scoreRows;
                std::fill_n(scores.begin() + static_cast<std::ptrdiff_t>(currentRow * querySpan),
                            querySpan, std::uint16_t{});

                Index head{};
                Index tail{};
                for (Index queryOffset = 1U;
                     queryOffset <= maximumQueryOffset;
                     ++queryOffset) {
                    const Index added = queryOffset - 1U;
                    while (tail > head &&
                           vertical[deque[tail - 1U]] <= vertical[added]) {
                        --tail;
                    }
                    deque[tail++] = added;
                    const Index firstPreviousQuery =
                        queryOffset > queryTransition
                            ? queryOffset - queryTransition
                            : 0U;
                    while (head < tail && deque[head] < firstPreviousQuery) {
                        ++head;
                    }
                    if (head == tail || vertical[deque[head]] == 0U) continue;

                    const Index targetEnd = targetBegin + targetOffset;
                    const Index queryBegin =
                        queryDomain.end - (queryBeginReverse + queryOffset);
                    if (!compactCanPair(target.str()[targetEnd], query.str()[queryBegin])) continue;
                    const auto pairCount = static_cast<std::uint16_t>(
                        vertical[deque[head]] + 1U);
                    scores[currentRow * querySpan + queryOffset] = pairCount;
                    retain({
                        targetBegin, targetEnd, queryBegin, queryEnd, pairCount,
                    });
                }
            }
        }
    }

    if (retained.empty()) return result;
    const auto minimumEnergy = totalEnergy(retained.front().pairCount);
    result.interactions.reserve(retained.size());
    for (const auto& site : retained) {
        if (totalEnergy(site.pairCount) >
            minimumEnergy + config.output.deltaEnergy + 1e-9) {
            continue;
        }
        Interaction interaction;
        interaction.targetId = target.id();
        interaction.queryId = query.id();
        interaction.pairs.push_back({site.targetBegin, site.queryEnd});
        if (site.targetBegin != site.targetEnd) {
            interaction.pairs.push_back({site.targetEnd, site.queryBegin});
        }
        interaction.energy.initiation = -1.0;
        interaction.energy.loops =
            -static_cast<Energy>(site.pairCount - 1U);
        interaction.energy.additive = config.additiveEnergy;
        interaction.ensembleFreeEnergy = interaction.energy.total();
        result.interactions.push_back(std::move(interaction));
    }
    return result;
}

} // namespace

auto parseIntervals(const Sequence& sequence, const std::string& specification)
    -> std::expected<std::vector<Interval>, std::string> {
    if (sequence.empty()) return std::unexpected("cannot parse regions for an empty sequence");
    if (specification.empty()) return std::vector<Interval>{{0U, sequence.size() - 1U}};
    std::vector<Interval> result;
    std::size_t start{};
    while (start < specification.size()) {
        const auto comma = specification.find(',', start);
        auto token = std::string_view(specification).substr(
            start, comma == std::string::npos ? std::string::npos : comma - start);
        while (!token.empty() && std::isspace(static_cast<unsigned char>(token.front())) != 0) {
            token.remove_prefix(1U);
        }
        while (!token.empty() && std::isspace(static_cast<unsigned char>(token.back())) != 0) {
            token.remove_suffix(1U);
        }
        if (token.empty()) return std::unexpected("region list contains an empty entry");
        std::optional<std::size_t> separator;
        for (std::size_t position = 1U; position < token.size(); ++position) {
            if (token[position] == '-') {
                separator = position;
                break;
            }
        }
        if (!separator) return std::unexpected("region must use FROM-TO encoding");
        auto first = parseSigned(token.substr(0U, *separator));
        auto last = parseSigned(token.substr(*separator + 1U));
        if (!first || !last || *first > *last) return std::unexpected("invalid region coordinate range");
        auto internalFirst = sequence.internalIndex(*first);
        auto internalLast = sequence.internalIndex(*last);
        if (!internalFirst || !internalLast) return std::unexpected("region is outside the sequence");
        result.push_back({*internalFirst, *internalLast});
        if (comma == std::string::npos) break;
        start = comma + 1U;
    }
    if (start == specification.size()) return std::unexpected("region list contains an empty entry");
    std::ranges::sort(result, [](const Interval left, const Interval right) {
        return std::tie(left.begin, left.end) < std::tie(right.begin, right.end);
    });
    result.erase(std::unique(result.begin(), result.end()), result.end());
    for (Index index = 1U; index < result.size(); ++index) {
        if (result[index - 1U].overlaps(result[index])) {
            return std::unexpected("overlapping regions are ambiguous; use disjoint intervals");
        }
    }
    return result;
}

auto configuredRegions(const Sequence& sequence, const SideConfig& config)
    -> std::expected<std::vector<Interval>, std::string> {
    if (!config.regions.empty() && config.regionLengthMax != 0U) {
        return std::unexpected("explicit and automatic regions are mutually exclusive");
    }
    return parseIntervals(sequence, config.regions);
}

auto decomposeAccessibleRegions(
    const std::span<const Interval> input,
    const std::size_t regionLengthMax,
    const std::size_t seedLength,
    const AccessibilityProvider& accessibility) -> std::vector<Interval> {
    if (regionLengthMax == 0U) return {input.begin(), input.end()};
    if (seedLength == 0U) throw std::invalid_argument("automatic region decomposition needs a positive seed length");

    std::vector<Interval> pending(input.begin(), input.end());
    std::vector<Interval> result;
    while (!pending.empty()) {
        const auto current = pending.back();
        pending.pop_back();
        if (current.size() < seedLength) continue;
        if (current.size() <= regionLengthMax) {
            result.push_back(current);
            continue;
        }

        Index cutBegin = current.begin;
        Energy highestOpening = -infinity;
        for (Index candidate = current.begin;
             candidate <= current.end + 1U - seedLength;
             ++candidate) {
            const auto opening = accessibility.openingEnergy(
                {candidate, candidate + seedLength - 1U});
            if (opening > highestOpening + 1e-12) {
                highestOpening = opening;
                cutBegin = candidate;
            }
        }
        const Index cutEnd = cutBegin + seedLength - 1U;
        if (cutEnd < current.end) pending.push_back({cutEnd + 1U, current.end});
        if (cutBegin > current.begin) pending.push_back({current.begin, cutBegin - 1U});
    }
    std::ranges::sort(result, {}, &Interval::begin);
    return result;
}

auto decomposeWindows(
    const Interval parent,
    const std::size_t width,
    const std::size_t overlap) -> std::vector<Interval> {
    if (parent.size() == 0U) return {};
    if (width == 0U || parent.size() <= width) return {parent};
    if (overlap >= width) throw std::invalid_argument("window overlap has to be smaller than its width");

    const Index step = width - overlap;
    std::vector<Interval> result;
    for (Index begin = parent.begin;;) {
        const auto end = std::min(parent.end, begin + width - 1U);
        result.push_back({begin, end});
        if (end == parent.end) break;
        begin += step;
    }
    return result;
}

auto reducePredictions(
    const std::span<const PredictionResult> predictions,
    const Config& config) -> PredictionResult {
    PredictionResult result;
    if (predictions.empty()) return result;
    result.rt = predictions.front().rt;
    result.targetLogPartition = predictions.front().targetLogPartition;
    result.queryLogPartition = predictions.front().queryLogPartition;
    result.targetEnsembleFreeEnergy = predictions.front().targetEnsembleFreeEnergy;
    result.queryEnsembleFreeEnergy = predictions.front().queryEnsembleFreeEnergy;

    std::unordered_map<std::string, Interaction> unique;
    for (const auto& prediction : predictions) {
        if (std::abs(prediction.rt - result.rt) > 1e-12) {
            throw std::invalid_argument("cannot reduce predictions with different RT values");
        }
        for (const auto& interaction : prediction.ensembleSites) {
            const auto signature = interactionSiteSignature(interaction);
            const auto found = unique.find(signature);
            if (found == unique.end()) {
                unique.insert_or_assign(signature, interaction);
            } else {
                const auto siteEnergy = std::min(
                    found->second.ensembleFreeEnergy, interaction.ensembleFreeEnergy);
                if (interactionLess(interaction, found->second)) found->second = interaction;
                found->second.ensembleFreeEnergy = siteEnergy;
            }
        }
    }

    result.ensembleSites.reserve(unique.size());
    for (auto& [signature, interaction] : unique) {
        static_cast<void>(signature);
        result.ensembleSites.push_back(std::move(interaction));
    }
    unique.clear();
    std::ranges::sort(result.ensembleSites, interactionLess);

    result.logPartition = -infinity;
    result.ensembleFreeEnergy = infinity;
    const bool seededExtensionWithoutPartition =
        config.model == InteractionModel::seedExtension &&
        (config.seed.required || !config.seed.explicitSeeds.empty());
    if (!seededExtensionWithoutPartition) {
        for (const auto& interaction : result.ensembleSites) {
            if (std::isfinite(interaction.ensembleFreeEnergy)) {
                result.logPartition = logAdd(
                    result.logPartition, -interaction.ensembleFreeEnergy / result.rt);
            }
        }
    }
    if (std::isfinite(result.logPartition)) {
        result.ensembleFreeEnergy = -result.rt * result.logPartition;
        for (auto& interaction : result.ensembleSites) {
            interaction.probability = std::isfinite(interaction.ensembleFreeEnergy)
                ? std::exp(-interaction.ensembleFreeEnergy / result.rt - result.logPartition)
                : 0.0;
        }
    }
    if (result.ensembleSites.empty() || config.output.number == 0U) return result;

    std::vector<const Interaction*> ranked;
    ranked.reserve(result.ensembleSites.size());
    for (const auto& interaction : result.ensembleSites) ranked.push_back(&interaction);
    orderForOutput(ranked, config);
    const auto minimumEnergy = ranked.front()->energy.total();
    for (const auto* interaction : ranked) {
        if (interaction->energy.total() > minimumEnergy + config.output.deltaEnergy + 1e-9) continue;
        if (overlapsSelected(*interaction, result.interactions, config.output.overlap)) continue;
        result.interactions.push_back(*interaction);
        if (result.interactions.size() >= config.output.number) break;
    }
    return result;
}

Predictor::Predictor(Config config)
    : config_(std::move(config)), energy_(makeEnergyModel(config_)) {}

auto Predictor::predict(
    const Sequence& target,
    const Sequence& query,
    const AccessibilityProvider& targetAccessibility,
    const AccessibilityProvider& queryAccessibility) const -> PredictionResult {
    if (target.empty() || query.empty()) {
        PredictionResult empty;
        empty.rt = energy_->rt();
        initializeMonomerEnsembles(empty, targetAccessibility, queryAccessibility);
        return empty;
    }
    auto targetRegionsResult = configuredRegions(target, config_.target);
    auto queryRegionsResult = configuredRegions(query, config_.query);
    if (!targetRegionsResult) throw std::invalid_argument(targetRegionsResult.error());
    if (!queryRegionsResult) throw std::invalid_argument(queryRegionsResult.error());
    auto targetRegions = decomposeAccessibleRegions(
        *targetRegionsResult, config_.target.regionLengthMax, config_.seed.basePairs,
        targetAccessibility);
    auto queryRegions = decomposeAccessibleRegions(
        *queryRegionsResult, config_.query.regionLengthMax, config_.seed.basePairs,
        queryAccessibility);

    std::vector<PredictionResult> predictions;
    predictions.reserve(targetRegions.size() * queryRegions.size());
    for (const auto targetRegion : targetRegions) {
        for (const auto queryRegion : queryRegions) {
            predictions.push_back(predict(target, query, targetAccessibility, queryAccessibility,
                                          targetRegion, queryRegion));
        }
    }
    if (predictions.empty()) {
        PredictionResult empty;
        empty.rt = energy_->rt();
        initializeMonomerEnsembles(empty, targetAccessibility, queryAccessibility);
        return empty;
    }
    if (predictions.size() == 1U) return std::move(predictions.front());
    return reducePredictions(predictions, config_);
}

auto Predictor::predict(
    const Sequence& target,
    const Sequence& query,
    const AccessibilityProvider& targetAccessibility,
    const AccessibilityProvider& queryAccessibility,
    const Interval targetDomain,
    const Interval queryDomain) const -> PredictionResult {
    PredictionResult result;
    result.rt = energy_->rt();
    initializeMonomerEnsembles(result, targetAccessibility, queryAccessibility);
    if (target.empty() || query.empty()) return result;
    if (targetDomain.begin > targetDomain.end || targetDomain.end >= target.size() ||
        queryDomain.begin > queryDomain.end || queryDomain.end >= query.size()) {
        throw std::out_of_range("prediction domain is outside the input sequence");
    }

    if (supportsCompactExactBasePair(
            config_, targetAccessibility, queryAccessibility,
            targetDomain, queryDomain)) {
        return compactExactBasePairPrediction(
            config_, target, query, targetDomain, queryDomain, std::move(result));
    }
    auto explicitSeedResult = explicitSeeds(target, query, config_.seed.explicitSeeds);
    if (!explicitSeedResult) throw std::invalid_argument(explicitSeedResult.error());
    const auto& anchoredSeeds = *explicitSeedResult;

    const Index targetLengthLimit = effectiveLength(config_.target);
    const Index queryLengthLimit = effectiveLength(config_.query);
    const Index targetDomainLength = targetDomain.size();
    const Index queryLength = queryDomain.size();
    if (targetDomainLength > std::numeric_limits<Index>::max() / queryLength) {
        throw std::length_error("prediction domain matrix size overflows addressable memory");
    }
    const Index cellCount = targetDomainLength * queryLength;
    std::vector<std::vector<PathState>> states(cellCount);
    const auto cell = [queryLength, targetBegin = targetDomain.begin](
                          const Index targetIndex, const Index queryReverse) {
        return (targetIndex - targetBegin) * queryLength + queryReverse;
    };

    const auto seedsForPath = [&](const std::vector<BasePair>& path) -> std::vector<SeedMatch> {
        if (!anchoredSeeds.empty()) {
            auto matches = explicitSeedMatches(path, anchoredSeeds);
            for (auto& match : matches) {
                auto breakdown = energy_->evaluate(target, query, std::span(path).subspan(
                    match.firstPair, match.lastPair - match.firstPair + 1U));
                const auto& first = path[match.firstPair];
                const auto& last = path[match.lastPair];
                const Interval targetInterval{first.target, last.target};
                const Interval queryInterval{last.query, first.query};
                breakdown.openingTarget = targetAccessibility.openingEnergy(targetInterval);
                breakdown.openingQuery = queryAccessibility.openingEnergy(queryInterval);
                match.energy = breakdown.total();
                match.openingTarget = breakdown.openingTarget;
                match.openingQuery = breakdown.openingQuery;
                match.unpairedTarget = targetAccessibility.unpairedProbability(targetInterval);
                match.unpairedQuery = queryAccessibility.unpairedProbability(queryInterval);
            }
            std::vector<SeedMatch> normalized;
            mergeSeedMatches(normalized, std::move(matches));
            return normalized;
        }
        if (!config_.seed.required || path.size() < config_.seed.basePairs) return {};
        const Index begin = path.size() - config_.seed.basePairs;
        const Index end = path.size() - 1U;
        const auto targetGap = path[end].target - path[begin].target + 1U - config_.seed.basePairs;
        const auto queryGap = path[begin].query - path[end].query + 1U - config_.seed.basePairs;
        const auto queryMax = config_.seed.queryMaxUnpaired < 0
            ? config_.seed.maxUnpaired
            : static_cast<Index>(config_.seed.queryMaxUnpaired);
        const auto targetMax = config_.seed.targetMaxUnpaired < 0
            ? config_.seed.maxUnpaired
            : static_cast<Index>(config_.seed.targetMaxUnpaired);
        if (targetGap > targetMax || queryGap > queryMax ||
            targetGap + queryGap > config_.seed.maxUnpaired) return {};
        for (Index index = begin; index <= end; ++index) {
            if (config_.seed.noGu && isGuPair(target[path[index].target], query[path[index].query])) {
                return {};
            }
        }
        if (config_.seed.noGuAtEnds &&
            (isGuPair(target[path[begin].target], query[path[begin].query]) ||
             isGuPair(target[path[end].target], query[path[end].query]))) return {};

        if (!config_.seed.targetRanges.empty()) {
            auto ranges = parseIntervals(target, config_.seed.targetRanges);
            if (!ranges || !std::ranges::any_of(*ranges, [&](const Interval range) {
                return range.contains(path[begin].target) && range.contains(path[end].target);
            })) return {};
        }
        if (!config_.seed.queryRanges.empty()) {
            auto ranges = parseIntervals(query, config_.seed.queryRanges);
            if (!ranges || !std::ranges::any_of(*ranges, [&](const Interval range) {
                return range.contains(path[end].query) && range.contains(path[begin].query);
            })) return {};
        }

        const auto slice = std::span(path).subspan(begin, config_.seed.basePairs);
        auto breakdown = energy_->evaluate(target, query, slice);
        const Interval targetInterval{path[begin].target, path[end].target};
        const Interval queryInterval{path[end].query, path[begin].query};
        breakdown.openingTarget = targetAccessibility.openingEnergy(targetInterval);
        breakdown.openingQuery = queryAccessibility.openingEnergy(queryInterval);
        const auto targetProbability = targetAccessibility.unpairedProbability(targetInterval);
        const auto queryProbability = queryAccessibility.unpairedProbability(queryInterval);
        if (breakdown.hybrid() > config_.seed.maxHybridEnergy + 1e-9 ||
            breakdown.total() > config_.seed.maxEnergy + 1e-9 ||
            targetProbability + 1e-15 < config_.seed.minUnpairedProbability ||
            queryProbability + 1e-15 < config_.seed.minUnpairedProbability) {
            return {};
        }
        return {SeedMatch{begin, end, breakdown.total(), breakdown.openingTarget,
                          breakdown.openingQuery, targetProbability, queryProbability}};
    };

    const bool seedRequired = config_.seed.required || !anchoredSeeds.empty();
    const auto maxStates = config_.mode == PredictionMode::heuristic
        ? 96U
        : std::numeric_limits<std::size_t>::max();
    const bool heuristicSeedExtension =
        config_.mode == PredictionMode::heuristic &&
        config_.model == InteractionModel::seedExtension && seedRequired;
    std::vector<PathState> heuristicSeedStates;
    std::set<std::string> heuristicSeedPaths;
    const auto preserveSeedState = [&](const PathState& state) {
        if (anchoredSeeds.empty()) return false;
        if (std::ranges::any_of(state.seeds, [&](const SeedMatch& seed) {
                return seed.firstPair == 0U ||
                       seed.lastPair + 1U == state.path.size();
            })) {
            return true;
        }
        if (state.seedMatched) return false;
        const auto& start = state.path.front();
        const auto& currentPair = state.path.back();
        return std::ranges::any_of(anchoredSeeds, [&](const ExplicitSeed& seed) {
            const auto& seedFirst = seed.pairs.front();
            const auto& seedLast = seed.pairs.back();
            if (start.target > seedFirst.target || start.query < seedFirst.query ||
                currentPair.target > seedLast.target || currentPair.query < seedLast.query) {
                return false;
            }
            return (targetLengthLimit == 0U ||
                    seedLast.target - start.target + 1U <= targetLengthLimit) &&
                   (queryLengthLimit == 0U ||
                    start.query - seedLast.query + 1U <= queryLengthLimit);
        });
    };

    for (Index targetIndex = targetDomain.begin; targetIndex <= targetDomain.end; ++targetIndex) {
        if (targetAccessibility.blocked(targetIndex)) continue;
        for (Index queryReverse = 0U; queryReverse < queryLength; ++queryReverse) {
            const Index queryIndex = queryDomain.end - queryReverse;
            if (queryAccessibility.blocked(queryIndex) || !canPair(target[targetIndex], query[queryIndex])) continue;

            std::unordered_map<StateKey, PathState, StateKeyHash> current;
            PathState initial{{{targetIndex, queryIndex}}, energy_->initiationEnergy(),
                              -energy_->initiationEnergy() / energy_->rt(), !seedRequired, false, {}};
            initial.seeds = seedsForPath(initial.path);
            if (!initial.seeds.empty()) {
                initial.seedMatched = true;
            }
            if (!config_.output.noGuAtEnds ||
                !isGuPair(target[targetIndex], query[queryIndex])) {
                mergeState(current, std::move(initial), config_, anchoredSeeds);
            }

            const Index targetDistanceMax = std::min(
                targetIndex - targetDomain.begin, config_.target.interactionLoopMax + 1U);
            const Index queryDistanceMax = std::min(queryReverse, config_.query.interactionLoopMax + 1U);
            for (Index targetDistance = 1U; targetDistance <= targetDistanceMax; ++targetDistance) {
                const Index previousTarget = targetIndex - targetDistance;
                for (Index queryDistance = 1U; queryDistance <= queryDistanceMax; ++queryDistance) {
                    const Index previousQueryReverse = queryReverse - queryDistance;
                    const Index previousQuery = queryDomain.end - previousQueryReverse;
                    const auto transition = energy_->transitionEnergy(
                        target, query, {previousTarget, previousQuery}, {targetIndex, queryIndex});
                    for (const auto& previous : states[cell(previousTarget, previousQueryReverse)]) {
                        if (targetLengthLimit != 0U && targetIndex - previous.path.front().target + 1U > targetLengthLimit) continue;
                        if (queryLengthLimit != 0U && previous.path.front().query - queryIndex + 1U > queryLengthLimit) continue;
                        const bool stacked = isStack(previous.path.back(), {targetIndex, queryIndex});
                        if (config_.output.noGuAtEnds && !stacked &&
                            (isGuPair(target[previous.path.back().target], query[previous.path.back().query]) ||
                             isGuPair(target[targetIndex], query[queryIndex]))) continue;
                        if (config_.output.noLonelyPairs && !stacked &&
                            !previous.lastHasLeftStack) continue;

                        PathState next = previous;
                        next.path.push_back({targetIndex, queryIndex});
                        next.hybridEnergy += transition;
                        next.logWeight -= transition / energy_->rt();
                        next.lastHasLeftStack = stacked;
                        auto seeds = seedsForPath(next.path);
                        if (!seeds.empty()) {
                            next.seedMatched = true;
                            mergeSeedMatches(next.seeds, std::move(seeds));
                        }
                        mergeState(current, std::move(next), config_, anchoredSeeds);
                    }
                }
            }

            auto& destination = states[cell(targetIndex, queryReverse)];
            destination.reserve(current.size());
            for (auto& [key, state] : current) {
                static_cast<void>(key);
                destination.push_back(std::move(state));
            }
            if (heuristicSeedExtension) {
                for (const auto& state : destination) {
                    const bool oneSided = std::ranges::any_of(
                        state.seeds, [&](const SeedMatch& seed) {
                            return seed.firstPair == 0U ||
                                   seed.lastPair + 1U == state.path.size();
                        });
                    if (oneSided && heuristicSeedPaths.insert(pathSignature(state.path)).second) {
                        heuristicSeedStates.push_back(state);
                    }
                }
            }
            if (destination.size() > maxStates) {
                std::vector<PathState> protectedStates;
                std::vector<PathState> candidates;
                std::map<std::tuple<Index, Index, bool>, PathState> boundaryBest;
                protectedStates.reserve(destination.size());
                candidates.reserve(destination.size());
                for (auto& state : destination) {
                    if (preserveSeedState(state)) {
                        protectedStates.push_back(std::move(state));
                        continue;
                    }
                    if (!anchoredSeeds.empty()) {
                        candidates.push_back(std::move(state));
                        continue;
                    }
                    const auto boundary = std::tuple{
                        state.path.front().target,
                        state.path.front().query,
                        state.seedMatched,
                    };
                    auto found = boundaryBest.find(boundary);
                    if (found == boundaryBest.end() ||
                        state.hybridEnergy < found->second.hybridEnergy - 1e-12 ||
                        (std::abs(state.hybridEnergy - found->second.hybridEnergy) <= 1e-12 &&
                         pathSignatureLess(state.path, found->second.path))) {
                        if (found != boundaryBest.end()) {
                            candidates.push_back(std::move(found->second));
                            found->second = std::move(state);
                        } else {
                            boundaryBest.emplace(boundary, std::move(state));
                        }
                    } else {
                        candidates.push_back(std::move(state));
                    }
                }
                for (auto& [boundary, state] : boundaryBest) {
                    static_cast<void>(boundary);
                    protectedStates.push_back(std::move(state));
                }
                const auto retained = std::min<std::size_t>(candidates.size(), maxStates);
                std::ranges::partial_sort(candidates,
                    candidates.begin() + static_cast<std::ptrdiff_t>(retained),
                    [](const PathState& left, const PathState& right) {
                        if (left.seedMatched != right.seedMatched) return left.seedMatched > right.seedMatched;
                        return left.hybridEnergy < right.hybridEnergy;
                    });
                candidates.resize(retained);
                destination = std::move(protectedStates);
                destination.insert(destination.end(),
                    std::make_move_iterator(candidates.begin()),
                    std::make_move_iterator(candidates.end()));
            }
        }
    }

    if (heuristicSeedExtension) {
        struct RightExtension {
            Energy transitionEnergy{infinity};
            std::vector<BasePair> suffix;
            std::vector<SeedMatch> seeds;
        };
        std::map<std::string, RightExtension> bestRight;
        for (const auto& state : heuristicSeedStates) {
            for (const auto& seed : state.seeds) {
                if (seed.firstPair != 0U || seed.lastPair >= state.path.size()) continue;
                const auto seedPath = std::span(state.path).subspan(
                    seed.firstPair, seed.lastPair - seed.firstPair + 1U);
                const auto key = pathSignature(seedPath);
                Energy extension{};
                for (Index index = seed.lastPair + 1U; index < state.path.size(); ++index) {
                    extension += energy_->transitionEnergy(
                        target, query, state.path[index - 1U], state.path[index]);
                }
                std::vector<BasePair> suffix(
                    state.path.begin() + static_cast<std::ptrdiff_t>(seed.lastPair),
                    state.path.end());
                auto [found, inserted] = bestRight.try_emplace(
                    key, RightExtension{extension, suffix, state.seeds});
                if (!inserted &&
                    (extension < found->second.transitionEnergy - 1e-12 ||
                     (std::abs(extension - found->second.transitionEnergy) <= 1e-12 &&
                      pathSignatureLess(suffix, found->second.suffix)))) {
                    found->second = RightExtension{
                        extension, std::move(suffix), state.seeds};
                }
            }
        }

        std::vector<PathState> synthesized;
        for (const auto& left : heuristicSeedStates) {
            for (const auto& seed : left.seeds) {
                if (seed.lastPair + 1U != left.path.size()) continue;
                const auto seedPath = std::span(left.path).subspan(
                    seed.firstPair, seed.lastPair - seed.firstPair + 1U);
                const auto found = bestRight.find(pathSignature(seedPath));
                if (found == bestRight.end() || found->second.suffix.size() <= 1U) continue;

                std::vector<BasePair> path = left.path;
                path.insert(path.end(), found->second.suffix.begin() + 1,
                            found->second.suffix.end());
                if ((targetLengthLimit != 0U &&
                     path.back().target - path.front().target + 1U > targetLengthLimit) ||
                    (queryLengthLimit != 0U &&
                     path.front().query - path.back().query + 1U > queryLengthLimit)) {
                    continue;
                }
                heuristicSeedPaths.insert(pathSignature(path));

                PathState combined;
                combined.path = std::move(path);
                combined.hybridEnergy = energy_->initiationEnergy();
                for (Index index = 1U; index < combined.path.size(); ++index) {
                    combined.hybridEnergy += energy_->transitionEnergy(
                        target, query, combined.path[index - 1U], combined.path[index]);
                }
                combined.logWeight = -combined.hybridEnergy / energy_->rt();
                combined.seedMatched = true;
                combined.lastHasLeftStack = combined.path.size() > 1U &&
                    isStack(combined.path[combined.path.size() - 2U], combined.path.back());
                combined.seeds = left.seeds;
                auto rightSeeds = found->second.seeds;
                for (auto& rightSeed : rightSeeds) {
                    rightSeed.firstPair += seed.firstPair;
                    rightSeed.lastPair += seed.firstPair;
                }
                mergeSeedMatches(combined.seeds, std::move(rightSeeds));
                if (!combined.seeds.empty()) synthesized.push_back(std::move(combined));
            }
        }
        heuristicSeedStates.insert(
            heuristicSeedStates.end(),
            std::make_move_iterator(synthesized.begin()),
            std::make_move_iterator(synthesized.end()));

        // Only the independently generated seed-extension states belong to H
        // output. Release the propagation matrix before evaluating sites.
        states.clear();
        states.resize(1U);
        states.front() = std::move(heuristicSeedStates);
    }

    std::map<SiteKey, SiteAccumulator> sites;
    std::set<std::string> seedOnlySeen;
    for (auto& cellStates : states) {
        for (auto& state : cellStates) {
            if (seedRequired && !state.seedMatched) continue;
            if (config_.output.noLonelyPairs && !state.lastHasLeftStack) continue;
            auto path = std::move(state.path);
            auto seeds = std::move(state.seeds);
            double stateLogWeight = state.logWeight;
            if (config_.mode == PredictionMode::seedOnly) {
                if (seeds.empty()) continue;
                SeedMatch seed = seeds.front();
                std::vector<BasePair> seedPath(
                    path.begin() + static_cast<std::ptrdiff_t>(seed.firstPair),
                    path.begin() + static_cast<std::ptrdiff_t>(seed.lastPair + 1U));
                path = std::move(seedPath);
                const auto signature = pathSignature(path);
                if (!seedOnlySeen.insert(signature).second) continue;
                seed.firstPair = 0U;
                seed.lastPair = path.size() - 1U;
                seeds = {seed};
                const auto seedEnergy = energy_->evaluate(target, query, path);
                stateLogWeight = -(seedEnergy.initiation + seedEnergy.loops) / energy_->rt();
            }
            if (config_.model == InteractionModel::helixBlocks &&
                decomposeHelixBlocks(
                    target, query, path, config_.helix, *energy_,
                    targetAccessibility, queryAccessibility).empty()) continue;
            if (config_.output.noGuAtEnds && guEndViolation(target, query, path)) continue;

            Interaction interaction;
            interaction.targetId = target.id();
            interaction.queryId = query.id();
            interaction.pairs = std::move(path);
            interaction.energy = energy_->evaluate(target, query, interaction.pairs);
            // Terminal dangles require one additional exterior position on
            // both interaction sites. If either configured site-length bound
            // is already saturated, that exterior context is outside the DP
            // domain and neither terminal dangle contributes. Raising a bound
            // by one makes the same traceback eligible again.
            if ((targetLengthLimit != 0U &&
                 interaction.targetRange().size() == targetLengthLimit) ||
                (queryLengthLimit != 0U &&
                 interaction.queryRange().size() == queryLengthLimit)) {
                interaction.energy.dangleLeft = 0.0;
                interaction.energy.dangleRight = 0.0;
            }
            interaction.energy.openingTarget = targetAccessibility.openingEnergy(interaction.targetRange());
            interaction.energy.openingQuery = queryAccessibility.openingEnergy(interaction.queryRange());
            interaction.energy.additive = config_.additiveEnergy;
            interaction.unpairedTarget =
                targetAccessibility.unpairedProbability(interaction.targetRange());
            interaction.unpairedQuery =
                queryAccessibility.unpairedProbability(interaction.queryRange());
            interaction.seeds = std::move(seeds);
            if (!std::isfinite(interaction.energy.total()) ||
                interaction.energy.total() >= config_.output.maxEnergy - 1e-12) continue;

            bool probabilityRejected = false;
            for (Index position = interaction.targetRange().begin; position <= interaction.targetRange().end; ++position) {
                if (targetAccessibility.positionUnpairedProbability(position) + 1e-15 <
                    config_.output.minUnpairedProbability) probabilityRejected = true;
            }
            for (Index position = interaction.queryRange().begin; position <= interaction.queryRange().end; ++position) {
                if (queryAccessibility.positionUnpairedProbability(position) + 1e-15 <
                    config_.output.minUnpairedProbability) probabilityRejected = true;
            }
            if (probabilityRejected) continue;

            const auto targetRange = interaction.targetRange();
            const auto queryRange = interaction.queryRange();
            const SiteKey key{targetRange.begin, targetRange.end, queryRange.begin, queryRange.end};
            // DP state weights contain initiation and interior transitions.
            // Site-specific exterior ends/dangles, accessibility and additive
            // energy complete every structure's Boltzmann exponent exactly once.
            const auto completeLogWeight = stateLogWeight -
                (interaction.energy.dangleLeft + interaction.energy.dangleRight +
                 interaction.energy.endLeft + interaction.energy.endRight +
                 interaction.energy.openingTarget + interaction.energy.openingQuery +
                 interaction.energy.additive) / energy_->rt();
            auto [iterator, inserted] = sites.try_emplace(key);
            if (inserted) {
                iterator->second.representative = std::move(interaction);
                iterator->second.logWeight = completeLogWeight;
            } else {
                iterator->second.logWeight = logAdd(iterator->second.logWeight, completeLogWeight);
                if (interactionLess(interaction, iterator->second.representative)) {
                    iterator->second.representative = std::move(interaction);
                }
            }
        }
        std::vector<PathState>{}.swap(cellStates);
    }
    std::vector<std::vector<PathState>>{}.swap(states);

    // The older single-site strategy keeps one interaction for each
    // left boundary. In antiparallel coordinates that boundary is identified
    // by (target begin, query end). Applying the same total ordering as final
    // output makes tied real-interaction fixtures deterministic.
    if (config_.mode == PredictionMode::heuristic &&
        config_.model == InteractionModel::singleSite) {
        std::map<std::pair<Index, Index>, SiteKey> bestSite;
        for (const auto& [siteKey, site] : sites) {
            const auto boundary = std::pair{siteKey.targetBegin, siteKey.queryEnd};
            const auto found = bestSite.find(boundary);
            if (found == bestSite.end() ||
                interactionLess(site.representative, sites.at(found->second).representative)) {
                bestSite.insert_or_assign(boundary, siteKey);
            }
        }
        std::map<SiteKey, SiteAccumulator> filtered;
        for (const auto& [boundary, siteKey] : bestSite) {
            static_cast<void>(boundary);
            filtered.emplace(siteKey, std::move(sites.at(siteKey)));
        }
        sites = std::move(filtered);
    }

    if (config_.mode == PredictionMode::heuristic &&
        config_.model == InteractionModel::seedExtension && seedRequired) {
        // Seed extension grows every seed in both directions independently.
        // All one-sided extensions are candidates, while two-sided candidates
        // use only the single best right extension of their seed. A shorter
        // right extension is intentionally not substituted when the best one
        // violates an outer interaction-length bound; this is the documented
        // speed/coverage tradeoff that distinguishes H from exact mode.
        struct RightExtension {
            Energy transitionEnergy{infinity};
            std::vector<BasePair> suffix;
        };
        std::map<std::string, RightExtension> bestRight;
        for (const auto& [siteKey, site] : sites) {
            static_cast<void>(siteKey);
            const auto& path = site.representative.pairs;
            for (const auto& seed : site.representative.seeds) {
                if (seed.firstPair != 0U || seed.lastPair >= path.size()) continue;
                const auto seedPath = std::span(path).subspan(
                    seed.firstPair, seed.lastPair - seed.firstPair + 1U);
                const auto key = pathSignature(seedPath);
                Energy extension{};
                for (Index index = seed.lastPair + 1U; index < path.size(); ++index) {
                    extension += energy_->transitionEnergy(
                        target, query, path[index - 1U], path[index]);
                }
                std::vector<BasePair> suffix(
                    path.begin() + static_cast<std::ptrdiff_t>(seed.lastPair), path.end());
                auto [found, inserted] = bestRight.try_emplace(
                    key, RightExtension{extension, suffix});
                if (!inserted &&
                    (extension < found->second.transitionEnergy - 1e-12 ||
                     (std::abs(extension - found->second.transitionEnergy) <= 1e-12 &&
                      pathSignatureLess(suffix, found->second.suffix)))) {
                    found->second = RightExtension{extension, std::move(suffix)};
                }
            }
        }

        std::map<SiteKey, SiteAccumulator> filtered;
        for (auto& [siteKey, site] : sites) {
            const auto& path = site.representative.pairs;
            const bool retain = std::ranges::any_of(
                site.representative.seeds, [&](const SeedMatch& seed) {
                    if (seed.lastPair >= path.size()) return false;
                    if (seed.firstPair == 0U || seed.lastPair + 1U == path.size()) return true;
                    const auto seedPath = std::span(path).subspan(
                        seed.firstPair, seed.lastPair - seed.firstPair + 1U);
                    const auto found = bestRight.find(pathSignature(seedPath));
                    if (found == bestRight.end()) return false;
                    const auto suffix = std::span(path).subspan(seed.lastPair);
                    return std::ranges::equal(suffix, found->second.suffix);
                });
            if (retain) filtered.emplace(siteKey, std::move(site));
        }
        sites = std::move(filtered);
    }

    result.ensembleSites.reserve(sites.size());
    const bool seededExtensionWithoutPartition =
        config_.model == InteractionModel::seedExtension && seedRequired;
    for (auto& [key, site] : sites) {
        static_cast<void>(key);
        site.representative.ensembleFreeEnergy =
            config_.model == InteractionModel::singleSite ||
                    config_.model == InteractionModel::seedExtension
                ? site.representative.energy.total()
                : -energy_->rt() * site.logWeight;
        if (config_.model == InteractionModel::ensemble) {
            const auto targetRange = site.representative.targetRange();
            const auto queryRange = site.representative.queryRange();
            site.representative.pairs = {{targetRange.begin, queryRange.end}};
            if (targetRange.begin != targetRange.end || queryRange.begin != queryRange.end) {
                site.representative.pairs.push_back({targetRange.end, queryRange.begin});
            }
            site.representative.seeds.clear();
            auto& energy = site.representative.energy;
            const auto centikcalSiteEnergy =
                std::trunc(site.representative.ensembleFreeEnergy * 100.0) / 100.0;
            energy.loops = centikcalSiteEnergy - energy.openingTarget -
                           energy.openingQuery - energy.additive - energy.initiation -
                           energy.dangleLeft - energy.dangleRight - energy.endLeft -
                           energy.endRight;
        }
        const auto reportedLogWeight =
            config_.model == InteractionModel::singleSite ||
                    config_.model == InteractionModel::seedExtension
            ? -site.representative.ensembleFreeEnergy / energy_->rt()
            : site.logWeight;
        if (!seededExtensionWithoutPartition) {
            result.logPartition = logAdd(result.logPartition, reportedLogWeight);
        }
        result.ensembleSites.push_back(std::move(site.representative));
    }
    sites.clear();
    if (std::isfinite(result.logPartition)) {
        result.ensembleFreeEnergy = -energy_->rt() * result.logPartition;
        for (auto& site : result.ensembleSites) {
            site.probability = std::exp(-site.ensembleFreeEnergy / energy_->rt() - result.logPartition);
        }
    }

    std::vector<const Interaction*> ranked;
    ranked.reserve(result.ensembleSites.size());
    for (const auto& interaction : result.ensembleSites) ranked.push_back(&interaction);
    orderForOutput(ranked, config_);
    if (ranked.empty() || config_.output.number == 0U) return result;
    const auto minimumEnergy = ranked.front()->energy.total();
    for (const auto* interaction : ranked) {
        if (interaction->energy.total() > minimumEnergy + config_.output.deltaEnergy + 1e-9) continue;
        if (overlapsSelected(*interaction, result.interactions, config_.output.overlap)) continue;
        result.interactions.push_back(*interaction);
        if (result.interactions.size() >= config_.output.number) break;
    }
    return result;
}

} // namespace intarnanew
