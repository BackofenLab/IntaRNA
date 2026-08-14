#include "intarnanew/helix_blocks.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>

namespace intarnanew {
namespace {

[[nodiscard]] auto gapSize(const BasePair left, const BasePair right) -> Index {
    if (left.target >= right.target || left.query <= right.query) {
        throw std::invalid_argument("helix block path is not monotone antiparallel");
    }
    return right.target - left.target - 1U + left.query - right.query - 1U;
}

[[nodiscard]] auto blockConstraintEnergy(
    const Sequence& target,
    const Sequence& query,
    const std::span<const BasePair> block,
    const HelixConfig& config,
    const HybridEnergyModel& energy,
    const AccessibilityProvider& targetAccessibility,
    const AccessibilityProvider& queryAccessibility) -> std::optional<Energy> {
    auto breakdown = energy.evaluate(target, query, block);
    const auto targetRange = Interval{block.front().target, block.back().target};
    const auto queryRange = Interval{block.back().query, block.front().query};
    if (targetAccessibility.unpairedProbability(targetRange) + 1e-15 <
            config.minUnpairedProbability ||
        queryAccessibility.unpairedProbability(queryRange) + 1e-15 <
            config.minUnpairedProbability) {
        return std::nullopt;
    }
    Energy constraintEnergy = breakdown.loops;
    if (config.useFullEnergy) {
        breakdown.openingTarget = targetAccessibility.openingEnergy(targetRange);
        breakdown.openingQuery = queryAccessibility.openingEnergy(queryRange);
        constraintEnergy = breakdown.total();
    }
    // The documented bound is exclusive: E(helix) < helixMaxE.
    if (!std::isfinite(constraintEnergy) ||
        constraintEnergy >= config.maxEnergy - 1e-12) {
        return std::nullopt;
    }
    return constraintEnergy;
}

} // namespace

auto decomposeHelixBlocks(
    const Sequence& target,
    const Sequence& query,
    const std::span<const BasePair> path,
    const HelixConfig& config,
    const HybridEnergyModel& energy,
    const AccessibilityProvider& targetAccessibility,
    const AccessibilityProvider& queryAccessibility) -> std::vector<HelixBlock> {
    if (path.empty() || config.minBasePairs == 0U ||
        config.minBasePairs > config.maxBasePairs) {
        return {};
    }

    struct Prefix {
        Energy score{infinity};
        Index blocks{std::numeric_limits<Index>::max()};
        std::optional<Index> previous;
        Energy lastConstraintEnergy{};
    };
    const auto decomposePrefix = [&](const Index pathEnd) {
        std::vector<Prefix> best(pathEnd + 1U);
        best.front().score = 0.0;
        best.front().blocks = 0U;
        for (Index end = 1U; end <= pathEnd; ++end) {
            const auto maximum = std::min(config.maxBasePairs, end);
            for (Index count = config.minBasePairs; count <= maximum; ++count) {
                const Index begin = end - count;
                bool withinBlock = true;
                for (Index index = begin + 1U; index < end; ++index) {
                    if (gapSize(path[index - 1U], path[index]) > config.maxInternalLoop) {
                        withinBlock = false;
                        break;
                    }
                }
                if (!withinBlock) continue;
                if (begin != 0U) {
                    if (!best[begin].previous) continue;
                    if (gapSize(path[begin - 1U], path[begin]) <= config.maxInternalLoop) continue;
                }
                const auto constraintEnergy = blockConstraintEnergy(
                    target, query, path.subspan(begin, count), config, energy,
                    targetAccessibility, queryAccessibility);
                if (!constraintEnergy) continue;

                const auto score = best[begin].score + *constraintEnergy;
                const auto blocks = best[begin].blocks + 1U;
                const bool improves = score < best[end].score - 1e-12 ||
                    (std::abs(score - best[end].score) <= 1e-12 && blocks < best[end].blocks);
                if (improves) {
                    best[end] = Prefix{score, blocks, begin, *constraintEnergy};
                }
            }
        }
        return best;
    };

    auto best = decomposePrefix(path.size());
    Index tracebackEnd = path.size();
    if (!best.back().previous && path.size() > 1U &&
        gapSize(path[path.size() - 2U], path.back()) > config.maxInternalLoop) {
        best = decomposePrefix(path.size() - 1U);
        tracebackEnd = path.size() - 1U;
    }
    if (!best.back().previous) return {};

    std::vector<HelixBlock> result;
    for (Index end = tracebackEnd; end != 0U;) {
        const auto& prefix = best[end];
        if (!prefix.previous) return {};
        const auto begin = *prefix.previous;
        result.push_back({begin, end - 1U, prefix.lastConstraintEnergy});
        end = begin;
    }
    std::ranges::reverse(result);
    return result;
}

} // namespace intarnanew
