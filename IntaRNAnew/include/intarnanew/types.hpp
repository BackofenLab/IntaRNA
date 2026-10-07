#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace intarnanew {

using Index = std::size_t;
using Energy = double;

inline constexpr Energy infinity = std::numeric_limits<Energy>::infinity();
// ViennaRNA/IntaRNA thermodynamic convention in kcal mol^-1 K^-1.
inline constexpr Energy gasConstantKcal = 0.00198717;

struct Interval {
    Index begin{};
    Index end{};

    [[nodiscard]] constexpr auto size() const noexcept -> Index {
        return end >= begin ? end - begin + 1U : 0U;
    }

    [[nodiscard]] constexpr auto contains(const Index value) const noexcept -> bool {
        return begin <= value && value <= end;
    }

    [[nodiscard]] constexpr auto overlaps(const Interval& other) const noexcept -> bool {
        return begin <= other.end && other.begin <= end;
    }

    friend constexpr auto operator==(const Interval&, const Interval&) -> bool = default;
};

struct BasePair {
    Index target{};
    Index query{};

    friend constexpr auto operator==(const BasePair&, const BasePair&) -> bool = default;
};

struct EnergyBreakdown {
    Energy openingTarget{};
    Energy openingQuery{};
    Energy initiation{};
    Energy loops{};
    Energy dangleLeft{};
    Energy dangleRight{};
    Energy endLeft{};
    Energy endRight{};
    Energy additive{};

    [[nodiscard]] constexpr auto hybrid() const noexcept -> Energy {
        return initiation + loops + dangleLeft + dangleRight + endLeft + endRight + additive;
    }

    [[nodiscard]] constexpr auto total() const noexcept -> Energy {
        return openingTarget + openingQuery + hybrid();
    }
};

struct SeedMatch {
    Index firstPair{};
    Index lastPair{};
    Energy energy{};
    Energy openingTarget{};
    Energy openingQuery{};
    double unpairedTarget{1.0};
    double unpairedQuery{1.0};

    friend auto operator==(const SeedMatch&, const SeedMatch&) -> bool = default;
};

[[nodiscard]] inline auto seedMatchLess(
    const SeedMatch& left,
    const SeedMatch& right) noexcept -> bool {
    if (left.energy != right.energy) return left.energy < right.energy;
    if (left.firstPair != right.firstPair) return left.firstPair < right.firstPair;
    return left.lastPair < right.lastPair;
}

struct Interaction {
    std::string targetId;
    std::string queryId;
    std::vector<BasePair> pairs;
    EnergyBreakdown energy;
    // Preserve the accessibility provider's native probabilities. They must
    // not be reconstructed from ED with the interaction energy model's RT:
    // the pedagogical base-pair model deliberately uses RT=1 while folding
    // accessibility remains a physical-temperature ensemble.
    double unpairedTarget{1.0};
    double unpairedQuery{1.0};
    std::vector<SeedMatch> seeds;
    Energy ensembleFreeEnergy{infinity};
    Energy probability{};

    [[nodiscard]] auto bestSeed() const noexcept -> const SeedMatch* {
        const auto best = std::ranges::min_element(seeds, seedMatchLess);
        return best == seeds.end() ? nullptr : &*best;
    }

    [[nodiscard]] auto targetRange() const -> Interval {
        if (pairs.empty()) {
            return {};
        }
        const auto [low, high] = std::ranges::minmax_element(
            pairs, {}, &BasePair::target);
        return {low->target, high->target};
    }

    [[nodiscard]] auto queryRange() const -> Interval {
        if (pairs.empty()) {
            return {};
        }
        const auto [low, high] = std::ranges::minmax_element(
            pairs, {}, &BasePair::query);
        return {low->query, high->query};
    }
};

} // namespace intarnanew
