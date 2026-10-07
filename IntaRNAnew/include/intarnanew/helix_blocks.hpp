#pragma once

#include "intarnanew/accessibility.hpp"
#include "intarnanew/config.hpp"
#include "intarnanew/energy.hpp"
#include "intarnanew/sequence.hpp"
#include "intarnanew/types.hpp"

#include <span>
#include <vector>

namespace intarnanew {

struct HelixBlock {
    Index firstPair{};
    Index lastPair{};
    Energy constraintEnergy{};

    [[nodiscard]] constexpr auto basePairCount() const noexcept -> Index {
        return lastPair - firstPair + 1U;
    }

    friend constexpr auto operator==(const HelixBlock&, const HelixBlock&) -> bool = default;
};

// Partitions an antiparallel interaction path into admissible helix blocks.
// A block contains between helixMinBP and helixMaxBP pairs. Consecutive pairs
// inside a block may enclose at most helixMaxIL unpaired bases in total. Two
// blocks must be separated by a larger loop, matching the published numeric
// contract. The right-most initiation pair may terminate a composition after
// at least one admissible block; it is not itself a helix block. Every returned
// block passes the strict helixMaxE threshold and the per-sequence helixMinPu
// threshold. The returned partition is deterministic.
[[nodiscard]] auto decomposeHelixBlocks(
    const Sequence& target,
    const Sequence& query,
    std::span<const BasePair> path,
    const HelixConfig& config,
    const HybridEnergyModel& energy,
    const AccessibilityProvider& targetAccessibility,
    const AccessibilityProvider& queryAccessibility) -> std::vector<HelixBlock>;

} // namespace intarnanew
