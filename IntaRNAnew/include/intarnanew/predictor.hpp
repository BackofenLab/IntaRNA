#pragma once

#include "intarnanew/accessibility.hpp"
#include "intarnanew/config.hpp"
#include "intarnanew/energy.hpp"
#include "intarnanew/sequence.hpp"
#include "intarnanew/types.hpp"

#include <expected>
#include <memory>
#include <span>
#include <vector>

namespace intarnanew {

struct PredictionResult {
    std::vector<Interaction> interactions;
    std::vector<Interaction> ensembleSites;
    double rt{};
    double logPartition{-infinity};
    Energy ensembleFreeEnergy{infinity};
    double targetLogPartition{-infinity};
    double queryLogPartition{-infinity};
    Energy targetEnsembleFreeEnergy{infinity};
    Energy queryEnsembleFreeEnergy{infinity};
};

class Predictor {
public:
    explicit Predictor(Config config);

    [[nodiscard]] auto predict(
        const Sequence& target,
        const Sequence& query,
        const AccessibilityProvider& targetAccessibility,
        const AccessibilityProvider& queryAccessibility) const -> PredictionResult;

    [[nodiscard]] auto predict(
        const Sequence& target,
        const Sequence& query,
        const AccessibilityProvider& targetAccessibility,
        const AccessibilityProvider& queryAccessibility,
        Interval targetDomain,
        Interval queryDomain) const -> PredictionResult;

    [[nodiscard]] auto rt() const noexcept -> double { return energy_->rt(); }

private:
    Config config_;
    std::unique_ptr<HybridEnergyModel> energy_;
};

[[nodiscard]] auto parseIntervals(const Sequence& sequence, const std::string& specification)
    -> std::expected<std::vector<Interval>, std::string>;

// Returns explicitly configured regions, or the whole sequence if no regions
// were specified. External coordinates are mapped with Sequence's skip-zero
// convention; malformed, empty, overlapping, and out-of-range intervals are
// rejected before prediction starts.
[[nodiscard]] auto configuredRegions(const Sequence& sequence, const SideConfig& config)
    -> std::expected<std::vector<Interval>, std::string>;

// Clean-room implementation of the documented accessible-region heuristic.
// Each range longer than regionLengthMax is split by removing the seed-length
// interval with the largest opening energy. Ties are resolved by the lowest
// internal index. Fragments shorter than seedLength are discarded.
[[nodiscard]] auto decomposeAccessibleRegions(
    std::span<const Interval> input,
    std::size_t regionLengthMax,
    std::size_t seedLength,
    const AccessibilityProvider& accessibility) -> std::vector<Interval>;

// Creates deterministic overlapping windows over a parent interval. The final
// window is clipped at the parent's end. Every interval of length <= overlap is
// contained in at least one window, including interactions crossing a step.
[[nodiscard]] auto decomposeWindows(
    Interval parent,
    std::size_t width,
    std::size_t overlap) -> std::vector<Interval>;

// Reduces independently predicted domains into one deterministic result.
// Interaction sites are deduplicated before the partition function is rebuilt
// and one global energy, delta-energy, count, and overlap-policy reduction is
// applied.
[[nodiscard]] auto reducePredictions(
    std::span<const PredictionResult> predictions,
    const Config& config) -> PredictionResult;

} // namespace intarnanew
