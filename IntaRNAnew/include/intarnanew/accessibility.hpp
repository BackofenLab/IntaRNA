#pragma once

#include "intarnanew/config.hpp"
#include "intarnanew/folding.hpp"
#include "intarnanew/sequence.hpp"

#include <expected>
#include <memory>
#include <mutex>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace intarnanew {

class AccessibilityProvider {
public:
    virtual ~AccessibilityProvider() = default;

    [[nodiscard]] virtual auto openingEnergy(Interval interval) const -> Energy = 0;
    [[nodiscard]] virtual auto unpairedProbability(Interval interval) const -> double = 0;
    [[nodiscard]] virtual auto positionUnpairedProbability(Index position) const -> double = 0;
    [[nodiscard]] virtual auto blocked(Index position) const -> bool = 0;
    [[nodiscard]] virtual auto ensembleLogPartition() const noexcept -> std::optional<double> {
        return std::nullopt;
    }
    [[nodiscard]] virtual auto ensembleFreeEnergy() const noexcept -> std::optional<Energy> {
        return std::nullopt;
    }
};

class TurnerAccessibility final : public AccessibilityProvider {
public:
    // Interval probabilities use either the configured global fold or the
    // average over every full local folding window containing the interval.
    // The ensemble accessors always describe a separate, unrestricted
    // whole-sequence monomer ensemble, independent of accessibility windows.
    static constexpr std::string_view modelName{"turner-nearest-neighbor-v1"};

    TurnerAccessibility(
        const Sequence& sequence,
        const SideConfig& config,
        const Config& globalConfig);

    [[nodiscard]] auto openingEnergy(Interval interval) const -> Energy override;
    [[nodiscard]] auto unpairedProbability(Interval interval) const -> double override;
    [[nodiscard]] auto positionUnpairedProbability(Index position) const -> double override;
    [[nodiscard]] auto blocked(Index position) const -> bool override;
    [[nodiscard]] auto ensembleLogPartition() const noexcept -> std::optional<double> override;
    [[nodiscard]] auto ensembleFreeEnergy() const noexcept -> std::optional<Energy> override;

private:
    struct Window {
        Index begin{};
        Index end{};
        std::shared_ptr<FoldingEnsemble> ensemble;
    };

    Index length_{};
    Index stride_{};
    Index windowLength_{};
    double rt_{};
    double probabilityExponent_{1.0};
    std::vector<bool> blocked_;
    std::shared_ptr<FoldingEnsemble> summary_;
    std::vector<Window> windows_;
    mutable std::vector<double> intervalProbabilities_;
    mutable std::mutex cacheMutex_;
};

class DisabledAccessibility final : public AccessibilityProvider {
public:
    explicit DisabledAccessibility(const Sequence& sequence, std::string_view constraint = {});
    DisabledAccessibility(const Sequence& sequence, const SideConfig& side, const Config& globalConfig);

    [[nodiscard]] auto openingEnergy(Interval interval) const -> Energy override;
    [[nodiscard]] auto unpairedProbability(Interval interval) const -> double override;
    [[nodiscard]] auto positionUnpairedProbability(Index position) const -> double override;
    [[nodiscard]] auto blocked(Index position) const -> bool override;
    [[nodiscard]] auto unconstrained(Interval interval) const noexcept -> bool;
    [[nodiscard]] auto ensembleLogPartition() const noexcept -> std::optional<double> override;
    [[nodiscard]] auto ensembleFreeEnergy() const noexcept -> std::optional<Energy> override;

private:
    Index length_{};
    std::vector<bool> blocked_;
    std::unique_ptr<FoldingEnsemble> ensemble_;
};

class NativeAccessibility final : public AccessibilityProvider {
public:
    // Exact for the documented native model: a pseudoknot-free noncrossing
    // matching ensemble with additive pair energies. This is deliberately not
    // a Turner nearest-neighbour, RNAplfold-window, or ViennaRNA-compatible
    // folding model. accessibilitySpan caps pair distance; accessibilityWindow
    // does not alter this global ensemble. One inside/outside calculation and
    // two log-semiring contractions generate every joint interval probability
    // in O(n^3) time/O(n^2) memory; each query is then O(1).
    static constexpr std::string_view modelName{"native-noncrossing-pair-v1"};

    NativeAccessibility(const Sequence& sequence, const SideConfig& config, double temperatureCelsius);

    [[nodiscard]] auto openingEnergy(Interval interval) const -> Energy override;
    [[nodiscard]] auto unpairedProbability(Interval interval) const -> double override;
    [[nodiscard]] auto positionUnpairedProbability(Index position) const -> double override;
    [[nodiscard]] auto blocked(Index position) const -> bool override;

private:
    [[nodiscard]] auto cacheOffset(Interval interval) const -> Index;

    Index length_{};
    Index stride_{};
    Index maxPairSpan_{};
    double rt_{};
    double logPartition_{};
    std::string sequence_;
    std::vector<char> constraints_;
    std::vector<bool> blocked_;
    std::vector<double> intervalProbabilities_;
};

class TableAccessibility final : public AccessibilityProvider {
public:
    TableAccessibility(
        const Sequence& sequence,
        const SideConfig& config,
        double temperatureCelsius,
        bool tableContainsProbabilities);
    TableAccessibility(
        const Sequence& sequence,
        const SideConfig& config,
        const Config& globalConfig,
        bool tableContainsProbabilities);

    [[nodiscard]] auto openingEnergy(Interval interval) const -> Energy override;
    [[nodiscard]] auto unpairedProbability(Interval interval) const -> double override;
    [[nodiscard]] auto positionUnpairedProbability(Index position) const -> double override;
    [[nodiscard]] auto blocked(Index position) const -> bool override;
    [[nodiscard]] auto ensembleLogPartition() const noexcept -> std::optional<double> override;
    [[nodiscard]] auto ensembleFreeEnergy() const noexcept -> std::optional<Energy> override;

private:
    [[nodiscard]] auto offset(Interval interval) const -> Index;

    Index length_{};
    Index stride_{};
    double rt_{};
    std::vector<double> probabilities_;
    std::vector<bool> blocked_;
    std::unique_ptr<FoldingEnsemble> ensemble_;
};

[[nodiscard]] auto makeAccessibility(
    const Sequence& sequence,
    const SideConfig& config,
    double temperatureCelsius) -> std::expected<std::unique_ptr<AccessibilityProvider>, std::string>;

[[nodiscard]] auto makeAccessibility(
    const Sequence& sequence,
    const SideConfig& config,
    const Config& globalConfig) -> std::expected<std::unique_ptr<AccessibilityProvider>, std::string>;

} // namespace intarnanew
