#pragma once

#include <cstddef>
#include <expected>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace intarnanew::tools {

enum class DistributionKind { gaussian, gumbel, gev };
enum class ProbabilityTail { lower, upper };
enum class AdjustmentMethod { none, bonferroni, holm, hochberg, benjaminiHochberg, benjaminiYekutieli };

struct DistributionFit {
    DistributionKind kind{DistributionKind::gev};
    double location{};
    double scale{1.0};
    // The GEV shape xi. It is exactly zero for Gaussian and Gumbel fits.
    double shape{};
    double negativeLogLikelihood{};
    std::size_t observations{};
    bool converged{};
    std::size_t iterations{};
};

[[nodiscard]] auto parseDistribution(std::string_view name)
    -> std::expected<DistributionKind, std::string>;
[[nodiscard]] auto distributionName(DistributionKind kind) noexcept -> std::string_view;

// Fits finite observations with deterministic maximum-likelihood estimation.
// Gaussian uses the population MLE. Gumbel and GEV use a bounded, deterministic
// Nelder-Mead search. GEV shape is constrained to [-1,1].
[[nodiscard]] auto fitDistribution(
    std::span<const double> observations,
    DistributionKind kind) -> std::expected<DistributionFit, std::string>;

[[nodiscard]] auto cumulativeProbability(double value, const DistributionFit& fit)
    -> std::expected<double, std::string>;
[[nodiscard]] auto tailProbability(
    double value,
    const DistributionFit& fit,
    ProbabilityTail tail = ProbabilityTail::lower) -> std::expected<double, std::string>;

// For interaction energies, a smaller value is better. This convenience
// function therefore returns the fitted lower-tail probability P(X <= energy).
[[nodiscard]] auto interactionEnergyPValue(double energy, const DistributionFit& fit)
    -> std::expected<double, std::string>;

// Exact empirical lower-tail estimate with the standard plus-one correction:
// (#{sample <= observed} + 1) / (sample size + 1).
[[nodiscard]] auto empiricalInteractionPValue(
    double observed,
    std::span<const double> randomScores) -> std::expected<double, std::string>;

[[nodiscard]] auto parseAdjustment(std::string_view name)
    -> std::expected<AdjustmentMethod, std::string>;
[[nodiscard]] auto adjustmentName(AdjustmentMethod method) noexcept -> std::string_view;

// Adjusts p-values in their original order. Ties are ordered by original index,
// making every method bit-for-bit deterministic for a given standard library.
[[nodiscard]] auto adjustPValues(
    std::span<const double> pValues,
    AdjustmentMethod method) -> std::expected<std::vector<double>, std::string>;

} // namespace intarnanew::tools
