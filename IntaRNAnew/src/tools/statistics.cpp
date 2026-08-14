#include "intarnanew/tools/statistics.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <numbers>
#include <utility>

namespace intarnanew::tools {
namespace {

constexpr double eulerMascheroni = 0.577215664901532860606512090082402431;
constexpr double rootTwo = 1.414213562373095048801688724209698079;

[[nodiscard]] auto finiteObservations(const std::span<const double> observations)
    -> std::expected<void, std::string> {
    if (observations.size() < 2U) {
        return std::unexpected("at least two observations are required");
    }
    if (std::ranges::any_of(observations, [](const double value) {
            return !std::isfinite(value);
        })) {
        return std::unexpected("all observations must be finite");
    }
    const auto [minimum, maximum] = std::ranges::minmax(observations);
    if (minimum == maximum) {
        return std::unexpected("distribution fitting requires non-constant observations");
    }
    return {};
}

struct Moments {
    double mean{};
    double standardDeviation{};
};

[[nodiscard]] auto moments(const std::span<const double> values) noexcept -> Moments {
    // Welford's recurrence avoids catastrophic cancellation for shifted data.
    double mean{};
    double sumSquares{};
    std::size_t count{};
    for (const double value : values) {
        ++count;
        const double delta = value - mean;
        mean += delta / static_cast<double>(count);
        sumSquares += delta * (value - mean);
    }
    return {mean, std::sqrt(sumSquares / static_cast<double>(count))};
}

[[nodiscard]] auto gaussianNll(
    const std::span<const double> observations,
    const double location,
    const double scale) noexcept -> double {
    if (!(scale > 0.0) || !std::isfinite(location) || !std::isfinite(scale)) {
        return std::numeric_limits<double>::infinity();
    }
    long double squares{};
    for (const double value : observations) {
        const long double normalized =
            (static_cast<long double>(value) - location) / scale;
        squares += normalized * normalized;
    }
    return static_cast<double>(
        static_cast<long double>(observations.size()) *
            std::log(scale * std::sqrt(2.0 * std::numbers::pi)) +
        0.5L * squares);
}

[[nodiscard]] auto gumbelNll(
    const std::span<const double> observations,
    const double location,
    const double logScale) noexcept -> double {
    const double scale = std::exp(logScale);
    if (!(scale > 0.0) || !std::isfinite(location) || !std::isfinite(scale)) {
        return std::numeric_limits<double>::infinity();
    }
    long double sum{};
    for (const double value : observations) {
        const double normalized = (value - location) / scale;
        if (!std::isfinite(normalized) || normalized < -700.0) {
            return std::numeric_limits<double>::infinity();
        }
        sum += normalized + std::exp(-normalized);
    }
    return static_cast<double>(
        static_cast<long double>(observations.size()) * logScale + sum);
}

[[nodiscard]] auto gevNll(
    const std::span<const double> observations,
    const double location,
    const double logScale,
    const double shape) noexcept -> double {
    const double scale = std::exp(logScale);
    if (!(scale > 0.0) || !std::isfinite(location) || !std::isfinite(scale) ||
        !std::isfinite(shape) || std::abs(shape) > 1.0) {
        return std::numeric_limits<double>::infinity();
    }
    if (std::abs(shape) < 1.0e-7) {
        return gumbelNll(observations, location, logScale);
    }
    long double sum{};
    for (const double value : observations) {
        const double support = 1.0 + shape * ((value - location) / scale);
        if (!(support > 0.0) || !std::isfinite(support)) {
            return std::numeric_limits<double>::infinity();
        }
        const double logSupport = std::log(support);
        const double exponential = std::exp(-logSupport / shape);
        if (!std::isfinite(exponential)) {
            return std::numeric_limits<double>::infinity();
        }
        sum += (1.0 + 1.0 / shape) * logSupport + exponential;
    }
    return static_cast<double>(
        static_cast<long double>(observations.size()) * logScale + sum);
}

template <std::size_t Dimensions, typename Objective>
struct SearchResult {
    std::array<double, Dimensions> point{};
    double value{std::numeric_limits<double>::infinity()};
    std::size_t iterations{};
    bool converged{};
};

template <std::size_t Dimensions, typename Objective>
[[nodiscard]] auto nelderMead(
    const std::array<double, Dimensions>& initial,
    const std::array<double, Dimensions>& step,
    Objective objective) -> SearchResult<Dimensions, Objective> {
    struct Vertex {
        std::array<double, Dimensions> point{};
        double value{};
    };
    std::array<Vertex, Dimensions + 1U> simplex{};
    simplex[0].point = initial;
    simplex[0].value = objective(initial);
    for (std::size_t index = 0; index < Dimensions; ++index) {
        simplex[index + 1U].point = initial;
        simplex[index + 1U].point[index] += step[index];
        simplex[index + 1U].value = objective(simplex[index + 1U].point);
    }

    constexpr std::size_t maximumIterations = 4000U;
    constexpr double relativeTolerance = 1.0e-11;
    std::size_t iteration{};
    bool converged{};
    for (; iteration < maximumIterations; ++iteration) {
        std::ranges::sort(simplex, {}, &Vertex::value);
        double coordinateSpread{};
        for (std::size_t vertex = 1U; vertex < simplex.size(); ++vertex) {
            for (std::size_t dimension = 0; dimension < Dimensions; ++dimension) {
                coordinateSpread = std::max(
                    coordinateSpread,
                    std::abs(simplex[vertex].point[dimension] - simplex[0].point[dimension]) /
                        (1.0 + std::abs(simplex[0].point[dimension])));
            }
        }
        const double valueSpread =
            std::abs(simplex.back().value - simplex.front().value) /
            (1.0 + std::abs(simplex.front().value));
        if (coordinateSpread <= relativeTolerance && valueSpread <= relativeTolerance) {
            converged = true;
            break;
        }

        std::array<double, Dimensions> centroid{};
        for (std::size_t vertex = 0U; vertex < Dimensions; ++vertex) {
            for (std::size_t dimension = 0U; dimension < Dimensions; ++dimension) {
                centroid[dimension] += simplex[vertex].point[dimension] /
                                       static_cast<double>(Dimensions);
            }
        }
        const auto moved = [&](const Vertex& from, const double factor) {
            Vertex result;
            for (std::size_t dimension = 0U; dimension < Dimensions; ++dimension) {
                result.point[dimension] = centroid[dimension] +
                    factor * (centroid[dimension] - from.point[dimension]);
            }
            result.value = objective(result.point);
            return result;
        };

        auto reflected = moved(simplex.back(), 1.0);
        if (reflected.value < simplex.front().value) {
            auto expanded = moved(simplex.back(), 2.0);
            simplex.back() = expanded.value < reflected.value ? expanded : reflected;
        } else if (reflected.value < simplex[Dimensions - 1U].value) {
            simplex.back() = reflected;
        } else {
            Vertex contracted;
            if (reflected.value < simplex.back().value) {
                for (std::size_t dimension = 0U; dimension < Dimensions; ++dimension) {
                    contracted.point[dimension] = centroid[dimension] +
                        0.5 * (reflected.point[dimension] - centroid[dimension]);
                }
            } else {
                for (std::size_t dimension = 0U; dimension < Dimensions; ++dimension) {
                    contracted.point[dimension] = centroid[dimension] +
                        0.5 * (simplex.back().point[dimension] - centroid[dimension]);
                }
            }
            contracted.value = objective(contracted.point);
            if (contracted.value < simplex.back().value) {
                simplex.back() = contracted;
            } else {
                for (std::size_t vertex = 1U; vertex < simplex.size(); ++vertex) {
                    for (std::size_t dimension = 0U; dimension < Dimensions; ++dimension) {
                        simplex[vertex].point[dimension] = simplex.front().point[dimension] +
                            0.5 * (simplex[vertex].point[dimension] -
                                   simplex.front().point[dimension]);
                    }
                    simplex[vertex].value = objective(simplex[vertex].point);
                }
            }
        }
    }
    std::ranges::sort(simplex, {}, &Vertex::value);
    return {simplex.front().point, simplex.front().value, iteration, converged};
}

[[nodiscard]] auto validFit(const DistributionFit& fit) -> std::expected<void, std::string> {
    if (!(fit.scale > 0.0) || !std::isfinite(fit.location) || !std::isfinite(fit.scale) ||
        !std::isfinite(fit.shape)) {
        return std::unexpected("distribution parameters are not finite and valid");
    }
    if (fit.kind == DistributionKind::gev && std::abs(fit.shape) > 1.0) {
        return std::unexpected("GEV shape must be in [-1,1]");
    }
    return {};
}

[[nodiscard]] auto sortedIndices(const std::span<const double> values)
    -> std::vector<std::size_t> {
    std::vector<std::size_t> indices(values.size());
    std::iota(indices.begin(), indices.end(), 0U);
    std::ranges::stable_sort(indices, [&](const auto left, const auto right) {
        return values[left] < values[right];
    });
    return indices;
}

} // namespace

auto parseDistribution(const std::string_view name)
    -> std::expected<DistributionKind, std::string> {
    if (name == "gauss" || name == "gaussian" || name == "normal") {
        return DistributionKind::gaussian;
    }
    if (name == "gumbel") {
        return DistributionKind::gumbel;
    }
    if (name == "gev") {
        return DistributionKind::gev;
    }
    return std::unexpected("unknown distribution '" + std::string{name} + "'");
}

auto distributionName(const DistributionKind kind) noexcept -> std::string_view {
    switch (kind) {
        case DistributionKind::gaussian: return "gaussian";
        case DistributionKind::gumbel: return "gumbel";
        case DistributionKind::gev: return "gev";
    }
    return "unknown";
}

auto fitDistribution(
    const std::span<const double> observations,
    const DistributionKind kind) -> std::expected<DistributionFit, std::string> {
    auto validation = finiteObservations(observations);
    if (!validation) {
        return std::unexpected(validation.error());
    }
    if (kind == DistributionKind::gev && observations.size() < 3U) {
        return std::unexpected("GEV fitting requires at least three observations");
    }
    const auto sample = moments(observations);
    if (kind == DistributionKind::gaussian) {
        return DistributionFit{
            kind,
            sample.mean,
            sample.standardDeviation,
            0.0,
            gaussianNll(observations, sample.mean, sample.standardDeviation),
            observations.size(),
            true,
            0U};
    }

    const double initialScale = sample.standardDeviation * std::sqrt(6.0) / std::numbers::pi;
    const double initialLocation = sample.mean - eulerMascheroni * initialScale;
    if (kind == DistributionKind::gumbel) {
        const std::array<double, 2U> initial{initialLocation, std::log(initialScale)};
        const std::array<double, 2U> step{
            std::max(sample.standardDeviation * 0.1, 1.0e-3), 0.1};
        auto result = nelderMead<2U>(initial, step, [&](const auto& point) {
            return gumbelNll(observations, point[0], point[1]);
        });
        if (!std::isfinite(result.value)) {
            return std::unexpected("Gumbel optimization did not find valid parameters");
        }
        return DistributionFit{
            kind,
            result.point[0],
            std::exp(result.point[1]),
            0.0,
            result.value,
            observations.size(),
            result.converged,
            result.iterations};
    }

    // Multiple deterministic starts reduce local-minimum and support-boundary
    // sensitivity without introducing a random source.
    constexpr std::array<double, 5U> starts{-0.25, -0.1, 0.0, 0.1, 0.25};
    SearchResult<3U, decltype([&](const std::array<double, 3U>&) { return 0.0; })> best;
    // SearchResult's Objective template argument is only a tag and is not stored.
    best.value = std::numeric_limits<double>::infinity();
    for (const double shape : starts) {
        const std::array<double, 3U> initial{
            initialLocation, std::log(initialScale), shape};
        const std::array<double, 3U> step{
            std::max(sample.standardDeviation * 0.1, 1.0e-3), 0.1, 0.05};
        auto result = nelderMead<3U>(initial, step, [&](const auto& point) {
            return gevNll(observations, point[0], point[1], point[2]);
        });
        if (result.value < best.value) {
            best.point = result.point;
            best.value = result.value;
            best.iterations = result.iterations;
            best.converged = result.converged;
        }
    }
    if (!std::isfinite(best.value)) {
        return std::unexpected("GEV optimization did not find valid parameters");
    }
    return DistributionFit{
        kind,
        best.point[0],
        std::exp(best.point[1]),
        best.point[2],
        best.value,
        observations.size(),
        best.converged,
        best.iterations};
}

auto cumulativeProbability(const double value, const DistributionFit& fit)
    -> std::expected<double, std::string> {
    if (!std::isfinite(value)) {
        return std::unexpected("probability evaluation value must be finite");
    }
    auto validation = validFit(fit);
    if (!validation) {
        return std::unexpected(validation.error());
    }
    double cumulative{};
    if (fit.kind == DistributionKind::gaussian) {
        cumulative = 0.5 * std::erfc(-(value - fit.location) / (fit.scale * rootTwo));
    } else if (fit.kind == DistributionKind::gumbel || std::abs(fit.shape) < 1.0e-7) {
        cumulative = std::exp(-std::exp(-(value - fit.location) / fit.scale));
    } else {
        const double support = 1.0 + fit.shape * ((value - fit.location) / fit.scale);
        if (!(support > 0.0)) {
            cumulative = fit.shape > 0.0 ? 0.0 : 1.0;
        } else {
            cumulative = std::exp(-std::pow(support, -1.0 / fit.shape));
        }
    }
    return std::clamp(cumulative, 0.0, 1.0);
}

auto tailProbability(
    const double value,
    const DistributionFit& fit,
    const ProbabilityTail tail) -> std::expected<double, std::string> {
    auto cumulative = cumulativeProbability(value, fit);
    if (!cumulative) {
        return std::unexpected(cumulative.error());
    }
    if (tail == ProbabilityTail::lower) {
        return *cumulative;
    }
    // Avoid cancellation in 1-CDF for values far into the upper tail.
    if (fit.kind == DistributionKind::gaussian) {
        return 0.5 * std::erfc((value - fit.location) / (fit.scale * rootTwo));
    }
    if (fit.kind == DistributionKind::gumbel || std::abs(fit.shape) < 1.0e-7) {
        const double exponent = std::exp(-(value - fit.location) / fit.scale);
        return std::clamp(-std::expm1(-exponent), 0.0, 1.0);
    }
    const double support = 1.0 + fit.shape * ((value - fit.location) / fit.scale);
    if (!(support > 0.0)) {
        return fit.shape > 0.0 ? 1.0 : 0.0;
    }
    const double exponent = std::pow(support, -1.0 / fit.shape);
    return std::clamp(-std::expm1(-exponent), 0.0, 1.0);
}

auto interactionEnergyPValue(const double energy, const DistributionFit& fit)
    -> std::expected<double, std::string> {
    return tailProbability(energy, fit, ProbabilityTail::lower);
}

auto empiricalInteractionPValue(
    const double observed,
    const std::span<const double> randomScores) -> std::expected<double, std::string> {
    if (!std::isfinite(observed)) {
        return std::unexpected("observed interaction energy must be finite");
    }
    if (randomScores.empty()) {
        return std::unexpected("at least one random score is required");
    }
    std::size_t atLeastAsGood{};
    for (const double score : randomScores) {
        if (!std::isfinite(score)) {
            return std::unexpected("random scores must be finite");
        }
        atLeastAsGood += score <= observed;
    }
    return static_cast<double>(atLeastAsGood + 1U) /
           static_cast<double>(randomScores.size() + 1U);
}

auto parseAdjustment(const std::string_view name)
    -> std::expected<AdjustmentMethod, std::string> {
    if (name == "none") return AdjustmentMethod::none;
    if (name == "bonferroni" || name == "bonf") return AdjustmentMethod::bonferroni;
    if (name == "holm") return AdjustmentMethod::holm;
    if (name == "hochberg") return AdjustmentMethod::hochberg;
    if (name == "bh" || name == "fdr" || name == "benjamini-hochberg") {
        return AdjustmentMethod::benjaminiHochberg;
    }
    if (name == "by" || name == "benjamini-yekutieli") {
        return AdjustmentMethod::benjaminiYekutieli;
    }
    return std::unexpected("unknown p-value adjustment method '" + std::string{name} + "'");
}

auto adjustmentName(const AdjustmentMethod method) noexcept -> std::string_view {
    switch (method) {
        case AdjustmentMethod::none: return "none";
        case AdjustmentMethod::bonferroni: return "bonferroni";
        case AdjustmentMethod::holm: return "holm";
        case AdjustmentMethod::hochberg: return "hochberg";
        case AdjustmentMethod::benjaminiHochberg: return "benjamini-hochberg";
        case AdjustmentMethod::benjaminiYekutieli: return "benjamini-yekutieli";
    }
    return "unknown";
}

auto adjustPValues(
    const std::span<const double> pValues,
    const AdjustmentMethod method) -> std::expected<std::vector<double>, std::string> {
    if (pValues.empty()) {
        return std::unexpected("at least one p-value is required");
    }
    for (const double value : pValues) {
        if (!std::isfinite(value) || value < 0.0 || value > 1.0) {
            return std::unexpected("all p-values must be finite and in [0,1]");
        }
    }
    std::vector<double> adjusted(pValues.begin(), pValues.end());
    const double count = static_cast<double>(pValues.size());
    if (method == AdjustmentMethod::none) {
        return adjusted;
    }
    if (method == AdjustmentMethod::bonferroni) {
        for (auto& value : adjusted) {
            value = std::min(1.0, value * count);
        }
        return adjusted;
    }

    const auto order = sortedIndices(pValues);
    if (method == AdjustmentMethod::holm) {
        double running{};
        for (std::size_t rank = 0; rank < order.size(); ++rank) {
            running = std::max(
                running,
                (count - static_cast<double>(rank)) * pValues[order[rank]]);
            adjusted[order[rank]] = std::min(1.0, running);
        }
        return adjusted;
    }

    double multiplier = count;
    if (method == AdjustmentMethod::benjaminiYekutieli) {
        double harmonic{};
        for (std::size_t index = 1U; index <= pValues.size(); ++index) {
            harmonic += 1.0 / static_cast<double>(index);
        }
        multiplier *= harmonic;
    }
    double running = 1.0;
    for (std::size_t reverse = order.size(); reverse-- > 0U;) {
        const double rank = static_cast<double>(reverse + 1U);
        double candidate{};
        if (method == AdjustmentMethod::hochberg) {
            candidate = (count - static_cast<double>(reverse)) * pValues[order[reverse]];
        } else {
            candidate = multiplier * pValues[order[reverse]] / rank;
        }
        running = std::min(running, candidate);
        adjusted[order[reverse]] = std::clamp(running, 0.0, 1.0);
    }
    return adjusted;
}

} // namespace intarnanew::tools
