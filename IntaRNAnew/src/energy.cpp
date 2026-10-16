#include "intarnanew/energy.hpp"

#include "thermo_parameters.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <optional>
#include <stdexcept>
#include <string>

namespace intarnanew {
namespace {

// ViennaRNA uses GASCONST = 1.98717 cal/(K mol).
constexpr Energy viennaGasConstantKcal = 0.00198717;

enum class PairType : std::size_t { cg, gc, gu, ug, au, ua, ambiguous };

[[nodiscard]] auto pairType(const char target, const char query) noexcept -> PairType {
    const auto left = static_cast<char>(std::toupper(static_cast<unsigned char>(target)));
    const auto right = static_cast<char>(std::toupper(static_cast<unsigned char>(query)));
    if (left == 'C' && right == 'G') return PairType::cg;
    if (left == 'G' && right == 'C') return PairType::gc;
    if (left == 'G' && right == 'U') return PairType::gu;
    if (left == 'U' && right == 'G') return PairType::ug;
    if (left == 'A' && right == 'U') return PairType::au;
    if (left == 'U' && right == 'A') return PairType::ua;
    return PairType::ambiguous;
}

[[nodiscard]] constexpr auto pairIndex(const PairType type) noexcept -> std::size_t {
    return static_cast<std::size_t>(type);
}

[[nodiscard]] constexpr auto reversePair(const PairType type) noexcept -> PairType {
    constexpr std::array reversed{
        PairType::gc, PairType::cg, PairType::ug, PairType::gu,
        PairType::ua, PairType::au, PairType::ambiguous,
    };
    return reversed[pairIndex(type)];
}

[[nodiscard]] auto baseIndex(const char nucleotide) noexcept -> std::size_t {
    switch (static_cast<char>(std::toupper(static_cast<unsigned char>(nucleotide)))) {
        case 'A': return 1U;
        case 'C': return 2U;
        case 'G': return 3U;
        case 'U': return 4U;
        default: return 0U;
    }
}

[[nodiscard]] auto canonicalBaseIndex(const char nucleotide) noexcept -> std::optional<std::size_t> {
    switch (static_cast<char>(std::toupper(static_cast<unsigned char>(nucleotide)))) {
        case 'A': return 0U;
        case 'C': return 1U;
        case 'G': return 2U;
        case 'U': return 3U;
        default: return std::nullopt;
    }
}

[[nodiscard]] auto stackValue(const detail::NearestNeighborParameters& parameters,
                              const PairType outer, const PairType inner) noexcept -> Energy {
    return parameters.stack[pairIndex(outer) * 7U + pairIndex(inner)];
}

[[nodiscard]] auto mismatchValue(const std::vector<Energy>& table, const PairType type,
                                 const char first, const char second) noexcept -> Energy {
    return table[(pairIndex(type) * 5U + baseIndex(first)) * 5U + baseIndex(second)];
}

[[nodiscard]] auto int11Value(const detail::NearestNeighborParameters& parameters,
                              const PairType outer, const PairType inner,
                              const char targetBase, const char queryBase) noexcept -> Energy {
    const auto index = (((pairIndex(outer) * 7U + pairIndex(inner)) * 5U +
                         baseIndex(targetBase)) * 5U + baseIndex(queryBase));
    return parameters.int11[index];
}

[[nodiscard]] auto int21Value(const detail::NearestNeighborParameters& parameters,
                              const PairType outer, const PairType inner,
                              const char first, const char second, const char third) noexcept -> Energy {
    const auto index = ((((pairIndex(outer) * 7U + pairIndex(inner)) * 5U +
                          baseIndex(first)) * 5U + baseIndex(second)) * 5U + baseIndex(third));
    return parameters.int21[index];
}

[[nodiscard]] auto int22Value(const detail::NearestNeighborParameters& parameters,
                              const PairType outer, const PairType inner,
                              const char first, const char second,
                              const char third, const char fourth) noexcept -> std::optional<Energy> {
    if (pairIndex(outer) >= 6U || pairIndex(inner) >= 6U) return std::nullopt;
    const auto firstIndex = canonicalBaseIndex(first);
    const auto secondIndex = canonicalBaseIndex(second);
    const auto thirdIndex = canonicalBaseIndex(third);
    const auto fourthIndex = canonicalBaseIndex(fourth);
    if (!firstIndex || !secondIndex || !thirdIndex || !fourthIndex) return std::nullopt;
    const auto index = (((((pairIndex(outer) * 6U + pairIndex(inner)) * 4U + *firstIndex) * 4U +
                           *secondIndex) * 4U + *thirdIndex) * 4U + *fourthIndex);
    return parameters.int22[index];
}

[[nodiscard]] auto dangleValue(const std::vector<Energy>& table, const PairType type,
                               const char nucleotide) noexcept -> Energy {
    const auto value = table[pairIndex(type) * 5U + baseIndex(nucleotide)];
    return std::isfinite(value) ? std::min(0.0, value) : 0.0;
}

[[nodiscard]] auto loopInitiation(const std::vector<Energy>& table, const Index size,
                                  const Energy logarithmicSlope) noexcept -> Energy {
    if (size < table.size()) return table[size];
    return table.back() + logarithmicSlope *
        std::log(static_cast<double>(size) / static_cast<double>(table.size() - 1U));
}

[[nodiscard]] auto validPath(const std::span<const BasePair> pairs) noexcept -> bool {
    if (pairs.empty()) return false;
    for (Index index = 1U; index < pairs.size(); ++index) {
        if (pairs[index - 1U].target >= pairs[index].target ||
            pairs[index - 1U].query <= pairs[index].query) return false;
    }
    return true;
}

[[nodiscard]] constexpr auto weakPair(const PairType type) noexcept -> bool {
    return type == PairType::gu || type == PairType::ug ||
           type == PairType::au || type == PairType::ua;
}

} // namespace

BasePairEnergyModel::BasePairEnergyModel(const double temperatureCelsius) {
    static_cast<void>(temperatureCelsius);
}

auto BasePairEnergyModel::evaluate(
    const Sequence&,
    const Sequence&,
    const std::span<const BasePair> pairs) const -> EnergyBreakdown {
    if (!validPath(pairs)) throw std::invalid_argument("interaction base pairs are not a monotone antiparallel path");
    EnergyBreakdown result;
    result.initiation = initiationEnergy();
    result.loops = -static_cast<Energy>(pairs.size() - 1U);
    return result;
}

auto BasePairEnergyModel::transitionEnergy(
    const Sequence&,
    const Sequence&,
    const BasePair,
    const BasePair) const -> Energy {
    return -1.0;
}

NearestNeighborEnergyModel::NearestNeighborEnergyModel(
    const double temperatureCelsius,
    const std::string_view parameterSet,
    const bool includeDangles)
    : rt_(viennaGasConstantKcal * (temperatureCelsius + 273.15)),
      includeDangles_(includeDangles),
      parameters_(detail::loadNearestNeighborParameters(temperatureCelsius, parameterSet)) {
    initiation_ = parameters_->duplexInit;
}

auto NearestNeighborEnergyModel::transitionEnergy(
    const Sequence& target,
    const Sequence& query,
    const BasePair left,
    const BasePair right) const -> Energy {
    if (left.target >= right.target || left.query <= right.query) {
        throw std::invalid_argument("invalid interaction transition orientation");
    }
    const auto targetGap = right.target - left.target - 1U;
    const auto queryGap = left.query - right.query - 1U;
    const auto outer = pairType(target[left.target], query[left.query]);
    const auto inner = reversePair(pairType(target[right.target], query[right.query]));
    if (targetGap == 0U && queryGap == 0U) {
        return stackValue(*parameters_, outer, inner);
    }

    const auto total = targetGap + queryGap;
    if (targetGap == 0U || queryGap == 0U) {
        auto energy = loopInitiation(parameters_->bulge, total, parameters_->logarithmicLoopSlope);
        if (total == 1U) {
            energy += stackValue(*parameters_, outer, inner);
        } else {
            if (weakPair(outer)) energy += parameters_->terminalAu;
            if (weakPair(inner)) energy += parameters_->terminalAu;
        }
        return energy;
    }

    const auto targetOuter = target[left.target + 1U];
    const auto targetInner = target[right.target - 1U];
    const auto queryOuter = query[left.query - 1U];
    const auto queryInner = query[right.query + 1U];
    if (targetGap == 1U && queryGap == 1U) {
        return int11Value(*parameters_, outer, inner, targetOuter, queryOuter);
    }
    if (targetGap == 1U && queryGap == 2U) {
        return int21Value(*parameters_, outer, inner, targetOuter, queryInner, queryOuter);
    }
    if (targetGap == 2U && queryGap == 1U) {
        return int21Value(*parameters_, inner, outer, queryInner, targetOuter, targetInner);
    }
    if (targetGap == 2U && queryGap == 2U) {
        if (const auto exact = int22Value(*parameters_, outer, inner,
                                         targetOuter, targetInner, queryInner, queryOuter)) {
            return *exact;
        }
    }

    auto energy = loopInitiation(parameters_->internal, total, parameters_->logarithmicLoopSlope);
    const auto asymmetry = targetGap > queryGap ? targetGap - queryGap : queryGap - targetGap;
    energy += std::min(parameters_->ninioMaximum,
                       parameters_->ninioSlope * static_cast<Energy>(asymmetry));
    const std::vector<Energy>* mismatch = &parameters_->mismatchInternal;
    if (targetGap == 1U || queryGap == 1U) {
        mismatch = &parameters_->mismatchInternal1n;
    } else if ((targetGap == 2U && queryGap == 3U) ||
               (targetGap == 3U && queryGap == 2U)) {
        mismatch = &parameters_->mismatchInternal23;
    }
    energy += mismatchValue(*mismatch, outer, targetOuter, queryOuter);
    energy += mismatchValue(*mismatch, inner, queryInner, targetInner);
    return energy;
}

auto NearestNeighborEnergyModel::terminalPenalty(const char target, const char query) const noexcept -> Energy {
    return weakPair(pairType(target, query)) ? parameters_->terminalAu : 0.0;
}

auto NearestNeighborEnergyModel::evaluate(
    const Sequence& target,
    const Sequence& query,
    const std::span<const BasePair> pairs) const -> EnergyBreakdown {
    if (!validPath(pairs)) throw std::invalid_argument("interaction base pairs are not a monotone antiparallel path");
    EnergyBreakdown result;
    result.initiation = initiation_;
    for (Index index = 1U; index < pairs.size(); ++index) {
        result.loops += transitionEnergy(target, query, pairs[index - 1U], pairs[index]);
    }
    const auto leftType = pairType(target[pairs.front().target], query[pairs.front().query]);
    const auto rightType = reversePair(pairType(target[pairs.back().target], query[pairs.back().query]));
    result.endLeft = terminalPenalty(target[pairs.front().target], query[pairs.front().query]);
    result.endRight = terminalPenalty(target[pairs.back().target], query[pairs.back().query]);

    if (includeDangles_) {
        const auto& left = pairs.front();
        const auto& right = pairs.back();
        const bool hasLeftTarget = left.target > 0U;
        const bool hasLeftQuery = left.query + 1U < query.size();
        const bool hasRightTarget = right.target + 1U < target.size();
        const bool hasRightQuery = right.query > 0U;
        if (hasLeftTarget && hasLeftQuery) {
            result.dangleLeft = std::min(0.0,
                mismatchValue(parameters_->mismatchExterior, leftType,
                    target[left.target - 1U], query[left.query + 1U]));
        } else if (hasLeftTarget) {
            result.dangleLeft = std::min(result.dangleLeft,
                dangleValue(parameters_->dangle5, leftType, target[left.target - 1U]));
        } else if (hasLeftQuery) {
            result.dangleLeft = std::min(result.dangleLeft,
                dangleValue(parameters_->dangle3, leftType, query[left.query + 1U]));
        }
        if (hasRightTarget && hasRightQuery) {
            result.dangleRight = std::min(0.0,
                mismatchValue(parameters_->mismatchExterior, rightType,
                    query[right.query - 1U], target[right.target + 1U]));
        } else if (hasRightTarget) {
            result.dangleRight = std::min(result.dangleRight,
                dangleValue(parameters_->dangle3, rightType, target[right.target + 1U]));
        } else if (hasRightQuery) {
            result.dangleRight = std::min(result.dangleRight,
                dangleValue(parameters_->dangle5, rightType, query[right.query - 1U]));
        }
    }
    return result;
}

auto makeEnergyModel(const Config& config) -> std::unique_ptr<HybridEnergyModel> {
    if (config.energy == EnergyKind::basePair) {
        return std::make_unique<BasePairEnergyModel>(config.temperatureCelsius);
    }
    return std::make_unique<NearestNeighborEnergyModel>(
        config.temperatureCelsius, config.energyParameters, !config.noDangles);
}

} // namespace intarnanew
