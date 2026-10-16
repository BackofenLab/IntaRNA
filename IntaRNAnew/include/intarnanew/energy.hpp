#pragma once

#include "intarnanew/config.hpp"
#include "intarnanew/sequence.hpp"
#include "intarnanew/types.hpp"

#include <memory>
#include <span>
#include <string_view>

namespace intarnanew {

class HybridEnergyModel {
public:
    virtual ~HybridEnergyModel() = default;

    [[nodiscard]] virtual auto evaluate(
        const Sequence& target,
        const Sequence& query,
        std::span<const BasePair> pairs) const -> EnergyBreakdown = 0;

    [[nodiscard]] virtual auto transitionEnergy(
        const Sequence& target,
        const Sequence& query,
        BasePair left,
        BasePair right) const -> Energy = 0;

    [[nodiscard]] virtual auto initiationEnergy() const noexcept -> Energy = 0;
    [[nodiscard]] virtual auto rt() const noexcept -> double = 0;
};

class BasePairEnergyModel final : public HybridEnergyModel {
public:
    explicit BasePairEnergyModel(double temperatureCelsius);

    [[nodiscard]] auto evaluate(
        const Sequence& target,
        const Sequence& query,
        std::span<const BasePair> pairs) const -> EnergyBreakdown override;
    [[nodiscard]] auto transitionEnergy(
        const Sequence& target,
        const Sequence& query,
        BasePair left,
        BasePair right) const -> Energy override;
    [[nodiscard]] auto initiationEnergy() const noexcept -> Energy override { return -1.0; }
    [[nodiscard]] auto rt() const noexcept -> double override { return 1.0; }
};

namespace detail { struct NearestNeighborParameters; }

class NearestNeighborEnergyModel final : public HybridEnergyModel {
public:
    NearestNeighborEnergyModel(
        double temperatureCelsius,
        std::string_view parameterSet,
        bool includeDangles);

    [[nodiscard]] auto evaluate(
        const Sequence& target,
        const Sequence& query,
        std::span<const BasePair> pairs) const -> EnergyBreakdown override;
    [[nodiscard]] auto transitionEnergy(
        const Sequence& target,
        const Sequence& query,
        BasePair left,
        BasePair right) const -> Energy override;
    [[nodiscard]] auto initiationEnergy() const noexcept -> Energy override { return initiation_; }
    [[nodiscard]] auto rt() const noexcept -> double override { return rt_; }

private:
    [[nodiscard]] auto terminalPenalty(char target, char query) const noexcept -> Energy;

    double rt_{};
    Energy initiation_{};
    bool includeDangles_{true};
    std::shared_ptr<const detail::NearestNeighborParameters> parameters_;
};

[[nodiscard]] auto makeEnergyModel(const Config& config) -> std::unique_ptr<HybridEnergyModel>;

} // namespace intarnanew
