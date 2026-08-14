#pragma once

#include "intarnanew/types.hpp"

#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace intarnanew::folding_detail {

struct FoldingParameters {
    std::vector<Energy> stack;
    std::vector<Energy> hairpin;
    std::vector<Energy> bulge;
    std::vector<Energy> internal;
    std::vector<Energy> mismatchHairpin;
    std::vector<Energy> mismatchInternal;
    std::vector<Energy> mismatchInternal1n;
    std::vector<Energy> mismatchInternal23;
    std::vector<Energy> mismatchMulti;
    std::vector<Energy> mismatchExterior;
    std::vector<Energy> dangle5;
    std::vector<Energy> dangle3;
    std::vector<Energy> int11;
    std::vector<Energy> int21;
    std::vector<Energy> int22;
    std::unordered_map<std::string, Energy> specialHairpins;
    Energy multiloopUnpaired{};
    Energy multiloopClosing{};
    Energy multiloopStem{};
    Energy terminalAu{};
    Energy ninioSlope{};
    Energy ninioMaximum{};
    Energy logarithmicLoopSlope{};
};

[[nodiscard]] auto loadFoldingParameters(
    double temperatureCelsius,
    std::string_view parameterSet) -> std::shared_ptr<const FoldingParameters>;

} // namespace intarnanew::folding_detail
