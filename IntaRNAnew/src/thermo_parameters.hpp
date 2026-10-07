#pragma once

// Internal, temperature-scaled representation of ViennaRNA v2 parameters.

#include "intarnanew/types.hpp"

#include <memory>
#include <string_view>
#include <vector>

namespace intarnanew::detail {

struct NearestNeighborParameters {
    std::vector<Energy> stack;
    std::vector<Energy> bulge;
    std::vector<Energy> internal;
    std::vector<Energy> mismatchInternal;
    std::vector<Energy> mismatchInternal1n;
    std::vector<Energy> mismatchInternal23;
    std::vector<Energy> mismatchExterior;
    std::vector<Energy> int11;
    std::vector<Energy> int21;
    std::vector<Energy> int22;
    std::vector<Energy> dangle5;
    std::vector<Energy> dangle3;
    Energy terminalAu{};
    Energy duplexInit{};
    Energy ninioSlope{};
    Energy ninioMaximum{};
    Energy logarithmicLoopSlope{};
};

[[nodiscard]] auto loadNearestNeighborParameters(
    double temperatureCelsius,
    std::string_view parameterSet) -> std::shared_ptr<const NearestNeighborParameters>;

} // namespace intarnanew::detail
