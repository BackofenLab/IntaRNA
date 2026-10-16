#pragma once

#include "intarnanew/accessibility.hpp"
#include "intarnanew/config.hpp"
#include "intarnanew/predictor.hpp"
#include "intarnanew/sequence.hpp"

#include <expected>
#include <string>
#include <string_view>

namespace intarnanew {

class OutputFormatter {
public:
    [[nodiscard]] static auto primary(
        const Config& config,
        const Sequence& target,
        const Sequence& query,
        const PredictionResult& result,
        bool includeCsvHeader = true) -> std::expected<std::string, std::string>;

    [[nodiscard]] static auto auxiliary(
        std::string_view descriptor,
        const Config& config,
        const Sequence& target,
        const Sequence& query,
        const AccessibilityProvider& targetAccessibility,
        const AccessibilityProvider& queryAccessibility,
        const PredictionResult& result) -> std::expected<std::string, std::string>;
};

[[nodiscard]] auto writeOutput(std::string_view destination, std::string_view content)
    -> std::expected<void, std::string>;

} // namespace intarnanew
