#pragma once

#include "intarnanew/accessibility.hpp"
#include "intarnanew/config.hpp"
#include "intarnanew/predictor.hpp"
#include "intarnanew/sequence.hpp"

#include <expected>
#include <memory>
#include <string>

namespace intarnanew {

// Complete in-process evaluation of one target/query pair. The providers are
// retained so native ED/Pu profiles can be consumed by companion tools without
// rebuilding either folding ensemble. The referenced Sequence objects must
// outlive the returned providers.
struct PairPrediction {
    PredictionResult prediction;
    std::unique_ptr<AccessibilityProvider> targetAccessibility;
    std::unique_ptr<AccessibilityProvider> queryAccessibility;
};

// Applies the same explicit/automatic region planning, bounded window
// decomposition, site deduplication, global output constraints, and monomer
// summaries as a single pair in the command-line application. Domain work is
// deliberately sequential: callers such as randomization and mutation tools
// parallelize independent pairs without creating nested worker pools.
[[nodiscard]] auto predictPair(
    const Config& config,
    const Sequence& target,
    const Sequence& query) -> std::expected<PairPrediction, std::string>;

} // namespace intarnanew
