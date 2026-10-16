#pragma once

#include <cstddef>
#include <cstdint>
#include <expected>
#include <functional>
#include <string>
#include <string_view>
#include <vector>

namespace intarnanew::tools {

enum class ShuffleMode { query, target, both };
enum class ShufflePreservation { mononucleotide, dinucleotide };

struct RandomScoreOptions {
    std::size_t cardinality{};
    ShuffleMode mode{ShuffleMode::both};
    ShufflePreservation preservation{ShufflePreservation::dinucleotide};
    std::uint64_t randomSeed{};
    // Zero selects hardware_concurrency, bounded to cardinality.
    std::size_t threads{};
};

using InteractionScoreEvaluator = std::function<
    std::expected<double, std::string>(std::string_view query, std::string_view target)>;

// Clean orchestration boundary for IntaRNApvalue-style sampling. It never
// starts an external process: callers inject the C++ prediction-library
// evaluator. Random inputs and result ordering are invariant to thread count.
// The evaluator must be safe for simultaneous calls when threads > 1.
[[nodiscard]] auto sampleRandomInteractionScores(
    std::string_view query,
    std::string_view target,
    const RandomScoreOptions& options,
    const InteractionScoreEvaluator& evaluator) -> std::expected<std::vector<double>, std::string>;

} // namespace intarnanew::tools
