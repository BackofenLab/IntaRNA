#include "intarnanew/tools/pvalue.hpp"

#include "intarnanew/tools/mutations.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <expected>
#include <thread>

namespace intarnanew::tools {
namespace {

[[nodiscard]] auto mixedSeed(const std::uint64_t seed, const std::uint64_t stream) noexcept
    -> std::uint64_t {
    auto value = seed + 0x9E3779B97F4A7C15ULL * (stream + 1U);
    value = (value ^ (value >> 30U)) * 0xBF58476D1CE4E5B9ULL;
    value = (value ^ (value >> 27U)) * 0x94D049BB133111EBULL;
    return value ^ (value >> 31U);
}

[[nodiscard]] auto shuffled(
    const std::string_view sequence,
    const std::uint64_t seed,
    const ShufflePreservation preservation) -> std::expected<std::string, std::string> {
    if (preservation == ShufflePreservation::dinucleotide) {
        return shuffleDinucleotides(sequence, seed);
    }
    return shuffleMononucleotides(sequence, seed);
}

} // namespace

auto sampleRandomInteractionScores(
    const std::string_view query,
    const std::string_view target,
    const RandomScoreOptions& options,
    const InteractionScoreEvaluator& evaluator) -> std::expected<std::vector<double>, std::string> {
    if (query.empty() || target.empty()) {
        return std::unexpected("query and target must not be empty");
    }
    if (options.cardinality == 0U) {
        return std::unexpected("random-score cardinality must be positive");
    }
    if (!evaluator) {
        return std::unexpected("an interaction-score evaluator is required");
    }

    std::vector<std::string> queries(options.cardinality, std::string{query});
    std::vector<std::string> targets(options.cardinality, std::string{target});
    for (std::size_t index = 0U; index < options.cardinality; ++index) {
        if (options.mode == ShuffleMode::query || options.mode == ShuffleMode::both) {
            auto generated = shuffled(
                query, mixedSeed(options.randomSeed, 2U * index), options.preservation);
            if (!generated) {
                return std::unexpected("query shuffle " + std::to_string(index) + ": " + generated.error());
            }
            queries[index] = std::move(*generated);
        }
        if (options.mode == ShuffleMode::target || options.mode == ShuffleMode::both) {
            auto generated = shuffled(
                target, mixedSeed(options.randomSeed, 2U * index + 1U), options.preservation);
            if (!generated) {
                return std::unexpected("target shuffle " + std::to_string(index) + ": " + generated.error());
            }
            targets[index] = std::move(*generated);
        }
    }

    std::vector<std::expected<double, std::string>> results(options.cardinality);
    std::atomic<std::size_t> next{};
    const auto requested = options.threads == 0U
        ? std::max(1U, std::thread::hardware_concurrency())
        : options.threads;
    const auto threadCount = std::min<std::size_t>(requested, options.cardinality);
    {
        std::vector<std::jthread> workers;
        workers.reserve(threadCount);
        for (std::size_t worker = 0U; worker < threadCount; ++worker) {
            workers.emplace_back([&] {
                for (;;) {
                    const auto index = next.fetch_add(1U, std::memory_order_relaxed);
                    if (index >= options.cardinality) {
                        break;
                    }
                    try {
                        results[index] = evaluator(queries[index], targets[index]);
                    } catch (const std::exception& error) {
                        results[index] = std::unexpected(
                            "score evaluator threw an exception: " + std::string{error.what()});
                    } catch (...) {
                        results[index] = std::unexpected("score evaluator threw an unknown exception");
                    }
                }
            });
        }
    }

    std::vector<double> scores;
    scores.reserve(options.cardinality);
    for (std::size_t index = 0U; index < results.size(); ++index) {
        if (!results[index]) {
            return std::unexpected(
                "score evaluation " + std::to_string(index) + ": " + results[index].error());
        }
        if (!std::isfinite(*results[index])) {
            return std::unexpected("score evaluation returned a non-finite value");
        }
        scores.push_back(*results[index]);
    }
    return scores;
}

} // namespace intarnanew::tools
