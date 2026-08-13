#include "intarnanew/cli.hpp"
#include "intarnanew/runner.hpp"
#include "intarnanew/sequence.hpp"
#include "intarnanew/tools/pvalue.hpp"
#include "intarnanew/tools/statistics.hpp"

#include <charconv>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <exception>
#include <iomanip>
#include <iostream>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using intarnanew::Config;
using intarnanew::Sequence;
using intarnanew::tools::DistributionKind;
using intarnanew::tools::RandomScoreOptions;
using intarnanew::tools::ShuffleMode;

struct Arguments {
    std::string query;
    std::string target;
    std::string parameterFile;
    std::size_t samples{};
    ShuffleMode shuffle{ShuffleMode::both};
    DistributionKind distribution{DistributionKind::gev};
    bool scores{};
    std::size_t threads{};
    std::uint64_t seed{};
};

[[nodiscard]] auto help() -> std::string_view {
    return R"(intarnanew-pvalue - native randomization-based interaction p-values

Usage:
  intarnanew-pvalue -q RNA_OR_FASTA -t RNA_OR_FASTA -s N -m q|t|b [options]

Options:
  -q, --query INPUT          query RNA or FASTA/gzip file (first record)
  -t, --target INPUT         target RNA or FASTA/gzip file (first record)
  -s, -c, --samples N        positive random-sample count
      --cardinality, --scores N
                              aliases for the random-sample count
  -m, --shuffle-mode q|t|b   shuffle query, target, or both
  -p, --parameterFile FILE   IntaRNAnew parameter file for scoring
  -d, --distribution KIND    gev, gumbel, or gauss (default: gev)
  -o, --output KIND          pvalue or scores (default: pvalue)
      --threads N            worker count; zero uses available hardware
      --randSeed N           deterministic unsigned random seed (default: 0)
  -h, --help                 show this help

All predictions run in-process through the C++ library; no executable or
script is launched. Randomized sequences and result order are thread invariant.
)";
}

template <class Integer>
[[nodiscard]] auto integer(
    const std::string_view text,
    const std::string_view option) -> std::expected<Integer, std::string> {
    Integer value{};
    const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
    if (text.empty() || error != std::errc{} || end != text.data() + text.size()) {
        return std::unexpected(std::string(option) + " requires an unsigned integer");
    }
    return value;
}

[[nodiscard]] auto parse(std::span<const std::string_view> values)
    -> std::expected<std::optional<Arguments>, std::string> {
    if (values.empty()) return std::unexpected("missing arguments; use --help");
    Arguments result;
    for (std::size_t index = 0U; index < values.size(); ++index) {
        auto option = values[index];
        std::optional<std::string_view> inlineValue;
        if (option.starts_with("--")) {
            if (const auto separator = option.find('='); separator != std::string_view::npos) {
                inlineValue = option.substr(separator + 1U);
                option = option.substr(0U, separator);
            }
        }
        if (option == "-h" || option == "--help") return std::optional<Arguments>{};
        const auto next = [&]() -> std::expected<std::string_view, std::string> {
            if (inlineValue) {
                const auto value = *inlineValue;
                inlineValue.reset();
                return value;
            }
            if (++index >= values.size()) {
                return std::unexpected("missing value for " + std::string(option));
            }
            return values[index];
        };
        if (option == "-q" || option == "--query") {
            auto value = next(); if (!value) return std::unexpected(value.error());
            result.query = *value;
        } else if (option == "-t" || option == "--target") {
            auto value = next(); if (!value) return std::unexpected(value.error());
            result.target = *value;
        } else if (option == "-s" || option == "-c" || option == "--samples" ||
                   option == "--cardinality" || option == "--scores") {
            auto value = next(); if (!value) return std::unexpected(value.error());
            auto parsed = integer<std::size_t>(*value, option);
            if (!parsed || *parsed == 0U) {
                return std::unexpected(parsed ? "--samples must be positive" : parsed.error());
            }
            result.samples = *parsed;
        } else if (option == "-m" || option == "--shuffle-mode") {
            auto value = next(); if (!value) return std::unexpected(value.error());
            if (*value == "q") result.shuffle = ShuffleMode::query;
            else if (*value == "t") result.shuffle = ShuffleMode::target;
            else if (*value == "b") result.shuffle = ShuffleMode::both;
            else return std::unexpected("shuffle mode must be q, t, or b");
        } else if (option == "-p" || option == "--parameterFile") {
            auto value = next(); if (!value) return std::unexpected(value.error());
            result.parameterFile = *value;
        } else if (option == "-d" || option == "--distribution") {
            auto value = next(); if (!value) return std::unexpected(value.error());
            auto parsed = intarnanew::tools::parseDistribution(*value);
            if (!parsed) return std::unexpected(parsed.error());
            result.distribution = *parsed;
        } else if (option == "-o" || option == "--output") {
            auto value = next(); if (!value) return std::unexpected(value.error());
            if (*value == "scores") result.scores = true;
            else if (*value == "pvalue") result.scores = false;
            else return std::unexpected("output must be pvalue or scores");
        } else if (option == "--threads") {
            auto value = next(); if (!value) return std::unexpected(value.error());
            auto parsed = integer<std::size_t>(*value, option);
            if (!parsed) return std::unexpected(parsed.error());
            result.threads = *parsed;
        } else if (option == "--randSeed") {
            auto value = next(); if (!value) return std::unexpected(value.error());
            auto parsed = integer<std::uint64_t>(*value, option);
            if (!parsed) return std::unexpected(parsed.error());
            result.seed = *parsed;
        } else {
            return std::unexpected("unknown option '" + std::string(option) + "'");
        }
    }
    if (result.query.empty()) return std::unexpected("--query is required");
    if (result.target.empty()) return std::unexpected("--target is required");
    if (result.samples == 0U) return std::unexpected("--samples is required");
    return std::optional<Arguments>{std::move(result)};
}

[[nodiscard]] auto scoringConfig(const Arguments& arguments)
    -> std::expected<Config, std::string> {
    std::vector<std::string> storage;
    storage.reserve(5U);
    if (!arguments.parameterFile.empty()) {
        storage.push_back("--parameterFile=" + arguments.parameterFile);
    }
    storage.emplace_back("--target=" + arguments.target);
    storage.emplace_back("--query=" + arguments.query);
    storage.emplace_back("--outNumber=1");
    storage.emplace_back("--threads=1");
    std::vector<std::string_view> views;
    views.reserve(storage.size());
    for (const auto& value : storage) views.push_back(value);
    return intarnanew::Cli::parse(views, "IntaRNAnew");
}

[[nodiscard]] auto score(
    const Config& base,
    const std::string_view queryText,
    const std::string_view targetText) -> std::expected<double, std::string> {
    try {
        const Sequence target(base.target.id, std::string(targetText), base.target.firstPosition);
        const Sequence query(base.query.id, std::string(queryText), base.query.firstPosition);
        auto evaluation = intarnanew::predictPair(base, target, query);
        if (!evaluation) return std::unexpected(evaluation.error());
        return evaluation->prediction.interactions.empty()
            ? base.output.maxEnergy
            : evaluation->prediction.interactions.front().energy.total();
    } catch (const std::exception& error) {
        return std::unexpected(error.what());
    }
}

} // namespace

auto main(const int argc, const char* const* argv) -> int {
    std::vector<std::string_view> values;
    values.reserve(static_cast<std::size_t>(argc > 0 ? argc - 1 : 0));
    for (int index = 1; index < argc; ++index) values.emplace_back(argv[index]);
    auto arguments = parse(values);
    if (!arguments) {
        std::cerr << "intarnanew-pvalue: " << arguments.error() << '\n';
        return 2;
    }
    if (!*arguments) {
        std::cout << help();
        return 0;
    }
    auto config = scoringConfig(**arguments);
    if (!config) {
        std::cerr << "intarnanew-pvalue: " << config.error() << '\n';
        return 2;
    }
    auto targetInput = intarnanew::SequenceReader::read(
        (*arguments)->target, config->target.id, config->target.firstPosition, std::cin);
    if (!targetInput) {
        std::cerr << "intarnanew-pvalue: target input: " << targetInput.error() << '\n';
        return 2;
    }
    auto queryInput = intarnanew::SequenceReader::read(
        (*arguments)->query, config->query.id, config->query.firstPosition, std::cin);
    if (!queryInput) {
        std::cerr << "intarnanew-pvalue: query input: " << queryInput.error() << '\n';
        return 2;
    }
    if (targetInput->empty() || queryInput->empty()) {
        std::cerr << "intarnanew-pvalue: query and target inputs must contain a sequence\n";
        return 2;
    }
    const auto& target = targetInput->front();
    const auto& query = queryInput->front();
    config->target.id = target.id();
    config->query.id = query.id();

    auto observed = score(*config, query.str(), target.str());
    if (!observed) {
        std::cerr << "intarnanew-pvalue: observed score: " << observed.error() << '\n';
        return 1;
    }
    RandomScoreOptions options;
    options.cardinality = (*arguments)->samples;
    options.mode = (*arguments)->shuffle;
    options.randomSeed = (*arguments)->seed;
    options.threads = (*arguments)->threads;
    auto sampled = intarnanew::tools::sampleRandomInteractionScores(
        query.str(), target.str(), options,
        [&](const std::string_view shuffledQuery, const std::string_view shuffledTarget) {
            return score(*config, shuffledQuery, shuffledTarget);
        });
    if (!sampled) {
        std::cerr << "intarnanew-pvalue: " << sampled.error() << '\n';
        return 1;
    }

    std::cout << std::setprecision(17);
    if ((*arguments)->scores) {
        for (const auto value : *sampled) std::cout << value << '\n';
        return std::cout ? 0 : 1;
    }
    auto fitted = intarnanew::tools::fitDistribution(
        *sampled, (*arguments)->distribution);
    if (!fitted) {
        std::cerr << "intarnanew-pvalue: " << fitted.error() << '\n';
        return 1;
    }
    auto probability = intarnanew::tools::interactionEnergyPValue(*observed, *fitted);
    if (!probability) {
        std::cerr << "intarnanew-pvalue: " << probability.error() << '\n';
        return 1;
    }
    std::cout << *probability << '\n';
    return std::cout ? 0 : 1;
}
