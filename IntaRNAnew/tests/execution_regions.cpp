#include "intarnanew/accessibility.hpp"
#include "intarnanew/cli.hpp"
#include "intarnanew/predictor.hpp"
#include "intarnanew/sequence.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <span>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

namespace {

int failureCount{};

void check(const bool condition, const std::string_view description) {
    if (!condition) {
        std::cerr << "FAILED: " << description << '\n';
        ++failureCount;
    }
}

[[nodiscard]] auto config(std::initializer_list<std::string_view> arguments)
    -> intarnanew::Config {
    const std::vector<std::string_view> values(arguments);
    auto parsed = intarnanew::Cli::parse(values);
    if (!parsed) {
        std::cerr << "configuration failed: " << parsed.error() << '\n';
        std::exit(2);
    }
    return std::move(*parsed);
}

[[nodiscard]] auto signatures(const intarnanew::PredictionResult& prediction)
    -> std::vector<std::tuple<intarnanew::Interval, intarnanew::Interval, double>> {
    std::vector<std::tuple<intarnanew::Interval, intarnanew::Interval, double>> result;
    result.reserve(prediction.interactions.size());
    for (const auto& interaction : prediction.interactions) {
        result.emplace_back(
            interaction.targetRange(), interaction.queryRange(), interaction.energy.total());
    }
    return result;
}

} // namespace

auto main() -> int {
    using namespace intarnanew;

    {
        const Sequence sequence("signed", "AAAAAAAA", -5);
        auto parsed = parseIntervals(sequence, "-5--3,-2-2");
        check(parsed && *parsed == std::vector<Interval>{{0U, 2U}, {3U, 6U}},
              "signed external regions map across the skipped zero coordinate");
        check(!parseIntervals(sequence, "-5--2,-3-1"),
              "overlapping manual regions are rejected independent of insertion order");
        check(!parseIntervals(sequence, "-5--3,"),
              "a trailing empty region is rejected");
        check(!parseIntervals(sequence, "0-1"),
              "external coordinate zero is rejected for a negative origin");
    }

    {
        const auto windows = decomposeWindows({0U, 13U}, 10U, 2U);
        check(windows == std::vector<Interval>{{0U, 9U}, {8U, 13U}},
              "window decomposition uses a fixed step and clipped final window");
        bool allCovered = true;
        for (Index begin = 0U; begin <= 13U; ++begin) {
            for (Index length = 1U; length <= 2U && begin + length <= 14U; ++length) {
                const Interval candidate{begin, begin + length - 1U};
                if (!std::ranges::any_of(windows, [&](const Interval window) {
                        return window.contains(candidate.begin) && window.contains(candidate.end);
                    })) {
                    allCovered = false;
                }
            }
        }
        check(allCovered, "every overlap-bounded interaction is contained in a window");
    }

    {
        const Sequence sequence("automatic", "CCCCCCCCCCCC");
        DisabledAccessibility accessibility(sequence);
        const std::vector<Interval> whole{{0U, sequence.size() - 1U}};
        const auto regions = decomposeAccessibleRegions(whole, 4U, 2U, accessibility);
        check(regions == std::vector<Interval>{{8U, 11U}},
              "automatic-region accessibility ties use the deterministic left-most cut");
    }

    {
        const std::vector<std::string_view> enabled{
            "--target=CCCCCCCCCCCCCC", "--query=GGGGGGGGGGGG", "--energy=B",
            "--acc=N", "--noSeed", "--intLenMax=2", "--windowWidth=10",
            "--windowOverlap=2", "--tRegionLenMax=4", "--outPerRegion",
        };
        auto parsed = Cli::parse(enabled);
        check(parsed && parsed->windowWidth == 10U && parsed->windowOverlap == 2U &&
                  parsed->target.regionLengthMax == 4U && parsed->output.perRegion,
              "region, window, and per-region CLI options are enabled");

        const std::vector<std::string_view> mixedRegions{
            "--target=CCCC", "--query=GG", "--energy=B", "--acc=N", "--noSeed",
            "--tRegion=1-4", "--tRegionLenMax=3",
        };
        check(!Cli::parse(mixedRegions), "manual and automatic regions are mutually exclusive");

        const std::vector<std::string_view> unboundedWindow{
            "--target=CCCCCCCCCCCC", "--query=GGGGGGGGGG", "--energy=B", "--acc=N",
            "--noSeed", "--windowWidth=10", "--windowOverlap=2",
        };
        check(!Cli::parse(unboundedWindow),
              "windowing rejects an unbounded interaction length");

        const std::vector<std::string_view> ensembleWindow{
            "--target=CCCCCCCCCCCC", "--query=GGGGGGGGGG", "--energy=B", "--acc=N",
            "--noSeed", "--model=P", "--intLenMax=2", "--windowWidth=10",
            "--windowOverlap=2",
        };
        check(!Cli::parse(ensembleWindow),
              "partition-site mode is rejected across overlapping windows");
    }

    {
        auto run = config({
            "--target=CCCCCCCCCCCCCC", "--query=GGGGGGGGGGGG", "--energy=B",
            "--acc=N", "--noSeed", "--mode=M", "--model=S", "--intLenMax=2",
            "--outNumber=1000", "--outOverlap=B", "--outDeltaE=100",
        });
        const Sequence target("target", "CCCCCCCCCCCCCC");
        const Sequence query("query", "GGGGGGGGGGGG");
        DisabledAccessibility targetAccessibility(target);
        DisabledAccessibility queryAccessibility(query);
        const Predictor predictor(run);
        const auto full = predictor.predict(
            target, query, targetAccessibility, queryAccessibility);

        std::vector<PredictionResult> domains;
        const auto targetWindows = decomposeWindows({0U, target.size() - 1U}, 10U, 2U);
        const auto queryWindows = decomposeWindows({0U, query.size() - 1U}, 10U, 2U);
        for (const auto targetWindow : targetWindows) {
            for (const auto queryWindow : queryWindows) {
                domains.push_back(predictor.predict(
                    target, query, targetAccessibility, queryAccessibility,
                    targetWindow, queryWindow));
            }
        }
        const auto reduced = reducePredictions(domains, run);
        check(signatures(full) == signatures(reduced),
              "windowed real-interaction sites equal the non-windowed prediction");
        check(std::abs(full.logPartition - reduced.logPartition) < 1e-12,
              "overlapping windows are deduplicated before partition reduction");
        check(std::ranges::any_of(reduced.interactions, [](const Interaction& interaction) {
            return interaction.targetRange() == Interval{9U, 10U};
        }), "an interaction crossing the first window end is retained");

        std::ranges::reverse(domains);
        const auto reverseReduced = reducePredictions(domains, run);
        check(signatures(reduced) == signatures(reverseReduced) &&
                  std::abs(reduced.logPartition - reverseReduced.logPartition) < 1e-12,
              "domain reduction is byte-order independent for parallel completion");
    }

    {
        auto run = config({
            "--target=CCCCCCCC", "--query=GGGG", "--energy=B", "--acc=N",
            "--noSeed", "--mode=M", "--model=S", "--outNumber=1",
            "--outOverlap=B", "--outDeltaE=100", "--outPerRegion",
        });
        const Sequence target("target", "CCCCCCCC");
        const Sequence query("query", "GGGG");
        DisabledAccessibility targetAccessibility(target);
        DisabledAccessibility queryAccessibility(query);
        const Predictor predictor(run);
        const std::vector<Interval> targetRegions{{0U, 1U}, {6U, 7U}};
        const std::vector<Interval> queryRegions{{0U, 1U}, {2U, 3U}};
        std::vector<PredictionResult> regionResults;
        for (const auto targetRegion : targetRegions) {
            for (const auto queryRegion : queryRegions) {
                auto prediction = predictor.predict(target, query,
                    targetAccessibility, queryAccessibility, targetRegion, queryRegion);
                check(prediction.interactions.size() == 1U &&
                          std::abs(prediction.interactions.front().energy.total() + 2.0) < 1e-12,
                      "each real region combination has an independent optimum");
                regionResults.push_back(std::move(prediction));
            }
        }
        const auto globallyReduced = reducePredictions(regionResults, run);
        check(globallyReduced.interactions.size() == 1U,
              "without per-region grouping the global output limit is applied once");
    }

    {
        auto run = config({
            "--target=CCCCCCCCCCCC", "--query=GG", "--energy=B", "--acc=N",
            "--noSeed", "--seedBP=2", "--mode=M", "--model=S",
            "--tRegionLenMax=4", "--qRegion=1-2", "--outNumber=1",
            "--outOverlap=B", "--outDeltaE=100",
        });
        const Sequence target("target", "CCCCCCCCCCCC");
        const Sequence query("query", "GG");
        DisabledAccessibility targetAccessibility(target);
        DisabledAccessibility queryAccessibility(query);
        const auto prediction = Predictor(run).predict(
            target, query, targetAccessibility, queryAccessibility);
        check(prediction.interactions.size() == 1U &&
                  prediction.interactions.front().targetRange() == Interval{8U, 9U},
              "automatic-region real interaction matches the legacy tie fixture");
    }

    {
        auto run = config({
            "--target=CCAACCCACC", "--query=GGGG", "--energy=B", "--acc=N",
            "--mode=H", "--model=X", "--seedBP=2", "--intLenMax=2",
            "--outNumber=10", "--outOverlap=Q", "--outDeltaE=0",
        });
        const Sequence target("target", "CCAACCCACC");
        const Sequence query("query", "GGGG");
        DisabledAccessibility targetAccessibility(target);
        DisabledAccessibility queryAccessibility(query);
        const auto prediction = Predictor(run).predict(
            target, query, targetAccessibility, queryAccessibility);
        std::vector<std::pair<Interval, Interval>> sites;
        for (const auto& interaction : prediction.interactions) {
            sites.emplace_back(interaction.targetRange(), interaction.queryRange());
        }
        check(sites == std::vector<std::pair<Interval, Interval>>{
                  {{0U, 1U}, {0U, 1U}},
                  {{4U, 5U}, {2U, 3U}},
                  {{8U, 9U}, {2U, 3U}},
              },
              "target-nonoverlap tie traversal matches the legacy RNA interaction corpus");
    }

    {
        auto run = config({
            "--target=GGGAAACCC", "--query=GGGAAACCC", "--energy=B", "--acc=N",
            "--noSeed", "--mode=M", "--model=P", "--outNumber=1",
        });
        const Sequence target("target", "GGGAAACCC");
        const Sequence query("query", "GGGAAACCC");
        DisabledAccessibility targetAccessibility(target);
        DisabledAccessibility queryAccessibility(query);
        const auto prediction = Predictor(run).predict(
            target, query, targetAccessibility, queryAccessibility);
        check(prediction.interactions.size() == 1U, "model-P fixture has a best site");
        if (!prediction.interactions.empty()) {
            const auto& interaction = prediction.interactions.front();
            check(interaction.pairs == std::vector<BasePair>{{0U, 8U}, {8U, 0U}},
                  "model P reports only its site boundary pairs");
            check(std::abs(interaction.energy.total() -
                           std::trunc(interaction.ensembleFreeEnergy * 100.0) / 100.0) < 1e-12,
                  "model-P reported energy is the centikcal-truncated site free energy");
            check(std::abs(interaction.energy.total() + 7.91) < 0.02,
                  "model-P base-pair fixture matches the legacy site energy");
        }
    }

    if (failureCount == 0) {
        std::cout << "All IntaRNAnew region/window execution tests passed.\n";
    }
    return failureCount == 0 ? 0 : 1;
}
