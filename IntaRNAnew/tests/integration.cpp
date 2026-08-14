#include "intarnanew/accessibility.hpp"
#include "intarnanew/cli.hpp"
#include "intarnanew/energy.hpp"
#include "intarnanew/predictor.hpp"
#include "intarnanew/sequence.hpp"

#include <cmath>
#include <cstdlib>
#include <initializer_list>
#include <iostream>
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace {

auto failureCount = 0;

void check(const bool condition, const std::string_view description) {
    if (!condition) {
        std::cerr << "FAILED: " << description << '\n';
        ++failureCount;
    }
}

[[nodiscard]] auto isInterval(
    const intarnanew::Interval interval,
    const intarnanew::Index begin,
    const intarnanew::Index end) -> bool {
    return interval.begin == begin && interval.end == end;
}

[[nodiscard]] auto baseConfig(const std::vector<std::string_view>& arguments) -> intarnanew::Config {
    auto parsed = intarnanew::Cli::parse(arguments);
    if (!parsed) {
        std::cerr << "parse failed: " << parsed.error() << '\n';
        std::exit(2);
    }
    return std::move(*parsed);
}

[[nodiscard]] auto predict(
    const intarnanew::Config& config,
    const intarnanew::Sequence& target,
    const intarnanew::Sequence& query) -> intarnanew::PredictionResult {
    auto targetAccessibility = intarnanew::makeAccessibility(target, config.target, config.temperatureCelsius);
    auto queryAccessibility = intarnanew::makeAccessibility(query, config.query, config.temperatureCelsius);
    if (!targetAccessibility || !queryAccessibility) {
        std::cerr << "accessibility setup failed\n";
        std::exit(2);
    }
    return intarnanew::Predictor(config).predict(target, query, **targetAccessibility, **queryAccessibility);
}

} // namespace

auto main() -> int {
    using namespace intarnanew;

    {
        const auto config = baseConfig({"--target=CCAACACC", "--query=GG", "--energy=B", "--acc=N",
            "--noSeed", "--mode=M", "--model=S", "--outMode=C", "--outNumber=20",
            "--outOverlap=B", "--outDeltaE=0"});
        const Sequence target("target", "CCAACACC");
        const Sequence query("query", "GG");
        const auto result = predict(config, target, query);
        check(result.interactions.size() == 10U, "legacy real-interaction base-pair fixture has ten optima");
        check(std::ranges::all_of(result.interactions, [](const auto& interaction) {
            return std::abs(interaction.energy.total() + 2.0) < 1e-9;
        }), "all base-pair fixture optima have -2 kcal/mol");
        check(isInterval(result.interactions.front().targetRange(), 0U, 1U),
              "first base-pair fixture interaction uses target positions 1..2");
    }

    {
        const auto config = baseConfig({"--target=CCAACCCACC", "--query=GGGG", "--energy=B", "--acc=N",
            "--seedBP=3", "--mode=M", "--model=S", "--outMode=C", "--outNumber=10",
            "--outOverlap=B", "--outDeltaE=0", "--outNoLP=true"});
        const Sequence target("target", "CCAACCCACC");
        const Sequence query("query", "GGGG");
        const auto result = predict(config, target, query);
        check(result.interactions.size() == 2U, "legacy no-lonely-pair seed fixture has two optima");
        check(result.interactions.front().pairs.size() == 3U, "no-lonely-pair seed uses three pairs");
        check(isInterval(result.interactions.front().targetRange(), 4U, 6U),
              "no-lonely-pair seed targets C stretch");
    }

    {
        const auto config = baseConfig({"--target=CCAACCCACC", "--query=GGGG", "--energy=B", "--acc=N",
            "--seedBP=2", "--mode=H", "--model=X", "--tIntLenMax=2", "--qIntLenMax=2",
            "--outMode=C", "--outNumber=10", "--outOverlap=N", "--outDeltaE=0"});
        const Sequence target("target", "CCAACCCACC");
        const Sequence query("query", "GGGG");
        const auto result = predict(config, target, query);
        check(result.interactions.size() == 2U, "non-overlap policy reproduces two disjoint fixtures");
        check(!result.interactions[0].targetRange().overlaps(result.interactions[1].targetRange()),
              "selected target intervals do not overlap");
        check(!result.interactions[0].queryRange().overlaps(result.interactions[1].queryRange()),
              "selected query intervals do not overlap");
    }

    {
        const auto config = baseConfig({"--target=UUUUUUCCCUAAAUCCUUCUUU",
            "--query=GGGGGGGGGGGGGGGGGGGGG", "--energy=B", "--acc=N", "--noSeed",
            "--tRegion=8-17", "--qRegion=4-14", "--outNoGUend=true", "--outNumber=1",
            "--outDeltaE=0"});
        const Sequence target("target", "UUUUUUCCCUAAAUCCUUCUUU");
        const Sequence query("query", "GGGGGGGGGGGGGGGGGGGGG");
        const auto result = predict(config, target, query);
        check(result.interactions.size() == 1U, "no-GU-end fixture has one best interaction");
        if (!result.interactions.empty()) {
            check(result.interactions.front().pairs.size() == 4U,
                  "no-GU-end fixture has four base pairs");
            check(result.interactions.front().targetRange() == Interval{7U, 15U},
                  "no-GU-end fixture spans external target positions 8..16");
            check(std::abs(result.interactions.front().energy.total() + 4.0) < 1e-9,
                  "no-GU-end fixture has -4 kcal/mol");
        }
    }

    {
        const auto config = baseConfig({"--target=CGC", "--query=GUG", "--energy=B",
            "--acc=N", "--noSeed", "--mode=M", "--model=S", "--outNoGUend=true",
            "--outNumber=1", "--outDeltaE=0"});
        const Sequence target("target", "CGC");
        const Sequence query("query", "GUG");
        const auto result = predict(config, target, query);
        check(result.interactions.size() == 1U,
              "no-GU-end permits a GU pair stacked inside a helix");
        if (!result.interactions.empty()) {
            check(result.interactions.front().pairs ==
                  std::vector<BasePair>{{0U, 2U}, {1U, 1U}, {2U, 0U}},
                  "internal-GU fixture retains all three stacked pairs");
            check(std::abs(result.interactions.front().energy.total() + 3.0) < 1e-9,
                  "internal-GU base-pair fixture has -3 kcal/mol");
        }
    }

    {
        const NearestNeighborEnergyModel model(37.0, "Turner04", false);
        const auto evaluate = [&](const std::string_view targetText,
                                  const std::string_view queryText,
                                  const std::initializer_list<BasePair> pairs) {
            const Sequence target("target", std::string(targetText));
            const Sequence query("query", std::string(queryText));
            const std::vector<BasePair> path(pairs);
            return model.evaluate(target, query, path);
        };
        check(std::abs(model.rt() - 0.616321) < 5e-7,
              "Turner04 RT matches the ViennaRNA constant");
        check(std::abs(evaluate("CG", "CG", {{0U, 1U}, {1U, 0U}}).loops + 2.4) < 1e-9,
              "Turner04 GC/GC stack oracle");
        check(std::abs(evaluate("GCG", "CGC", {{0U, 2U}, {1U, 1U}, {2U, 0U}}).loops + 5.8) < 1e-9,
              "Turner04 three-pair stack oracle");
        check(std::abs(evaluate("AG", "CU", {{0U, 1U}, {1U, 0U}}).loops + 2.1) < 1e-9,
              "Turner04 AU/GC stack orientation oracle");
        check(std::abs(evaluate("GA", "UC", {{0U, 1U}, {1U, 0U}}).loops + 2.4) < 1e-9,
              "Turner04 GC/AU stack orientation oracle");
        check(std::abs(evaluate("GG", "CU", {{0U, 1U}, {1U, 0U}}).loops + 2.1) < 1e-9,
              "Turner04 GU/GC stack orientation oracle");
        check(std::abs(evaluate("GG", "UC", {{0U, 1U}, {1U, 0U}}).loops + 1.5) < 1e-9,
              "Turner04 GC/GU stack orientation oracle");
        check(std::abs(evaluate("GAC", "GC", {{0U, 1U}, {2U, 0U}}).loops - 0.4) < 1e-9,
              "Turner04 one-nucleotide bulge oracle");
        check(std::abs(evaluate("GAC", "GAC", {{0U, 2U}, {2U, 0U}}).loops - 0.8) < 1e-9,
              "Turner04 int11 oracle");
        check(std::abs(evaluate("GAAC", "GAC", {{0U, 2U}, {3U, 0U}}).loops - 2.3) < 1e-9,
              "Turner04 int21 orientation oracle");
        check(std::abs(evaluate("GAC", "GAAC", {{0U, 3U}, {2U, 0U}}).loops - 2.3) < 1e-9,
              "Turner04 reversed int21 orientation oracle");
        check(std::abs(evaluate("UCGU", "AAA", {{0U, 2U}, {3U, 0U}}).loops - 3.7) < 1e-9,
              "Turner04 asymmetric int21 1x2 nucleotide-order oracle");
        check(std::abs(evaluate("AAU", "AAGU", {{0U, 3U}, {2U, 0U}}).loops - 3.7) < 1e-9,
              "Turner04 asymmetric int21 2x1 nucleotide-order oracle");
        check(std::abs(evaluate("GAAC", "GAAC", {{0U, 3U}, {3U, 0U}}).loops - 1.5) < 1e-9,
              "Turner04 int22 oracle");

        const NearestNeighborEnergyModel dangleModel(37.0, "Turner04", true);
        const Sequence target("target", "AAGC");
        const Sequence query("query", "ACUA");
        const std::vector<BasePair> path{{1U, 2U}, {2U, 1U}};
        const auto dangles = dangleModel.evaluate(target, query, path);
        check(std::abs(dangles.loops + 2.1) < 1e-9 &&
              std::abs(dangles.dangleLeft + 1.0) < 1e-9 &&
              std::abs(dangles.dangleRight + 1.1) < 1e-9 &&
              std::abs(dangles.hybrid() - 0.4) < 1e-9,
              "Turner04 exterior mismatch/dangle oracle");

        const NearestNeighborEnergyModel model20(20.0, "Turner04", false);
        const Sequence auTarget("target", "A");
        const Sequence auQuery("query", "U");
        const std::vector<BasePair> auPath{{0U, 0U}};
        const auto at20 = model20.evaluate(auTarget, auQuery, auPath);
        check(std::abs(at20.initiation - 4.07) < 1e-9 &&
              std::abs(at20.endLeft - 0.67) < 1e-9 &&
              std::abs(at20.endRight - 0.67) < 1e-9,
              "ViennaRNA centikcal temperature interpolation oracle");

        const NearestNeighborEnergyModel turner99(37.0, "Turner99", false);
        const NearestNeighborEnergyModel andronescu(37.0, "Andronescu07", false);
        check(std::abs(turner99.initiationEnergy() - 4.1) < 1e-9 &&
              std::abs(andronescu.initiationEnergy() - 4.1) < 1e-9,
              "named ViennaRNA parameter sets load natively");
        check(std::abs(BasePairEnergyModel(0.0).rt() - 1.0) < 1e-12,
              "base-pair-counting energy model uses RT=1");
    }

    {
        const Sequence sequence("iupac", "acgturyswkmbdhvn", -5);
        check(sequence.str() == "ACGUUNNNNNNNNNNN",
              "non-canonical IUPAC symbols normalize to non-pairing N");
        check(sequence.externalIndex(0U) == -5, "external coordinate origin");
        check(!canPair('R', 'Y'), "ambiguous IUPAC symbols do not pair");
        check(!canPair('N', 'A'), "N does not pair");
        check(!canPair('A', 'C'), "incompatible pair rejected");
    }

    if (failureCount == 0) {
        std::cout << "All IntaRNAnew real-interaction integration tests passed.\n";
    }
    return failureCount == 0 ? 0 : 1;
}
