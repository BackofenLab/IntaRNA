#include "intarnanew/accessibility.hpp"
#include "intarnanew/config.hpp"
#include "intarnanew/predictor.hpp"
#include "intarnanew/sequence.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <span>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace {

using intarnanew::BasePair;
using intarnanew::Config;
using intarnanew::DisabledAccessibility;
using intarnanew::Energy;
using intarnanew::Index;
using intarnanew::Interaction;
using intarnanew::InteractionModel;
using intarnanew::PredictionMode;
using intarnanew::PredictionResult;
using intarnanew::Predictor;
using intarnanew::Sequence;
using intarnanew::canPair;
using intarnanew::infinity;

struct SiteKey {
    Index targetBegin{};
    Index targetEnd{};
    Index queryBegin{};
    Index queryEnd{};

    friend auto operator<=>(const SiteKey&, const SiteKey&) = default;
};

struct OracleSite {
    Energy minimumEnergy{infinity};
    double partition{};
    std::size_t structureCount{};
};

using OracleSites = std::map<SiteKey, OracleSite>;

int failures{};

void check(const bool condition, const std::string_view description) {
    if (!condition) {
        std::cerr << "FAILED: " << description << '\n';
        ++failures;
    }
}

[[nodiscard]] auto close(
    const double left,
    const double right,
    const double relativeTolerance = 1e-11) -> bool {
    const auto scale = std::max({1.0, std::abs(left), std::abs(right)});
    return std::abs(left - right) <= relativeTolerance * scale;
}

[[nodiscard]] auto isStack(const BasePair left, const BasePair right) noexcept -> bool {
    return right.target == left.target + 1U && left.query == right.query + 1U;
}

[[nodiscard]] auto noLonelyPairs(const std::span<const BasePair> path) noexcept -> bool {
    if (path.size() < 2U) return false;
    for (Index index = 0U; index < path.size(); ++index) {
        const bool leftStack = index > 0U && isStack(path[index - 1U], path[index]);
        const bool rightStack = index + 1U < path.size() && isStack(path[index], path[index + 1U]);
        if (!leftStack && !rightStack) return false;
    }
    return true;
}

[[nodiscard]] auto siteKey(const std::span<const BasePair> path) -> SiteKey {
    return {
        path.front().target,
        path.back().target,
        path.back().query,
        path.front().query,
    };
}

// This enumerator intentionally does not reuse any production recurrence or
// energy implementation. For the documented base-pair backend each legal
// intermolecular pair contributes exactly -1 and RT is exactly 1.
[[nodiscard]] auto enumerate(
    const Sequence& target,
    const Sequence& query,
    const Config& config) -> OracleSites {
    OracleSites sites;
    std::vector<BasePair> path;

    const auto registerPath = [&] {
        if (config.output.noLonelyPairs && !noLonelyPairs(path)) return;
        const auto energy = -static_cast<Energy>(path.size());
        auto& site = sites[siteKey(path)];
        site.minimumEnergy = std::min(site.minimumEnergy, energy);
        site.partition += std::exp(-energy); // RT=1 for --energy=B
        ++site.structureCount;
    };

    std::function<void(BasePair)> extend;
    extend = [&](const BasePair previous) {
        registerPath();
        const auto targetLast = std::min(
            target.size(), previous.target + config.target.interactionLoopMax + 2U);
        const auto queryDistanceMax = std::min(
            previous.query, config.query.interactionLoopMax + 1U);
        for (Index nextTarget = previous.target + 1U; nextTarget < targetLast; ++nextTarget) {
            for (Index queryDistance = 1U; queryDistance <= queryDistanceMax; ++queryDistance) {
                const Index nextQuery = previous.query - queryDistance;
                if (!canPair(target[nextTarget], query[nextQuery])) continue;
                if (config.target.interactionLengthMax != 0U &&
                    nextTarget - path.front().target + 1U > config.target.interactionLengthMax) {
                    continue;
                }
                if (config.query.interactionLengthMax != 0U &&
                    path.front().query - nextQuery + 1U > config.query.interactionLengthMax) {
                    continue;
                }
                path.push_back({nextTarget, nextQuery});
                extend(path.back());
                path.pop_back();
            }
        }
    };

    for (Index targetIndex = 0U; targetIndex < target.size(); ++targetIndex) {
        for (Index queryIndex = 0U; queryIndex < query.size(); ++queryIndex) {
            if (!canPair(target[targetIndex], query[queryIndex])) continue;
            path.push_back({targetIndex, queryIndex});
            extend(path.back());
            path.pop_back();
        }
    }
    return sites;
}

[[nodiscard]] auto predictionSites(const PredictionResult& result)
    -> std::map<SiteKey, const Interaction*> {
    std::map<SiteKey, const Interaction*> sites;
    for (const auto& interaction : result.ensembleSites) {
        const auto targetRange = interaction.targetRange();
        const auto queryRange = interaction.queryRange();
        const auto [iterator, inserted] = sites.emplace(
            SiteKey{targetRange.begin, targetRange.end, queryRange.begin, queryRange.end},
            &interaction);
        static_cast<void>(iterator);
        check(inserted, "predictor emits each interaction site once");
    }
    return sites;
}

void compareCase(
    const std::string_view targetBases,
    const std::string_view queryBases,
    const Index targetLoopMax,
    const Index queryLoopMax,
    const bool rejectLonelyPairs,
    const Index targetLengthMax = 0U,
    const Index queryLengthMax = 0U,
    const InteractionModel model = InteractionModel::singleSite) {
    Config config;
    config.mode = PredictionMode::exact;
    config.model = model;
    config.energy = intarnanew::EnergyKind::basePair;
    config.seed.required = false;
    config.target.accessibility = intarnanew::AccessibilityKind::disabled;
    config.query.accessibility = intarnanew::AccessibilityKind::disabled;
    config.target.interactionLoopMax = targetLoopMax;
    config.query.interactionLoopMax = queryLoopMax;
    config.target.interactionLengthMax = targetLengthMax;
    config.query.interactionLengthMax = queryLengthMax;
    config.output.noLonelyPairs = rejectLonelyPairs;
    config.output.number = 1'000U;
    config.output.deltaEnergy = 100.0;
    config.output.overlap = intarnanew::OverlapPolicy::both;
    config.output.maxEnergy = 0.0;

    const Sequence target("target", std::string(targetBases));
    const Sequence query("query", std::string(queryBases));
    const DisabledAccessibility targetAccessibility(target);
    const DisabledAccessibility queryAccessibility(query);

    const auto oracle = enumerate(target, query, config);
    const auto prediction = Predictor(config).predict(
        target, query, targetAccessibility, queryAccessibility);
    const auto observed = predictionSites(prediction);

    check(observed.size() == oracle.size(), "exact predictor finds every enumerated site");
    double expectedPartition{};
    for (const auto& [key, expected] : oracle) {
        const auto expectedSitePartition = model == InteractionModel::ensemble
            ? expected.partition
            : std::exp(-expected.minimumEnergy);
        expectedPartition += expectedSitePartition;
        const auto found = observed.find(key);
        check(found != observed.end(), "enumerated site is present in predictor result");
        if (found == observed.end()) continue;
        const auto& interaction = *found->second;
        if (model == InteractionModel::ensemble) {
            const auto expectedFreeEnergy = -std::log(expected.partition);
            const auto expectedReportedEnergy = std::trunc(expectedFreeEnergy * 100.0) / 100.0;
            check(close(interaction.energy.total(), expectedReportedEnergy),
                  "ensemble-model reported energy is centikcal-truncated site free energy");
        } else {
            check(close(interaction.energy.total(), expected.minimumEnergy),
                  "single-site model energy matches independent MFE enumeration");
        }
        const auto observedSitePartition = std::exp(-interaction.ensembleFreeEnergy);
        check(close(observedSitePartition, expectedSitePartition),
              model == InteractionModel::ensemble
                  ? "ensemble-model site partition matches independent path enumeration"
                  : "single-site model contributes exactly one MFE weight per site");
    }

    if (expectedPartition == 0.0) {
        check(!std::isfinite(prediction.logPartition), "empty oracle has empty predictor partition");
        return;
    }
    check(close(std::exp(prediction.logPartition), expectedPartition),
          "global partition matches independent path enumeration");
    check(close(prediction.ensembleFreeEnergy, -std::log(expectedPartition)),
          "global ensemble energy is -RT log(Z) with RT=1");
    for (const auto& [key, expected] : oracle) {
        const auto found = observed.find(key);
        if (found == observed.end()) continue;
        const auto expectedSitePartition = model == InteractionModel::ensemble
            ? expected.partition
            : std::exp(-expected.minimumEnergy);
        check(close(found->second->probability, expectedSitePartition / expectedPartition),
              "site probability is normalized by the complete partition");
    }
}
void compareCompactKernel() {
    std::uint32_t randomState{0x6d2b79f5U};
    const auto nextRandom = [&]() -> std::uint32_t {
        randomState ^= randomState << 13U;
        randomState ^= randomState >> 17U;
        randomState ^= randomState << 5U;
        return randomState;
    };
    const auto randomSequence = [&](const Index length) -> std::string {
        static constexpr std::string_view alphabet{"ACGU"};
        std::string result(length, 'A');
        for (auto& nucleotide : result) {
            nucleotide = alphabet[nextRandom() % alphabet.size()];
        }
        return result;
    };
    const auto selectedSites = [](const PredictionResult& result) {
        std::vector<std::tuple<SiteKey, Energy>> sites;
        sites.reserve(result.interactions.size());
        for (const auto& interaction : result.interactions) {
            const auto targetRange = interaction.targetRange();
            const auto queryRange = interaction.queryRange();
            sites.emplace_back(
                SiteKey{
                    targetRange.begin, targetRange.end,
                    queryRange.begin, queryRange.end,
                },
                interaction.energy.total());
        }
        return sites;
    };

    for (Index caseIndex{}; caseIndex < 300U; ++caseIndex) {
        const Index targetLength = 3U + nextRandom() % 7U;
        const Index queryLength = 3U + nextRandom() % 7U;
        const Sequence target("target", randomSequence(targetLength));
        const Sequence query("query", randomSequence(queryLength));
        const DisabledAccessibility targetAccessibility(target);
        const DisabledAccessibility queryAccessibility(query);

        Config complete;
        complete.mode = PredictionMode::exact;
        complete.model = InteractionModel::singleSite;
        complete.energy = intarnanew::EnergyKind::basePair;
        complete.seed.required = false;
        complete.target.accessibility = intarnanew::AccessibilityKind::disabled;
        complete.query.accessibility = intarnanew::AccessibilityKind::disabled;
        complete.target.interactionLoopMax = nextRandom() % 5U;
        complete.query.interactionLoopMax = nextRandom() % 5U;
        if (caseIndex % 17U == 0U) {
            complete.target.interactionLengthMax = 0U;
            complete.query.interactionLengthMax = 0U;
        } else if (caseIndex % 19U == 0U) {
            complete.target.interactionLengthMax = 70'000U;
            complete.query.interactionLengthMax = 70'000U;
        } else {
            complete.target.interactionLengthMax =
                1U + nextRandom() % std::min<Index>(6U, targetLength);
            complete.query.interactionLengthMax =
                1U + nextRandom() % std::min<Index>(6U, queryLength);
        }
        complete.output.number = 1U + nextRandom() % 6U;
        complete.output.overlap = intarnanew::OverlapPolicy::both;
        complete.output.deltaEnergy = static_cast<Energy>(nextRandom() % 5U);
        complete.output.maxEnergy = -static_cast<Energy>(nextRandom() % 3U);
        complete.additiveEnergy =
            static_cast<Energy>(static_cast<int>(nextRandom() % 9U) - 4) * 0.25;

        auto compact = complete;
        compact.predictionRequirements.retainAllSites = false;
        compact.predictionRequirements.computeInteractionPartition = false;
        compact.predictionRequirements.traceback = false;

        const auto expected = Predictor(complete).predict(
            target, query, targetAccessibility, queryAccessibility);
        const auto observed = Predictor(compact).predict(
            target, query, targetAccessibility, queryAccessibility);
        const auto expectedSites = selectedSites(expected);
        const auto observedSites = selectedSites(observed);
        check(observedSites.size() == expectedSites.size(),
              "compact predictor returns the same selected-site count");
        const auto common = std::min(observedSites.size(), expectedSites.size());
        for (Index index{}; index < common; ++index) {
            check(std::get<0>(observedSites[index]) == std::get<0>(expectedSites[index]),
                  "compact predictor preserves selected-site ordering and boundaries");
            check(close(std::get<1>(observedSites[index]), std::get<1>(expectedSites[index])),
                  "compact predictor preserves selected-site energy");
        }
        check(observed.ensembleSites.empty(),
              "selected-only compact prediction does not materialize the site ensemble");
        check(!std::isfinite(observed.logPartition),
              "selected-only compact prediction skips the unrequested partition");
    }
    const auto completeConfig = [] {
        Config config;
        config.mode = PredictionMode::exact;
        config.model = InteractionModel::singleSite;
        config.energy = intarnanew::EnergyKind::basePair;
        config.seed.required = false;
        config.target.accessibility = intarnanew::AccessibilityKind::disabled;
        config.query.accessibility = intarnanew::AccessibilityKind::disabled;
        config.target.interactionLengthMax = 0U;
        config.query.interactionLengthMax = 0U;
        config.output.number = 20U;
        config.output.overlap = intarnanew::OverlapPolicy::both;
        config.output.deltaEnergy = 100.0;
        return config;
    };
    const auto selectedOnly = [](Config config) {
        config.predictionRequirements.retainAllSites = false;
        config.predictionRequirements.computeInteractionPartition = false;
        config.predictionRequirements.traceback = false;
        return config;
    };
    const auto sameSelection = [&](const PredictionResult& expected,
                                   const PredictionResult& observed) {
        const auto expectedSites = selectedSites(expected);
        const auto observedSites = selectedSites(observed);
        if (expectedSites.size() != observedSites.size()) return false;
        for (Index index{}; index < expectedSites.size(); ++index) {
            if (std::get<0>(expectedSites[index]) != std::get<0>(observedSites[index]) ||
                !close(std::get<1>(expectedSites[index]), std::get<1>(observedSites[index]))) {
                return false;
            }
        }
        return true;
    };

    {
        const Sequence target("target", "CCCC");
        const Sequence query("query", "GGGG");
        const DisabledAccessibility targetAccessibility(target, "..b.");
        const DisabledAccessibility queryAccessibility(query);
        const auto complete = completeConfig();
        const auto compact = selectedOnly(complete);
        const auto domain = intarnanew::Interval{0U, 3U};
        const auto expected = Predictor(complete).predict(
            target, query, targetAccessibility, queryAccessibility, domain, domain);
        const auto observed = Predictor(compact).predict(
            target, query, targetAccessibility, queryAccessibility, domain, domain);
        check(sameSelection(expected, observed),
              "compact eligibility honors the actual accessibility providers");
        check(!observed.ensembleSites.empty(),
              "blocked accessibility safely falls back to the full predictor");
    }

    {
        const Sequence target("target", "CCCCCCC");
        const Sequence query("query", "GGG");
        const DisabledAccessibility targetAccessibility(target);
        const DisabledAccessibility queryAccessibility(query);
        auto complete = completeConfig();
        complete.target.regions = "1-3,5-7";
        const auto compact = selectedOnly(complete);
        const auto expected = Predictor(complete).predict(
            target, query, targetAccessibility, queryAccessibility);
        const auto observed = Predictor(compact).predict(
            target, query, targetAccessibility, queryAccessibility);
        check(sameSelection(expected, observed),
              "selected-only requirements preserve multi-region reduction");
        check(!observed.ensembleSites.empty(),
              "multi-region prediction safely retains reducible site results");
    }

    {
        const Sequence target("target", "CCCC");
        const Sequence query("query", "GGGG");
        const DisabledAccessibility targetAccessibility(target);
        const DisabledAccessibility queryAccessibility(query);
        auto complete = completeConfig();
        complete.target.interactionLengthMax = 70'000U;
        complete.query.interactionLengthMax = 70'000U;
        complete.additiveEnergy = 1e300;
        complete.output.maxEnergy = infinity;
        complete.output.deltaEnergy = 0.0;
        complete.output.number = std::numeric_limits<Index>::max();
        const auto compact = selectedOnly(complete);
        const auto expected = Predictor(complete).predict(
            target, query, targetAccessibility, queryAccessibility);
        const auto observed = Predictor(compact).predict(
            target, query, targetAccessibility, queryAccessibility);
        check(sameSelection(expected, observed),
              "compact ranking matches generic floating-energy tie semantics");
        check(observed.ensembleSites.empty(),
              "actual finite domains enable compact prediction for oversized limits");
    }
}

} // namespace

auto main() -> int {
    compareCase("CCAC", "GGUG", 0U, 0U, false);
    compareCase("CCAUC", "GGAUG", 2U, 2U, false);
    compareCase("CCACCC", "GGUGG", 4U, 3U, false);
    compareCase("CCAUC", "GGAUG", 3U, 3U, true);
    compareCase("CCACCC", "GGUGG", 4U, 4U, true);
    compareCase("CCACCC", "GGUGG", 4U, 4U, false, 3U, 3U);

    compareCase("CCAC", "GGUG", 0U, 0U, false, 0U, 0U, InteractionModel::ensemble);
    compareCase("CCAUC", "GGAUG", 2U, 2U, false, 0U, 0U, InteractionModel::ensemble);
    compareCase("CCACCC", "GGUGG", 4U, 4U, true, 0U, 0U, InteractionModel::ensemble);

    compareCompactKernel();
    if (failures == 0) {
        std::cout << "All independent predictor-oracle tests passed.\n";
    }
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
