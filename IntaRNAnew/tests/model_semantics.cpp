#include "intarnanew/accessibility.hpp"
#include "intarnanew/config.hpp"
#include "intarnanew/energy.hpp"
#include "intarnanew/helix_blocks.hpp"
#include "intarnanew/predictor.hpp"
#include "intarnanew/sequence.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

using intarnanew::AccessibilityKind;
using intarnanew::AccessibilityProvider;
using intarnanew::BasePair;
using intarnanew::BasePairEnergyModel;
using intarnanew::Config;
using intarnanew::DisabledAccessibility;
using intarnanew::EnergyKind;
using intarnanew::InteractionModel;
using intarnanew::OverlapPolicy;
using intarnanew::PredictionMode;
using intarnanew::Predictor;
using intarnanew::Sequence;

class SeedOpeningAccessibility final : public AccessibilityProvider {
public:
    [[nodiscard]] auto openingEnergy(const intarnanew::Interval interval) const
        -> intarnanew::Energy override {
        return interval.begin == 0U ? 5.0 : 0.0;
    }
    [[nodiscard]] auto unpairedProbability(const intarnanew::Interval) const
        -> double override { return 1.0; }
    [[nodiscard]] auto positionUnpairedProbability(const intarnanew::Index) const
        -> double override { return 1.0; }
    [[nodiscard]] auto blocked(const intarnanew::Index) const -> bool override { return false; }
};

class NativeProbabilityAccessibility final : public AccessibilityProvider {
public:
    explicit NativeProbabilityAccessibility(const double probability)
        : probability_(probability) {}

    [[nodiscard]] auto openingEnergy(const intarnanew::Interval) const
        -> intarnanew::Energy override { return 0.0; }
    [[nodiscard]] auto unpairedProbability(const intarnanew::Interval) const
        -> double override { return probability_; }
    [[nodiscard]] auto positionUnpairedProbability(const intarnanew::Index) const
        -> double override { return 1.0; }
    [[nodiscard]] auto blocked(const intarnanew::Index) const -> bool override { return false; }

private:
    double probability_{};
};

int failures{};

void check(const bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "FAILED: " << message << '\n';
        ++failures;
    }
}

[[nodiscard]] auto baseConfig() -> Config {
    Config config;
    config.energy = EnergyKind::basePair;
    config.mode = PredictionMode::exact;
    config.model = InteractionModel::seedExtension;
    config.seed.required = false;
    config.target.accessibility = AccessibilityKind::disabled;
    config.query.accessibility = AccessibilityKind::disabled;
    config.output.number = 1'000U;
    config.output.overlap = OverlapPolicy::both;
    config.output.deltaEnergy = 100.0;
    config.output.maxEnergy = 0.0;
    return config;
}

void checkHelixContracts() {
    const Sequence target("target", "CCCCC");
    const Sequence query("query", "GGGGG");
    const DisabledAccessibility targetAccessibility(target);
    const DisabledAccessibility queryAccessibility(query);
    const BasePairEnergyModel energy(37.0);
    intarnanew::HelixConfig helix;
    helix.minBasePairs = 2U;
    helix.maxBasePairs = 2U;

    const std::vector<BasePair> twoPair{{0U, 1U}, {1U, 0U}};
    auto blocks = intarnanew::decomposeHelixBlocks(
        target, query, twoPair, helix, energy,
        targetAccessibility, queryAccessibility);
    check(blocks.size() == 1U && blocks.front().basePairCount() == 2U,
          "a two-pair canonical helix satisfies the default block contract");

    const std::vector<BasePair> overlong{{0U, 3U}, {1U, 2U}, {2U, 1U}, {3U, 0U}};
    check(intarnanew::decomposeHelixBlocks(
              target, query, overlong, helix, energy,
              targetAccessibility, queryAccessibility).empty(),
          "a canonical run cannot be split without a separating interior loop");

    const std::vector<BasePair> blockAndAnchor{{0U, 3U}, {1U, 2U}, {3U, 1U}};
    blocks = intarnanew::decomposeHelixBlocks(
        target, query, blockAndAnchor, helix, energy,
        targetAccessibility, queryAccessibility);
    check(blocks.size() == 1U && blocks.front().lastPair == 1U,
          "a separated right-most initiation pair can terminate a block composition");

    const std::vector<BasePair> twoBlocks{{0U, 4U}, {1U, 3U}, {3U, 2U}, {4U, 1U}};
    blocks = intarnanew::decomposeHelixBlocks(
        target, query, twoBlocks, helix, energy,
        targetAccessibility, queryAccessibility);
    check(blocks.size() == 2U,
          "two bounded helices separated by an interior loop compose an interaction");

    const std::vector<BasePair> oneBaseBulge{{0U, 2U}, {2U, 1U}};
    check(intarnanew::decomposeHelixBlocks(
              target, query, oneBaseBulge, helix, energy,
              targetAccessibility, queryAccessibility).empty(),
          "helixMaxIL zero rejects a one-base bulge inside a helix");
    helix.maxInternalLoop = 1U;
    check(intarnanew::decomposeHelixBlocks(
              target, query, oneBaseBulge, helix, energy,
              targetAccessibility, queryAccessibility).size() == 1U,
          "helixMaxIL bounds the total unpaired bases inside a helix transition");

    helix.maxInternalLoop = 0U;
    helix.maxEnergy = -1.0;
    check(intarnanew::decomposeHelixBlocks(
              target, query, twoPair, helix, energy,
              targetAccessibility, queryAccessibility).empty(),
          "helixMaxE is an exclusive loop-energy threshold");
    helix.useFullEnergy = true;
    check(intarnanew::decomposeHelixBlocks(
              target, query, twoPair, helix, energy,
              targetAccessibility, queryAccessibility).size() == 1U,
          "helixFullE switches the threshold from loop-only to full helix energy");
}

void checkModelPartitions() {
    const Sequence target("target", "CCCC");
    const Sequence query("query", "GGGG");
    const DisabledAccessibility targetAccessibility(target);
    const DisabledAccessibility queryAccessibility(query);

    auto singleSiteConfig = baseConfig();
    singleSiteConfig.model = InteractionModel::singleSite;
    const auto singleSite = Predictor(singleSiteConfig).predict(
        target, query, targetAccessibility, queryAccessibility);
    double expectedPartition{};
    for (const auto& interaction : singleSite.ensembleSites) {
        check(std::abs(interaction.ensembleFreeEnergy - interaction.energy.total()) < 1e-12,
              "model S assigns one MFE Boltzmann weight to each site");
        expectedPartition += std::exp(-interaction.energy.total());
    }
    check(std::abs(std::exp(singleSite.logPartition) - expectedPartition) < 1e-10,
          "model-S global partition is the sum of site-MFE weights");

    auto extensionConfig = baseConfig();
    const auto extension = Predictor(extensionConfig).predict(
        target, query, targetAccessibility, queryAccessibility);
    double extensionPartition{};
    for (const auto& interaction : extension.ensembleSites) {
        check(std::abs(interaction.ensembleFreeEnergy - interaction.energy.total()) < 1e-12,
              "unseeded model X reports one MFE weight per interaction site");
        extensionPartition += std::exp(-interaction.energy.total());
    }
    check(std::abs(std::exp(extension.logPartition) - extensionPartition) < 1e-10,
          "unseeded model X has the same valid site-MFE partition reduction as model S");

    auto seededExtensionConfig = baseConfig();
    seededExtensionConfig.seed.required = true;
    seededExtensionConfig.seed.basePairs = 2U;
    const auto seededExtension = Predictor(seededExtensionConfig).predict(
        target, query, targetAccessibility, queryAccessibility);
    check(!std::isfinite(seededExtension.logPartition) &&
              !std::isfinite(seededExtension.ensembleFreeEnergy),
          "seeded model X leaves its algorithmically incomplete global partition unavailable");
    for (const auto& interaction : seededExtension.ensembleSites) {
        check(std::abs(interaction.ensembleFreeEnergy - interaction.energy.total()) < 1e-12,
              "seeded model X still exposes scientifically valid site-MFE weights");
    }
}

void checkEnsembleCentikcalContract() {
    const Sequence target("target", "GGGAAACCC");
    const Sequence query("query", "GGGAAACCC");
    const DisabledAccessibility targetAccessibility(target);
    const DisabledAccessibility queryAccessibility(query);
    auto config = baseConfig();
    config.model = InteractionModel::ensemble;
    const auto prediction = Predictor(config).predict(
        target, query, targetAccessibility, queryAccessibility);
    check(!prediction.interactions.empty(), "model P produces a best site");
    if (prediction.interactions.empty()) return;
    const auto& best = prediction.interactions.front();
    const auto expectedReported = std::trunc(best.ensembleFreeEnergy * 100.0) / 100.0;
    check(std::abs(best.energy.total() - expectedReported) < 1e-12,
          "model P retains raw site free energy and reports centikcal truncation");
    check(best.ensembleFreeEnergy < best.energy.total(),
          "negative model-P raw energy is more precise than its reported value");
}

void checkEnsembleTerminalWeights() {
    auto config = baseConfig();
    config.energy = EnergyKind::nearestNeighbor;
    config.energyParameters = "Turner04";
    config.model = InteractionModel::ensemble;
    config.output.maxEnergy = 100.0;

    {
        const Sequence target("target", "A");
        const Sequence query("query", "U");
        const DisabledAccessibility targetAccessibility(target);
        const DisabledAccessibility queryAccessibility(query);
        const auto prediction = Predictor(config).predict(
            target, query, targetAccessibility, queryAccessibility);
        check(prediction.ensembleSites.size() == 1U,
              "single AU model-P ensemble contains one interaction site");
        if (!prediction.ensembleSites.empty()) {
            const auto& site = prediction.ensembleSites.front();
            check(std::abs(site.ensembleFreeEnergy - 5.10) < 1e-9,
                  "model-P site partition includes both terminal-AU penalties");
            check(std::abs(prediction.logPartition + 5.10 / prediction.rt) < 1e-9,
                  "model-P global partition includes terminal-AU penalties");
            check(std::abs(site.probability - 1.0) < 1e-12,
                  "single model-P site remains normalized");
        }
    }

    {
        const Sequence target("target", "AAG");
        const Sequence query("query", "AUA");
        const DisabledAccessibility targetAccessibility(target);
        const DisabledAccessibility queryAccessibility(query);
        const std::vector<BasePair> path{{1U, 1U}};
        const intarnanew::NearestNeighborEnergyModel energy(
            37.0, "Turner04", true);
        const auto expected = energy.evaluate(target, query, path).total();
        const auto prediction = Predictor(config).predict(
            target, query, targetAccessibility, queryAccessibility,
            {1U, 1U}, {1U, 1U});
        check(prediction.ensembleSites.size() == 1U,
              "single-cell model-P domain contains one interaction");
        if (!prediction.ensembleSites.empty()) {
            check(std::abs(prediction.ensembleSites.front().ensembleFreeEnergy - expected) < 1e-9,
                  "model-P site partition includes exterior dangle energy");
            check(std::abs(prediction.logPartition + expected / prediction.rt) < 1e-9,
                  "model-P global partition includes exterior dangle energy");
        }
    }
}

void checkSeedOnlyAndHelixPrediction() {
    const Sequence target("target", "CCCC");
    const Sequence query("query", "GGGG");
    const DisabledAccessibility targetAccessibility(target);
    const DisabledAccessibility queryAccessibility(query);

    auto seedConfig = baseConfig();
    seedConfig.mode = PredictionMode::seedOnly;
    seedConfig.seed.required = true;
    seedConfig.seed.basePairs = 2U;
    seedConfig.output.bestSeedOnly = true;
    const auto seeds = Predictor(seedConfig).predict(
        target, query, targetAccessibility, queryAccessibility);
    check(seeds.ensembleSites.size() == 9U,
          "seed-only mode enumerates every exact two-pair seed site");
    for (const auto& interaction : seeds.ensembleSites) {
        check(interaction.pairs.size() == 2U && interaction.seeds.size() == 1U,
              "seed-only interactions contain exactly their seed pairs");
    }

    auto bestSeedConfig = baseConfig();
    bestSeedConfig.seed.required = true;
    bestSeedConfig.seed.basePairs = 2U;
    bestSeedConfig.output.number = 1U;
    bestSeedConfig.output.bestSeedOnly = true;
    const auto bestSeed = Predictor(bestSeedConfig).predict(
        target, query, targetAccessibility, queryAccessibility);
    check(!bestSeed.interactions.empty() &&
              bestSeed.interactions.front().seeds.size() == 3U &&
              bestSeed.interactions.front().seeds[0U].firstPair == 0U &&
              bestSeed.interactions.front().seeds[1U].firstPair == 1U &&
              bestSeed.interactions.front().seeds[2U].firstPair == 2U,
          "discovered seed matches are retained once in deterministic energy/coordinate order");
    check(!bestSeed.interactions.empty() && bestSeed.interactions.front().bestSeed() != nullptr &&
              bestSeed.interactions.front().bestSeed()->firstPair == 0U &&
              bestSeed.interactions.front().bestSeed()->lastPair == 1U,
          "best-seed helper resolves equal-energy seeds by their interaction coordinates");

    auto explicitConfig = baseConfig();
    explicitConfig.seed.explicitSeeds = "1||&3||,3||&1||";
    explicitConfig.output.maxEnergy = 999.0;
    const SeedOpeningAccessibility seedOpeningAccessibility;
    const auto explicitPrediction = Predictor(explicitConfig).predict(
        target, query, seedOpeningAccessibility, queryAccessibility);
    bool explicitFullSiteChecked = false;
    for (const auto& interaction : explicitPrediction.ensembleSites) {
        if (interaction.targetRange() != intarnanew::Interval{0U, 3U} ||
            interaction.queryRange() != intarnanew::Interval{0U, 3U}) continue;
        explicitFullSiteChecked = true;
        check(interaction.seeds.size() == 2U && interaction.bestSeed() != nullptr &&
                  interaction.bestSeed()->firstPair == 2U &&
                  interaction.bestSeed()->lastPair == 3U,
              "best-seed selection evaluates every matching explicit seed, not input order");
    }
    check(explicitFullSiteChecked, "explicit-seed regression contains its full interaction site");

    auto helixConfig = baseConfig();
    helixConfig.mode = PredictionMode::heuristic;
    helixConfig.model = InteractionModel::helixBlocks;
    helixConfig.helix.maxBasePairs = 2U;
    const auto helixPrediction = Predictor(helixConfig).predict(
        target, query, targetAccessibility, queryAccessibility);
    check(!helixPrediction.interactions.empty() &&
              helixPrediction.interactions.front().pairs.size() == 3U &&
              std::abs(helixPrediction.interactions.front().energy.total() + 3.0) < 1e-12,
          "model B composes a bounded two-pair block with a separated anchor");
}

void checkNativeInteractionProbabilities() {
    const Sequence target("target", "CCCC");
    const Sequence query("query", "GGGG");
    const NativeProbabilityAccessibility targetAccessibility(0.37);
    const NativeProbabilityAccessibility queryAccessibility(0.61);

    auto config = baseConfig();
    for (const auto mode : {PredictionMode::exact, PredictionMode::seedOnly}) {
        config.mode = mode;
        config.seed.required = mode == PredictionMode::seedOnly;
        config.seed.basePairs = 2U;
        const auto prediction = Predictor(config).predict(
            target, query, targetAccessibility, queryAccessibility);
        check(!prediction.ensembleSites.empty(),
              "native-Pu regression produces interaction sites");
        for (const auto& interaction : prediction.ensembleSites) {
            check(std::abs(interaction.unpairedTarget - 0.37) < 1e-12 &&
                      std::abs(interaction.unpairedQuery - 0.61) < 1e-12,
                  "interaction preserves provider-native Pu independently of energy RT");
        }
    }

    config.mode = PredictionMode::exact;
    config.model = InteractionModel::ensemble;
    config.seed.required = false;
    const auto ensemble = Predictor(config).predict(
        target, query, targetAccessibility, queryAccessibility);
    check(!ensemble.ensembleSites.empty(), "model-P native-Pu regression produces sites");
    for (const auto& interaction : ensemble.ensembleSites) {
        check(std::abs(interaction.unpairedTarget - 0.37) < 1e-12 &&
                  std::abs(interaction.unpairedQuery - 0.61) < 1e-12,
              "model P preserves provider-native Pu during representative reduction");
    }
}

} // namespace

auto main() -> int {
    checkHelixContracts();
    checkModelPartitions();
    checkEnsembleCentikcalContract();
    checkEnsembleTerminalWeights();
    checkSeedOnlyAndHelixPrediction();
    checkNativeInteractionProbabilities();
    if (failures == 0) {
        std::cout << "All interaction-model semantic tests passed.\n";
    }
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
