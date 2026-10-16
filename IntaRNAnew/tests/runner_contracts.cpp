#include "intarnanew/runner.hpp"

#include <cmath>
#include <iostream>
#include <string_view>

namespace {

int failures{};

void check(const bool condition, const std::string_view description) {
    if (!condition) {
        std::cerr << "FAIL: " << description << '\n';
        ++failures;
    }
}

[[nodiscard]] auto basePairConfig() -> intarnanew::Config {
    intarnanew::Config config;
    config.energy = intarnanew::EnergyKind::basePair;
    config.mode = intarnanew::PredictionMode::exact;
    config.model = intarnanew::InteractionModel::singleSite;
    config.seed.required = false;
    config.target.accessibility = intarnanew::AccessibilityKind::disabled;
    config.query.accessibility = intarnanew::AccessibilityKind::disabled;
    config.output.number = 20U;
    config.output.overlap = intarnanew::OverlapPolicy::both;
    config.output.deltaEnergy = 100.0;
    return config;
}

} // namespace

auto main() -> int {
    const intarnanew::Sequence target{"target", "CCAACACC"};
    const intarnanew::Sequence query{"query", "GG"};
    auto config = basePairConfig();
    auto evaluated = intarnanew::predictPair(config, target, query);
    check(evaluated.has_value(), "whole-pair evaluation succeeds");
    if (evaluated) {
        check(!evaluated->prediction.interactions.empty(),
              "whole-pair evaluation returns favorable interactions");
        check(evaluated->prediction.rt == 1.0,
              "base-pair runner preserves the model RT contract");
        check(evaluated->targetAccessibility->openingEnergy({0U, 1U}) == 0.0 &&
                  evaluated->queryAccessibility->unpairedProbability({0U, 1U}) == 1.0,
              "runner retains its native accessibility providers");
        check(std::isfinite(evaluated->prediction.targetEnsembleFreeEnergy) &&
                  std::isfinite(evaluated->prediction.queryEnsembleFreeEnergy),
              "runner attaches whole-sequence monomer summaries");
    }

    config.target.regions = "2-7";
    auto restricted = intarnanew::predictPair(config, target, query);
    check(restricted.has_value(), "explicit-region pair evaluation succeeds");
    if (restricted) {
        bool withinRegion = true;
        for (const auto& interaction : restricted->prediction.interactions) {
            const auto range = interaction.targetRange();
            withinRegion = withinRegion && range.begin >= 1U && range.end <= 6U;
        }
        check(withinRegion, "runner applies configured external-coordinate regions");
    }

    config.target.regions = "100-101";
    auto invalid = intarnanew::predictPair(config, target, query);
    check(!invalid && invalid.error().find("target regions") != std::string::npos,
          "runner returns contextual region failures instead of throwing");

    if (failures == 0) {
        std::cout << "All in-process pair-runner contract tests passed.\n";
    }
    return failures == 0 ? 0 : 1;
}
