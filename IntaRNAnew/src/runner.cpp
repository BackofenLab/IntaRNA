#include "intarnanew/runner.hpp"

#include <exception>
#include <utility>
#include <vector>

namespace intarnanew {

auto predictPair(
    const Config& config,
    const Sequence& target,
    const Sequence& query) -> std::expected<PairPrediction, std::string> {
    try {
        auto targetAccessibility = makeAccessibility(target, config.target, config);
        if (!targetAccessibility) {
            return std::unexpected("target accessibility: " + targetAccessibility.error());
        }
        auto queryAccessibility = makeAccessibility(query, config.query, config);
        if (!queryAccessibility) {
            return std::unexpected("query accessibility: " + queryAccessibility.error());
        }

        auto targetRegions = configuredRegions(target, config.target);
        if (!targetRegions) return std::unexpected("target regions: " + targetRegions.error());
        auto queryRegions = configuredRegions(query, config.query);
        if (!queryRegions) return std::unexpected("query regions: " + queryRegions.error());
        auto plannedTarget = decomposeAccessibleRegions(
            *targetRegions, config.target.regionLengthMax, config.seed.basePairs,
            **targetAccessibility);
        auto plannedQuery = decomposeAccessibleRegions(
            *queryRegions, config.query.regionLengthMax, config.seed.basePairs,
            **queryAccessibility);

        Predictor predictor(config);
        std::vector<PredictionResult> domains;
        for (const auto targetRegion : plannedTarget) {
            const auto targetWindows = decomposeWindows(
                targetRegion, config.windowWidth, config.windowOverlap);
            for (const auto queryRegion : plannedQuery) {
                const auto queryWindows = decomposeWindows(
                    queryRegion, config.windowWidth, config.windowOverlap);
                for (const auto targetWindow : targetWindows) {
                    for (const auto queryWindow : queryWindows) {
                        domains.push_back(predictor.predict(
                            target, query, **targetAccessibility, **queryAccessibility,
                            targetWindow, queryWindow));
                    }
                }
            }
        }

        PredictionResult result;
        if (domains.empty()) {
            result.rt = predictor.rt();
        } else if (domains.size() == 1U) {
            result = std::move(domains.front());
        } else {
            result = reducePredictions(domains, config);
        }
        if (const auto value = (*targetAccessibility)->ensembleLogPartition()) {
            result.targetLogPartition = *value;
        }
        if (const auto value = (*queryAccessibility)->ensembleLogPartition()) {
            result.queryLogPartition = *value;
        }
        if (const auto value = (*targetAccessibility)->ensembleFreeEnergy()) {
            result.targetEnsembleFreeEnergy = *value;
        }
        if (const auto value = (*queryAccessibility)->ensembleFreeEnergy()) {
            result.queryEnsembleFreeEnergy = *value;
        }
        return PairPrediction{
            std::move(result), std::move(*targetAccessibility), std::move(*queryAccessibility)};
    } catch (const std::exception& error) {
        return std::unexpected(error.what());
    } catch (...) {
        return std::unexpected("unknown prediction failure");
    }
}

} // namespace intarnanew
