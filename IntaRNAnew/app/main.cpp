#include "intarnanew/accessibility.hpp"
#include "intarnanew/cli.hpp"
#include "intarnanew/output.hpp"
#include "intarnanew/output_plan.hpp"
#include "intarnanew/parallel.hpp"
#include "intarnanew/predictor.hpp"
#include "intarnanew/sequence.hpp"

#include <algorithm>
#include <cctype>
#include <exception>
#include <expected>
#include <iostream>
#include <limits>
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

struct OutputGroup {
    std::size_t targetIndex{};
    std::size_t queryIndex{};
    std::size_t targetRegionIndex{};
    std::size_t queryRegionIndex{};
    std::vector<std::size_t> workIndices;
};

struct WorkItem {
    std::size_t targetIndex{};
    std::size_t queryIndex{};
    intarnanew::Interval targetDomain;
    intarnanew::Interval queryDomain;
};

struct WorkResult {
    intarnanew::PredictionResult prediction;
};

[[nodiscard]] auto lowerAscii(const std::string_view text) -> std::string {
    std::string result(text);
    std::ranges::transform(result, result.begin(), [](const unsigned char character) {
        return static_cast<char>(std::tolower(character));
    });
    return result;
}

[[nodiscard]] auto outputPrefix(const std::string_view descriptor) -> std::string {
    return lowerAscii(descriptor.substr(0U, descriptor.find(':')));
}

[[nodiscard]] auto needsCompleteSiteEnsemble(const intarnanew::Config& config) -> bool {
    return std::ranges::any_of(config.output.destinations, [](const std::string& descriptor) {
        const auto prefix = outputPrefix(descriptor);
        return prefix == "qmine" || prefix == "tmine" || prefix == "qspotprob" ||
               prefix == "tspotprob" || prefix == "pmine" || prefix == "spotprob";
    });
}

[[nodiscard]] auto csvHasColumn(
    const std::string_view specification,
    const std::string_view requested) -> bool {
    if (specification.empty() || specification == "*") return true;
    std::size_t begin{};
    while (begin <= specification.size()) {
        const auto comma = specification.find(',', begin);
        const auto end = comma == std::string_view::npos ? specification.size() : comma;
        auto column = specification.substr(begin, end - begin);
        while (!column.empty() && std::isspace(static_cast<unsigned char>(column.front())) != 0) {
            column.remove_prefix(1U);
        }
        while (!column.empty() && std::isspace(static_cast<unsigned char>(column.back())) != 0) {
            column.remove_suffix(1U);
        }
        if (lowerAscii(column) == lowerAscii(requested)) return true;
        if (comma == std::string_view::npos) break;
        begin = comma + 1U;
    }
    return false;
}

[[nodiscard]] auto needsInteractionPartition(const intarnanew::Config& config) -> bool {
    if (config.output.mode == intarnanew::OutputMode::ensemble) return true;
    if (config.output.mode == intarnanew::OutputMode::csv &&
        (csvHasColumn(config.output.csvColumns, "Eall") ||
         csvHasColumn(config.output.csvColumns, "Zall") ||
         csvHasColumn(config.output.csvColumns, "EallTotal") ||
         csvHasColumn(config.output.csvColumns, "P_E"))) {
        return true;
    }
    return std::ranges::any_of(config.output.destinations, [](const std::string& descriptor) {
        const auto prefix = outputPrefix(descriptor);
        return prefix == "qspotprob" || prefix == "tspotprob" || prefix == "spotprob";
    });
}

[[nodiscard]] auto needsMonomerPartition(const intarnanew::Config& config) -> bool {
    if (config.output.mode == intarnanew::OutputMode::ensemble) return true;
    return config.output.mode == intarnanew::OutputMode::csv &&
        (csvHasColumn(config.output.csvColumns, "Etotal") ||
         csvHasColumn(config.output.csvColumns, "Eall1") ||
         csvHasColumn(config.output.csvColumns, "Eall2") ||
         csvHasColumn(config.output.csvColumns, "Zall1") ||
         csvHasColumn(config.output.csvColumns, "Zall2") ||
         csvHasColumn(config.output.csvColumns, "EallTotal"));
}

[[nodiscard]] auto needsInteractionTraceback(const intarnanew::Config& config) -> bool {
    if (config.output.mode == intarnanew::OutputMode::normal ||
        config.output.mode == intarnanew::OutputMode::detailed) {
        return true;
    }
    if (config.output.mode != intarnanew::OutputMode::csv) return false;
    static constexpr std::string_view structuralColumns[]{
        "hybridDP", "hybridDB", "hybridDPfull", "hybridDBfull", "bpList",
        "seedStart1", "seedEnd1", "seedStart2", "seedEnd2",
        "seedE", "seedED1", "seedED2", "seedPu1", "seedPu2",
    };
    return std::ranges::any_of(structuralColumns, [&](const std::string_view column) {
        return csvHasColumn(config.output.csvColumns, column);
    });
}
} // namespace

auto main(const int argc, char** argv) -> int {
    try {
        std::vector<std::string_view> arguments;
        arguments.reserve(static_cast<std::size_t>(std::max(0, argc - 1)));
        for (int index = 1; index < argc; ++index) arguments.emplace_back(argv[index]);

        auto configResult = intarnanew::Cli::parse(arguments, argv[0]);
        if (!configResult) {
            std::cerr << "IntaRNAnew: " << configResult.error() << "\nRun --help for usage.\n";
            return 2;
        }
        auto config = std::move(*configResult);
        if (config.action == intarnanew::RunAction::help) {
            std::cout << intarnanew::Cli::help(false);
            return 0;
        }
        if (config.action == intarnanew::RunAction::fullHelp) {
            std::cout << intarnanew::Cli::help(true);
            return 0;
        }
        if (config.action == intarnanew::RunAction::version) {
            std::cout << intarnanew::Cli::version();
            return 0;
        }

        config.predictionRequirements.computeMonomerPartition =
            needsMonomerPartition(config);

        auto targetResult = intarnanew::SequenceReader::read(
            config.target.input, config.target.id, config.target.firstPosition, std::cin);
        if (!targetResult) {
            std::cerr << "IntaRNAnew: target input: " << targetResult.error() << '\n';
            return 2;
        }
        auto queryResult = intarnanew::SequenceReader::read(
            config.query.input, config.query.id, config.query.firstPosition, std::cin);
        if (!queryResult) {
            std::cerr << "IntaRNAnew: query input: " << queryResult.error() << '\n';
            return 2;
        }
        auto targets = std::move(*targetResult);
        auto queries = std::move(*queryResult);
        auto selectedTargets = intarnanew::SequenceReader::select(
            std::move(targets), config.target.subset);
        if (!selectedTargets) {
            std::cerr << "IntaRNAnew: target set: " << selectedTargets.error() << '\n';
            return 2;
        }
        auto selectedQueries = intarnanew::SequenceReader::select(
            std::move(queries), config.query.subset);
        if (!selectedQueries) {
            std::cerr << "IntaRNAnew: query set: " << selectedQueries.error() << '\n';
            return 2;
        }
        targets = std::move(*selectedTargets);
        queries = std::move(*selectedQueries);
        if (config.output.pairwise && targets.size() != queries.size()) {
            std::cerr << "IntaRNAnew: --outPairwise requires equal target and query record counts\n";
            return 2;
        }

        std::vector<std::unique_ptr<intarnanew::AccessibilityProvider>> targetAccessibility;
        targetAccessibility.reserve(targets.size());
        for (std::size_t index = 0U; index < targets.size(); ++index) {
            auto accessibility = intarnanew::makeAccessibility(
                targets[index], config.target, config);
            if (!accessibility) {
                std::cerr << "IntaRNAnew: target accessibility " << index + 1U << ": "
                          << accessibility.error() << '\n';
                return 2;
            }
            targetAccessibility.push_back(std::move(*accessibility));
        }
        std::vector<std::unique_ptr<intarnanew::AccessibilityProvider>> queryAccessibility;
        queryAccessibility.reserve(queries.size());
        for (std::size_t index = 0U; index < queries.size(); ++index) {
            auto accessibility = intarnanew::makeAccessibility(
                queries[index], config.query, config);
            if (!accessibility) {
                std::cerr << "IntaRNAnew: query accessibility " << index + 1U << ": "
                          << accessibility.error() << '\n';
                return 2;
            }
            queryAccessibility.push_back(std::move(*accessibility));
        }

        const auto planRegions = [&](const std::vector<intarnanew::Sequence>& sequences,
                                     const intarnanew::SideConfig& side,
                                     const auto& accessibility,
                                     const std::string_view label)
            -> std::expected<std::vector<std::vector<intarnanew::Interval>>, std::string> {
            std::vector<std::vector<intarnanew::Interval>> all;
            all.reserve(sequences.size());
            for (std::size_t index = 0U; index < sequences.size(); ++index) {
                auto configured = intarnanew::configuredRegions(sequences[index], side);
                if (!configured) {
                    return std::unexpected(std::string(label) + " " + std::to_string(index + 1U) +
                                           " regions: " + configured.error());
                }
                try {
                    all.push_back(intarnanew::decomposeAccessibleRegions(
                        *configured, side.regionLengthMax, config.seed.basePairs,
                        *accessibility[index]));
                } catch (const std::exception& exception) {
                    return std::unexpected(std::string(label) + " " + std::to_string(index + 1U) +
                                           " automatic regions: " + exception.what());
                }
            }
            return all;
        };
        auto targetRegionsResult = planRegions(
            targets, config.target, targetAccessibility, "target");
        if (!targetRegionsResult) {
            std::cerr << "IntaRNAnew: " << targetRegionsResult.error() << '\n';
            return 2;
        }
        auto queryRegionsResult = planRegions(
            queries, config.query, queryAccessibility, "query");
        if (!queryRegionsResult) {
            std::cerr << "IntaRNAnew: " << queryRegionsResult.error() << '\n';
            return 2;
        }
        const auto& targetRegions = *targetRegionsResult;
        const auto& queryRegions = *queryRegionsResult;

        std::vector<std::pair<std::size_t, std::size_t>> pairs;
        if (config.output.pairwise) {
            for (std::size_t index = 0; index < targets.size(); ++index) pairs.emplace_back(index, index);
        } else {
            for (std::size_t targetIndex = 0; targetIndex < targets.size(); ++targetIndex) {
                for (std::size_t queryIndex = 0; queryIndex < queries.size(); ++queryIndex) {
                    pairs.emplace_back(targetIndex, queryIndex);
                }
            }
        }

        std::vector<OutputGroup> groups;
        std::vector<WorkItem> workItems;
        const auto appendDomains = [&](const std::size_t groupIndex,
                                       const std::size_t targetIndex,
                                       const std::size_t queryIndex,
                                       const intarnanew::Interval targetRegion,
                                       const intarnanew::Interval queryRegion) {
            const auto targetWindows = intarnanew::decomposeWindows(
                targetRegion, config.windowWidth, config.windowOverlap);
            const auto queryWindows = intarnanew::decomposeWindows(
                queryRegion, config.windowWidth, config.windowOverlap);
            for (const auto targetWindow : targetWindows) {
                for (const auto queryWindow : queryWindows) {
                    const auto workIndex = workItems.size();
                    workItems.push_back({targetIndex, queryIndex,
                                         targetWindow, queryWindow});
                    groups[groupIndex].workIndices.push_back(workIndex);
                }
            }
        };
        for (const auto& [targetIndex, queryIndex] : pairs) {
            if (!config.output.perRegion) {
                const auto groupIndex = groups.size();
                groups.push_back({targetIndex, queryIndex, 0U, 0U, {}});
                for (std::size_t targetRegionIndex = 0U;
                     targetRegionIndex < targetRegions[targetIndex].size(); ++targetRegionIndex) {
                    for (std::size_t queryRegionIndex = 0U;
                         queryRegionIndex < queryRegions[queryIndex].size(); ++queryRegionIndex) {
                        appendDomains(groupIndex, targetIndex, queryIndex,
                            targetRegions[targetIndex][targetRegionIndex],
                            queryRegions[queryIndex][queryRegionIndex]);
                    }
                }
                continue;
            }
            if (targetRegions[targetIndex].empty() || queryRegions[queryIndex].empty()) {
                groups.push_back({targetIndex, queryIndex, 0U, 0U, {}});
                continue;
            }
            for (std::size_t targetRegionIndex = 0U;
                 targetRegionIndex < targetRegions[targetIndex].size(); ++targetRegionIndex) {
                for (std::size_t queryRegionIndex = 0U;
                     queryRegionIndex < queryRegions[queryIndex].size(); ++queryRegionIndex) {
                    const auto groupIndex = groups.size();
                    groups.push_back({targetIndex, queryIndex,
                                      targetRegionIndex, queryRegionIndex, {}});
                    appendDomains(groupIndex, targetIndex, queryIndex,
                        targetRegions[targetIndex][targetRegionIndex],
                        queryRegions[queryIndex][queryRegionIndex]);
                }
            }
        }

        std::vector<intarnanew::OutputGroupKey> outputGroups;
        outputGroups.reserve(groups.size());
        for (const auto& group : groups) {
            outputGroups.push_back({
                group.targetIndex, group.queryIndex,
                group.targetRegionIndex, group.queryRegionIndex,
            });
        }
        auto outputPlan = intarnanew::planOutputs(
            config, targets.size(), queries.size(), outputGroups);
        if (!outputPlan) {
            std::cerr << "IntaRNAnew: output plan: " << outputPlan.error() << '\n';
            return 2;
        }

        // A complete site list is only needed by site/profile auxiliary
        // outputs or to merge overlapping work domains. Ordinary primary
        // output consumes the globally ranked interactions alone. This is the
        // same output-sensitive distinction made by Legacy's needZall/needBPs
        // flags, but direct Predictor users retain the conservative defaults.
        const bool oneDomainPerGroup = std::ranges::all_of(groups, [](const OutputGroup& group) {
            return group.workIndices.size() <= 1U;
        });
        if (oneDomainPerGroup &&
            config.output.overlap == intarnanew::OverlapPolicy::both &&
            !needsCompleteSiteEnsemble(config)) {
            config.predictionRequirements.retainAllSites = false;
        }
        config.predictionRequirements.computeInteractionPartition =
            needsInteractionPartition(config);
        config.predictionRequirements.traceback =
            needsInteractionTraceback(config);
        std::vector<WorkResult> workResults(workItems.size());
        const auto workerCount = std::min<std::size_t>(
            std::max<std::size_t>(1U, config.threads), std::max<std::size_t>(1U, workItems.size()));
        std::vector<std::unique_ptr<intarnanew::Predictor>> predictors;
        predictors.reserve(workerCount);
        try {
            for (std::size_t index = 0U; index < workerCount; ++index) {
                predictors.push_back(std::make_unique<intarnanew::Predictor>(config));
            }
        } catch (const std::exception& exception) {
            std::cerr << "IntaRNAnew: failed to initialize prediction workers: "
                      << exception.what() << '\n';
            return 1;
        }

        auto execution = intarnanew::runParallelIndexed(
            workItems.size(), workerCount,
            [&](const std::size_t workerIndex,
                const std::size_t taskIndex,
                const std::stop_token) {
                const auto& task = workItems[taskIndex];
                workResults[taskIndex].prediction = predictors[workerIndex]->predict(
                    targets[task.targetIndex], queries[task.queryIndex],
                    *targetAccessibility[task.targetIndex], *queryAccessibility[task.queryIndex],
                    task.targetDomain, task.queryDomain);
            });
        if (!execution) {
            if (execution.error().taskIndex == std::numeric_limits<std::size_t>::max()) {
                std::cerr << "IntaRNAnew: " << execution.error().message << '\n';
            } else {
                const auto& task = workItems[execution.error().taskIndex];
                std::cerr << "IntaRNAnew: prediction failed for target " << task.targetIndex + 1U
                          << ", query " << task.queryIndex + 1U << ": "
                          << execution.error().message << '\n';
            }
            return 1;
        }

        std::vector<intarnanew::PredictionResult> results;
        results.reserve(groups.size());
        for (const auto& group : groups) {
            std::vector<intarnanew::PredictionResult> predictions;
            predictions.reserve(group.workIndices.size());
            for (const auto workIndex : group.workIndices) {
                predictions.push_back(std::move(workResults[workIndex].prediction));
            }
            intarnanew::PredictionResult prediction;
            if (predictions.empty()) {
                prediction.rt = predictors.front()->rt();
            } else if (predictions.size() == 1U) {
                prediction = std::move(predictions.front());
            } else {
                prediction = intarnanew::reducePredictions(predictions, config);
            }
            if (const auto value = targetAccessibility[group.targetIndex]->ensembleLogPartition()) {
                prediction.targetLogPartition = *value;
            }
            if (const auto value = queryAccessibility[group.queryIndex]->ensembleLogPartition()) {
                prediction.queryLogPartition = *value;
            }
            if (const auto value = targetAccessibility[group.targetIndex]->ensembleFreeEnergy()) {
                prediction.targetEnsembleFreeEnergy = *value;
            }
            if (const auto value = queryAccessibility[group.queryIndex]->ensembleFreeEnergy()) {
                prediction.queryEnsembleFreeEnergy = *value;
            }

            results.push_back(std::move(prediction));
        }

        std::vector<intarnanew::OutputArtifact> artifacts;
        artifacts.reserve(outputPlan->publications.size());
        for (const auto& publication : outputPlan->publications) {
            std::string content;
            bool includeCsvHeader = true;
            for (const auto& part : publication.parts) {
                if (part.kind == intarnanew::OutputPartKind::primary) {
                    const auto& group = groups[part.groupIndex];
                    auto formatted = intarnanew::OutputFormatter::primary(
                        config, targets[group.targetIndex], queries[group.queryIndex],
                        results[part.groupIndex], includeCsvHeader);
                    if (!formatted) {
                        std::cerr << "IntaRNAnew: output failed for target "
                                  << group.targetIndex + 1U << ", query "
                                  << group.queryIndex + 1U << ": " << formatted.error() << '\n';
                        return 1;
                    }
                    content += *formatted;
                    includeCsvHeader = false;
                    continue;
                }

                const auto& descriptor = config.output.destinations[part.descriptorIndex];
                std::size_t targetIndex{};
                std::size_t queryIndex{};
                const intarnanew::PredictionResult* prediction{};
                intarnanew::PredictionResult emptyPrediction;
                if (part.kind == intarnanew::OutputPartKind::pairAuxiliary) {
                    const auto& group = groups[part.groupIndex];
                    targetIndex = group.targetIndex;
                    queryIndex = group.queryIndex;
                    prediction = &results[part.groupIndex];
                } else if (part.kind == intarnanew::OutputPartKind::queryAccessibility) {
                    queryIndex = part.sequenceIndex;
                    prediction = &emptyPrediction;
                } else {
                    targetIndex = part.sequenceIndex;
                    prediction = &emptyPrediction;
                }
                auto formatted = intarnanew::OutputFormatter::auxiliary(
                    descriptor, config, targets[targetIndex], queries[queryIndex],
                    *targetAccessibility[targetIndex], *queryAccessibility[queryIndex],
                    *prediction);
                if (!formatted) {
                    std::cerr << "IntaRNAnew: auxiliary output failed for target "
                              << targetIndex + 1U << ", query " << queryIndex + 1U
                              << ": " << formatted.error() << '\n';
                    return 1;
                }
                content += *formatted;
            }
            artifacts.push_back({publication.destination, std::move(content)});
        }
        if (auto status = intarnanew::publishOutputs(artifacts); !status) {
            std::cerr << "IntaRNAnew: " << status.error() << '\n';
            return 1;
        }
        return 0;
    } catch (const std::exception& exception) {
        std::cerr << "IntaRNAnew: unexpected failure: " << exception.what() << '\n';
        return 1;
    }
}
