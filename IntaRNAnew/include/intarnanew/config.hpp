#pragma once

#include "intarnanew/types.hpp"

#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace intarnanew {

enum class RunAction { predict, help, fullHelp, version };
enum class PredictionMode { heuristic, exact, seedOnly };
enum class InteractionModel { singleSite, seedExtension, helixBlocks, ensemble };
enum class EnergyKind { basePair, nearestNeighbor };
enum class AccessibilityKind { disabled, compute, probabilitiesFile, energiesFile };
enum class OutputMode { normal, detailed, csv, ensemble };
enum class OverlapPolicy { neither, target, query, both };
enum class ConfigSource { baseline, personality, parameterFile, commandLine };

struct ConfigOrigin {
    ConfigSource source{ConfigSource::baseline};
    std::string detail{"built-in default"};
    std::size_t position{};
};

struct ConfigAssignment {
    std::string option;
    std::string value;
    ConfigOrigin origin;
};

struct SideConfig {
    std::string input;
    std::string id;
    long long firstPosition{1};
    std::string subset;
    AccessibilityKind accessibility{AccessibilityKind::compute};
    std::size_t accessibilityWindow{150};
    std::size_t accessibilitySpan{100};
    std::string accessibilityConstraint;
    std::string accessibilityFile;
    std::size_t interactionLengthMax{};
    std::size_t interactionLoopMax{};
    std::string regions;
    std::size_t regionLengthMax{};
    double partitionScale{1.07};
    std::string shapeFile;
    std::string shapeMethod;
    std::string shapeConversion;
};

struct SeedConfig {
    bool required{true};
    std::string explicitSeeds;
    std::size_t basePairs{7};
    std::size_t maxUnpaired{};
    int queryMaxUnpaired{-1};
    int targetMaxUnpaired{-1};
    Energy maxEnergy{};
    Energy maxHybridEnergy{999.0};
    double minUnpairedProbability{};
    std::string queryRanges;
    std::string targetRanges;
    bool noGu{};
    bool noGuAtEnds{};
};

struct HelixConfig {
    std::size_t minBasePairs{2};
    std::size_t maxBasePairs{10};
    std::size_t maxInternalLoop{};
    double minUnpairedProbability{};
    Energy maxEnergy{};
    bool useFullEnergy{};
};

struct OutputConfig {
    std::vector<std::string> destinations{"STDOUT"};
    char separator{';'};
    OutputMode mode{OutputMode::normal};
    std::size_t number{1};
    OverlapPolicy overlap{OverlapPolicy::both};
    Energy maxEnergy{};
    double minUnpairedProbability{};
    Energy deltaEnergy{100.0};
    bool bestSeedOnly{};
    bool noLonelyPairs{};
    bool noGuAtEnds{};
    std::string csvColumns{"id1,start1,end1,id2,start2,end2,subseqDP,hybridDP,E"};
    std::string csvSort;
    bool perRegion{};
    bool pairwise{};
};

// Internal execution requirements are conservative for library callers: a
// directly constructed Config retains the complete site ensemble and computes
// interaction and monomer partition functions as well as tracebacks. The
// executable may relax requirements after it has inspected the complete output
// plan. Keeping these separate from the public output options lets the
// predictor avoid work that no requested output can observe without weakening
// PredictionResult's default contract.
struct PredictionRequirements {
    bool retainAllSites{true};
    bool computeInteractionPartition{true};
    bool computeMonomerPartition{true};
    bool traceback{true};
};

struct Config {
    RunAction action{RunAction::predict};
    SideConfig query = [] {
        SideConfig side{};
        side.id = "query";
        side.interactionLoopMax = 16U;
        return side;
    }();
    SideConfig target = [] {
        SideConfig side{};
        side.id = "target";
        side.interactionLoopMax = 10U;
        return side;
    }();
    SeedConfig seed;
    HelixConfig helix;
    PredictionMode mode{PredictionMode::heuristic};
    InteractionModel model{InteractionModel::seedExtension};
    EnergyKind energy{EnergyKind::nearestNeighbor};
    std::string energyParameters{"Turner04"};
    bool noDangles{};
    bool accessibilityNoLonelyPairs{};
    bool accessibilityNoGuAtEnds{};
    Energy additiveEnergy{};
    double temperatureCelsius{37.0};
    std::size_t windowWidth{};
    std::size_t windowOverlap{150};
    OutputConfig output;
    unsigned int threads{1};
    std::string personality{"IntaRNA"};
    std::string parameterFile;
    std::vector<std::string> parameterFiles;
    PredictionRequirements predictionRequirements;
    bool verbose{};
    std::string logFile;
    std::map<std::string, ConfigOrigin, std::less<>> provenance;
    std::vector<ConfigAssignment> assignmentHistory;
};

} // namespace intarnanew
