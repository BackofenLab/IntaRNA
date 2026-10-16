#include "intarnanew/accessibility.hpp"
#include "intarnanew/cli.hpp"
#include "intarnanew/output.hpp"
#include "intarnanew/sequence.hpp"

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <string_view>
#include <vector>

namespace {

int failureCount{};

void check(const bool condition, const std::string_view description) {
    if (!condition) {
        std::cerr << "FAILED: " << description << '\n';
        ++failureCount;
    }
}

[[nodiscard]] auto occurrences(
    const std::string_view text,
    const std::string_view needle) -> std::size_t {
    std::size_t count{};
    std::size_t position{};
    while ((position = text.find(needle, position)) != std::string_view::npos) {
        ++count;
        position += needle.size();
    }
    return count;
}

[[nodiscard]] auto outputInteraction(
    const std::size_t targetIndex,
    const std::string& targetId = "target") -> intarnanew::Interaction {
    intarnanew::Interaction interaction;
    interaction.targetId = targetId;
    interaction.queryId = "query";
    interaction.pairs.push_back({targetIndex, 0U});
    interaction.energy.initiation = -1.0;
    return interaction;
}

} // namespace

auto main() -> int {
    using namespace intarnanew;

    {
        const std::vector<std::string_view> arguments{
            "--TaRgEt=C", "--QuErY=G", "--EnErGy=B", "--AcC=N", "--NoSeEd",
            "--OuTmOdE=C", "--QiDxPoS0", "-2",
        };
        auto parsed = Cli::parse(arguments);
        check(parsed.has_value(), "long option names are case-insensitive");
        if (parsed) {
            check(parsed->query.firstPosition == -2, "space-separated negative numeric value is consumed");
            check(parsed->output.mode == OutputMode::csv, "mixed-case output option was applied");
            check(!parsed->seed.required, "mixed-case Boolean option was applied");
        }
    }

    {
        const auto rejected = [](const std::string_view option) {
            const std::vector<std::string_view> arguments{
                "--target=C", "--query=G", "--energy=B", "--acc=N", "--noSeed", option,
            };
            return Cli::parse(arguments);
        };
        check(!rejected("--windowWidth=4"), "prediction windows below the documented minimum are rejected");
        check(!rejected("--qShape=shape.txt"), "incomplete SHAPE configuration is rejected");
        check(rejected("--tRegionLenMax=20").has_value(),
              "implemented automatic-region configuration is accepted");
        check(rejected("--outPerRegion").has_value(),
              "implemented per-region output configuration is accepted");
        check(rejected("--qPfScale=1.2").has_value(),
              "implemented partition scaling is accepted");
        check(rejected("--accNoLP").has_value(),
              "implemented accessibility noLP control is accepted");
        check(rejected("--outBestSeedOnly").has_value(),
              "implemented best-seed filtering is accepted");
        check(rejected("--verbose").has_value(),
              "compatibility logging control is accepted and recorded");
    }

    {
        const Sequence sequence("coordinates", "AAAA", -2);
        check(sequence.externalIndex(0U) == -2, "negative origin first coordinate");
        check(sequence.externalIndex(1U) == -1, "coordinate before skipped zero");
        check(sequence.externalIndex(2U) == 1, "coordinate zero is skipped");
        check(sequence.externalIndex(3U) == 2, "coordinate after skipped zero");
        const auto internal = sequence.internalIndex(1);
        check(internal && *internal == 2U, "positive external coordinate maps across skipped zero");
        check(!sequence.internalIndex(0), "external coordinate zero is rejected for negative origins");
    }

    {
        const Sequence sequence("iupac", "ACGURYSWKMBDHVN");
        check(sequence.str() == "ACGUNNNNNNNNNNN", "non-canonical IUPAC input normalizes to N");
        check(!canPair('N', 'A'), "N is non-pairing");
        check(!canPair('R', 'Y'), "ambiguous public symbols are non-pairing");
        check(canPair('C', 'G'), "canonical complementary bases still pair");
    }

    {
        std::vector<Sequence> sequences;
        sequences.emplace_back("one", "A");
        sequences.emplace_back("two", "C");
        sequences.emplace_back("three", "G");
        auto selected = SequenceReader::select(std::move(sequences), "3, 1");
        check(selected && selected->size() == 2U, "sequence subset selects requested records");
        if (selected && selected->size() == 2U) {
            check((*selected)[0].id() == "one" && (*selected)[1].id() == "three",
                  "sequence subset preserves FASTA input order");
        }
        std::vector<Sequence> invalidSequences;
        invalidSequences.emplace_back("one", "A");
        check(!SequenceReader::select(std::move(invalidSequences), "0"),
              "sequence subset rejects zero index");
    }

    const Sequence target("target", "CCCCCCCCCCCC");
    const Sequence query("query", "G");
    PredictionResult result;
    result.rt = 1.0;
    result.interactions.push_back(outputInteraction(9U));
    result.interactions.push_back(outputInteraction(1U));
    result.interactions.push_back(outputInteraction(2U));

    {
        Config config;
        config.output.mode = OutputMode::csv;
        config.output.csvColumns = "start1,E";
        config.output.csvSort = "start1";
        auto formatted = OutputFormatter::primary(config, target, query, result);
        check(formatted && *formatted == "start1;E\n2;-1\n3;-1\n10;-1\n",
              "CSV coordinate sorting is numeric rather than lexical");

        auto withoutHeader = OutputFormatter::primary(config, target, query, result, false);
        check(withoutHeader && !withoutHeader->starts_with("start1;E"),
              "subsequent CSV documents can suppress their header");
        if (formatted && withoutHeader) {
            check(occurrences(*formatted + *withoutHeader, "start1;E\n") == 1U,
                  "combined multi-pair CSV has exactly one header");
        }
    }

    {
        Config config;
        config.output.mode = OutputMode::csv;
        config.output.csvColumns = "start1,P_E";
        config.output.csvSort = "P_E";
        PredictionResult probabilityResult;
        probabilityResult.rt = 1.0;
        probabilityResult.logPartition = 0.0;
        auto lowerProbability = outputInteraction(0U);
        lowerProbability.energy.initiation = -10.0;
        lowerProbability.probability = 0.2;
        auto higherProbability = outputInteraction(1U);
        higherProbability.probability = 0.8;
        probabilityResult.interactions.push_back(std::move(lowerProbability));
        probabilityResult.interactions.push_back(std::move(higherProbability));
        auto formatted = OutputFormatter::primary(
            config, target, query, probabilityResult);
        check(formatted && *formatted == "start1;P_E\n1;0.2\n2;0.8\n",
              "P_E formatting and sorting use the normalized site probability");
    }

    {
        Config config;
        config.output.mode = OutputMode::csv;
        config.output.csvColumns = "id1,E";
        PredictionResult escapedResult;
        escapedResult.rt = 1.0;
        escapedResult.interactions.push_back(outputInteraction(0U, "target;quoted"));
        auto formatted = OutputFormatter::primary(config, target, query, escapedResult);
        check(formatted && formatted->find("\"target;quoted\";-1") != std::string::npos,
              "CSV values containing the separator are quoted");
    }

    {
        Config config;
        config.output.mode = OutputMode::ensemble;
        PredictionResult empty;
        empty.rt = 0.6;
        empty.targetLogPartition = 0.0;
        empty.queryLogPartition = 0.0;
        empty.targetEnsembleFreeEnergy = -1.31;
        empty.queryEnsembleFreeEnergy = -0.25;
        auto formatted = OutputFormatter::primary(config, target, query, empty);
        check(formatted && *formatted ==
            "id1 target\nid2 query\nRT 0.6\nEall 0.00\nEall1 -1.31\n"
            "Eall2 -0.25\nEallTotal 0.00\n",
              "empty ensemble keeps monomer summaries but reports no formation total");
        check(formatted && formatted->find("inf") == std::string::npos,
              "non-finite implementation spelling is not emitted");

        PredictionResult populated;
        populated.rt = 1.0;
        populated.ensembleSites.push_back(outputInteraction(0U));
        populated.ensembleFreeEnergy = -5.889;
        populated.targetEnsembleFreeEnergy = -1.239;
        populated.queryEnsembleFreeEnergy = -2.349;
        auto populatedText = OutputFormatter::primary(config, target, query, populated);
        check(populatedText && *populatedText ==
            "id1 target\nid2 query\nRT 1\nEall -5.88\nEall1 -1.23\n"
            "Eall2 -2.34\nEallTotal -9.47\n",
              "ensemble energies use centikcal truncation and include monomer totals");
        PredictionResult signedZero;
        signedZero.rt = 1.0;
        signedZero.ensembleSites.push_back(outputInteraction(0U));
        signedZero.ensembleFreeEnergy = -1.0;
        signedZero.targetEnsembleFreeEnergy = -0.0;
        signedZero.queryEnsembleFreeEnergy = -0.0;
        const auto signedZeroText =
            OutputFormatter::primary(config, target, query, signedZero);
        check(signedZeroText &&
              signedZeroText->find("Eall1 0.00\nEall2 0.00\n") != std::string::npos,
              "fixed ensemble output canonicalizes negative zero");
    }

    {
        Config config;
        DisabledAccessibility accessibility(target);
        auto formatted = OutputFormatter::auxiliary(
            "QACC:ignored", config, target, query, accessibility, accessibility, result);
        check(formatted.has_value(), "auxiliary output descriptor names are case-insensitive");
    }

    {
        const Sequence encodedTarget("target", "ACG", 5);
        const Sequence encodedQuery("query", "UGC", -2);
        Interaction interaction;
        interaction.targetId = "target";
        interaction.queryId = "query";
        interaction.pairs = {{0U, 2U}, {2U, 0U}};
        interaction.energy.initiation = -2.0;
        interaction.energy.additive = 0.5;
        interaction.seeds = {SeedMatch{0U, 1U, -4.0, 0.5, 0.25, 0.6, 0.7}};
        PredictionResult encodedResult;
        encodedResult.rt = 1.0;
        encodedResult.interactions.push_back(interaction);
        encodedResult.logPartition = std::log(10.0);
        encodedResult.ensembleFreeEnergy = -3.0;
        encodedResult.targetLogPartition = std::log(2.0);
        encodedResult.queryLogPartition = std::log(3.0);
        encodedResult.targetEnsembleFreeEnergy = -1.0;
        encodedResult.queryEnsembleFreeEnergy = -2.0;

        Config config;
        config.output.mode = OutputMode::csv;
        config.output.csvColumns =
            "subseqDB,hybridDB,bpList,E,Etotal,E_hybrid,E_add,"
            "seedE,seedED1,seedED2,seedPu1,seedPu2,"
            "Eall1,Eall2,Zall1,Zall2,EallTotal";
        auto formatted = OutputFormatter::primary(
            config, encodedTarget, encodedQuery, encodedResult);
        check(formatted && *formatted ==
            "subseqDB;hybridDB;bpList;E;Etotal;E_hybrid;E_add;"
            "seedE;seedED1;seedED2;seedPu1;seedPu2;"
            "Eall1;Eall2;Zall1;Zall2;EallTotal\n"
            "5ACG&-2UGC;5|.|&-2|.|;(5,1):(7,-2);-1.5;-4.5;-1.5;0.5;"
            "-4;0.5;0.25;0.6;0.7;"
            "-1;-2;2;3;-6\n",
            "CSV sequence, bar, pair-list, and monomer-ensemble fields follow the public contract");

        encodedResult.interactions.front().seeds.clear();
        config.output.csvColumns =
            "seedStart1,seedEnd1,seedStart2,seedEnd2,seedE,"
            "seedED1,seedED2,seedPu1,seedPu2";
        auto seedless = OutputFormatter::primary(
            config, encodedTarget, encodedQuery, encodedResult);
        check(seedless && *seedless ==
            "seedStart1;seedEnd1;seedStart2;seedEnd2;seedE;"
            "seedED1;seedED2;seedPu1;seedPu2\n"
            "NAN;NAN;NAN;NAN;NAN;NAN;NAN;NAN;NAN\n",
            "seedless CSV fields use the public NAN convention");
    }

    {
        const Sequence chartTarget("target", "UUUUUUCCCUAAAUCCUUCUUU");
        const Sequence chartQuery("query", "GGGGGGGGGGGGGGGGGGGGG");
        Interaction interaction;
        interaction.targetId = "target";
        interaction.queryId = "query";
        interaction.pairs = {{7U, 6U}, {8U, 5U}, {14U, 4U}, {15U, 3U}};
        interaction.energy.initiation = -1.0;
        interaction.energy.loops = -3.0;
        PredictionResult chartResult;
        chartResult.rt = 1.0;
        chartResult.interactions.push_back(interaction);
        Config config;
        config.output.mode = OutputMode::normal;
        auto formatted = OutputFormatter::primary(
            config, chartTarget, chartQuery, chartResult);
        check(formatted && *formatted ==
            "\ntarget\n"
            "             8       16\n"
            "             |       |\n"
            "   5'-UUUUUUC  UAAAU  UUCUUU-3'\n"
            "             CC     CC\n"
            "             ||     ||\n"
            "             GG     GG\n"
            "3'-GGG...GGGG  -----  GGG-5'\n"
            "             |       |\n"
            "             7       4\n"
            "query\n\n"
            "interaction energy = -4 kcal/mol\n",
            "normal ASCII output is byte-compatible with the public no-GU-end fixture");
    }

    {
        const Sequence chartTarget("target", "UCCC", -2);
        const Sequence chartQuery("query", "GGGG", -2);
        Interaction interaction;
        interaction.targetId = "target";
        interaction.queryId = "query";
        interaction.pairs = {{0U, 3U}, {1U, 2U}, {2U, 1U}, {3U, 0U}};
        interaction.energy.initiation = -1.0;
        interaction.energy.loops = -3.0;
        interaction.energy.additive = 2.5;
        interaction.unpairedTarget = 0.35;
        interaction.unpairedQuery = 0.45;
        interaction.seeds = {SeedMatch{1U, 2U, -2.34, 0.5, 0.25, 0.6, 0.7}};
        PredictionResult chartResult;
        chartResult.rt = 1.0;
        chartResult.interactions.push_back(interaction);
        Config config;
        config.output.mode = OutputMode::detailed;
        auto formatted = OutputFormatter::primary(
            config, chartTarget, chartQuery, chartResult);
        check(formatted && *formatted ==
            "\ntarget\n"
            "            -2  +2\n"
            "             |  |\n"
            "          5'-    -3'\n"
            "             UCCC\n"
            "             :++|\n"
            "             GGGG\n"
            "          3'-    -5'\n"
            "             |  |\n"
            "            +2  -2\n"
            "query\n\n"
            "interaction seq1   = -2..2\n"
            "interaction seq2   = -2..2\n\n"
            "interaction energy = -1.5 kcal/mol\n"
            "  = E(init)        = -1\n"
            "  + E(loops)       = -3\n"
            "  + E(dangleLeft)  = 0\n"
            "  + E(dangleRight) = 0\n"
            "  + E(endLeft)     = 0\n"
            "  + E(endRight)    = 0\n"
            "    : E(hybrid)    = -1.5\n"
            "  + ED(seq1)       = 0\n"
            "    : Pu(seq1)     = 0.35\n"
            "  + ED(seq2)       = 0\n"
            "    : Pu(seq2)     = 0.45\n"
            "  + E(add)         = 2.5\n\n"
            "seed seq1   = -1..1\n"
            "seed seq2   = -1..1\n"
            "seed energy = -2.34\n"
            "seed ED1    = 0.5\n"
            "seed ED2    = 0.25\n"
            "seed Pu1    = 0.6\n"
            "seed Pu2    = 0.7\n",
            "detailed ASCII output renders stored native Pu values without RT reconstruction");
    }

    {
        const Sequence seedTarget("target", "CCCC");
        const Sequence seedQuery("query", "GGGG");
        Interaction interaction;
        interaction.targetId = "target";
        interaction.queryId = "query";
        interaction.pairs = {{0U, 3U}, {1U, 2U}, {2U, 1U}, {3U, 0U}};
        interaction.energy.loops = -4.0;
        interaction.seeds = {
            SeedMatch{2U, 3U, -2.0},
            SeedMatch{0U, 1U, -2.0},
            SeedMatch{1U, 2U, -2.0},
        };
        PredictionResult seedResult;
        seedResult.rt = 1.0;
        seedResult.interactions.push_back(interaction);

        Config csvConfig;
        csvConfig.output.mode = OutputMode::csv;
        csvConfig.output.csvColumns =
            "seedStart1,seedEnd1,seedStart2,seedEnd2,seedE";
        auto csvAll = OutputFormatter::primary(
            csvConfig, seedTarget, seedQuery, seedResult);
        check(csvAll && *csvAll ==
            "seedStart1;seedEnd1;seedStart2;seedEnd2;seedE\n"
            "1:2:3;2:3:4;3:2:1;4:3:2;-2:-2:-2\n",
            "CSV emits all seeds colon-delimited in energy/coordinate order");

        csvConfig.output.bestSeedOnly = true;
        auto csvBest = OutputFormatter::primary(
            csvConfig, seedTarget, seedQuery, seedResult);
        check(csvBest && *csvBest ==
            "seedStart1;seedEnd1;seedStart2;seedEnd2;seedE\n"
            "1;2;3;4;-2\n",
            "outBestSeedOnly restricts CSV seed fields to the deterministic best seed");

        Config detailedConfig;
        detailedConfig.output.mode = OutputMode::detailed;
        auto detailed = OutputFormatter::primary(
            detailedConfig, seedTarget, seedQuery, seedResult);
        check(detailed &&
              detailed->find("             ++++\n") != std::string::npos &&
              detailed->find("seed seq1   = 1..2 | 2..3 | 3..4\n") != std::string::npos &&
              detailed->find("seed seq2   = 3..4 | 2..3 | 1..2\n") != std::string::npos &&
              detailed->find("seed energy = -2 | -2 | -2\n") != std::string::npos,
            "detailed output joins all seeds with pipes and marks their chart union");

        detailedConfig.output.bestSeedOnly = true;
        auto detailedBest = OutputFormatter::primary(
            detailedConfig, seedTarget, seedQuery, seedResult);
        check(detailedBest &&
              detailedBest->find("             ++||\n") != std::string::npos &&
              detailedBest->find("seed seq1   = 1..2\n") != std::string::npos &&
              detailedBest->find("seed seq1   = 1..2 |") == std::string::npos,
            "outBestSeedOnly restricts detailed fields and chart markers to one seed");

        interaction.seeds.clear();
        seedResult.interactions = {interaction};
        detailedConfig.output.bestSeedOnly = false;
        auto noSeed = OutputFormatter::primary(
            detailedConfig, seedTarget, seedQuery, seedResult);
        check(noSeed && noSeed->find("seed seq1") == std::string::npos &&
              noSeed->find("++++") == std::string::npos,
            "interactions without seeds emit neither seed details nor seed chart markers");
    }

    {
        const Sequence trackedTarget("target", "CU", 5);
        const Sequence trackedQuery("query", "GA", -2);
        DisabledAccessibility targetAccessibility(trackedTarget);
        DisabledAccessibility queryAccessibility(trackedQuery);
        Interaction broad;
        broad.targetId = "target";
        broad.queryId = "query";
        broad.pairs = {{0U, 1U}, {1U, 0U}};
        broad.energy.initiation = -1.0;
        broad.probability = 0.6;
        Interaction focused;
        focused.targetId = "target";
        focused.queryId = "query";
        focused.pairs = {{0U, 0U}};
        focused.energy.initiation = -2.0;
        focused.probability = 0.2;
        PredictionResult tracked;
        tracked.rt = 1.0;
        tracked.ensembleSites = {broad, focused};
        Config config;

        auto qMin = OutputFormatter::auxiliary(
            "qMinE:ignored", config, trackedTarget, trackedQuery,
            targetAccessibility, queryAccessibility, tracked);
        check(qMin && *qMin == "idx;query;minE\n-2;G;-2\n-1;A;-1\n",
              "minimum-energy profiles include external coordinate, base, and legacy header");

        auto pairMin = OutputFormatter::auxiliary(
            "pMinE:ignored", config, trackedTarget, trackedQuery,
            targetAccessibility, queryAccessibility, tracked);
        check(pairMin && *pairMin ==
            "minE;G_-2;A_-1\nC_5;-2;-1\nU_6;-1;-1\n",
            "pair minimum energies use the documented target-row/query-column matrix");

        auto selectedSpots = OutputFormatter::auxiliary(
            "spotProb:5&-2,6&-1:ignored", config, trackedTarget, trackedQuery,
            targetAccessibility, queryAccessibility, tracked);
        check(selectedSpots && *selectedSpots ==
            "spot;probability\n4&-3;0.2\n5&-2;0.8\n6&-1;0.6\n",
            "specific spot output computes an overlap-safe union and coordinate-aware none label");

        auto allSpots = OutputFormatter::auxiliary(
            "spotProb::ignored", config, trackedTarget, trackedQuery,
            targetAccessibility, queryAccessibility, tracked);
        check(allSpots && *allSpots ==
            "spotProb;G_-2;A_-1\nC_5;0.8;0.6\nU_6;0.6;0.6\n",
            "empty spot list emits the documented complete probability matrix");

        auto malformedSpot = OutputFormatter::auxiliary(
            "spotProb:5x&-2:ignored", config, trackedTarget, trackedQuery,
            targetAccessibility, queryAccessibility, tracked);
        check(!malformedSpot, "spot coordinates reject trailing non-numeric characters");
    }

    {
        const Sequence tableTarget("target", "CCC", 5);
        const Sequence tableQuery("query", "GGG", -2);
        DisabledAccessibility targetAccessibility(tableTarget);
        DisabledAccessibility queryAccessibility(tableQuery);
        PredictionResult empty;
        Config config;
        config.query.accessibility = AccessibilityKind::disabled;
        config.query.interactionLengthMax = 2U;
        auto table = OutputFormatter::auxiliary(
            "qPu:ignored", config, tableTarget, tableQuery,
            targetAccessibility, queryAccessibility, empty);
        check(table && *table ==
            "#unpaired probabilities\n #i$\tl=1\t2\t\n"
            "-2\t1.000000e+00\tNA\t\n"
            "-1\t1.000000e+00\t1.000000e+00\t\n"
            "1\t1.000000e+00\t1.000000e+00\t\n",
            "accessibility output uses the public Pu matrix layout and effective length cap");
    }

    {
        const auto unique = std::to_string(
            std::chrono::steady_clock::now().time_since_epoch().count());
        const auto directory = std::filesystem::temp_directory_path() /
                               ("intarnanew-output-test-" + unique);
        std::error_code error;
        std::filesystem::create_directory(directory, error);
        check(!error, "transaction test directory was created");
        const auto destination = directory / "result.csv";
        {
            std::ofstream initial(destination);
            initial << "old";
        }
        auto status = writeOutput(destination.string(), "new\n");
        check(status.has_value(), "transactional file output commits successfully");
        std::ifstream input(destination);
        const std::string content(
            std::istreambuf_iterator<char>{input}, std::istreambuf_iterator<char>{});
        check(content == "new\n", "transactional output replaces the destination content");

        const auto directoryDestination = directory / "not-a-file";
        std::filesystem::create_directory(directoryDestination, error);
        auto rejected = writeOutput(directoryDestination.string(), "content");
        check(!rejected, "output refuses to replace a directory");
        check(std::filesystem::is_directory(directoryDestination),
              "failed output commit leaves the directory destination intact");

        std::filesystem::remove_all(directory, error);
    }

    if (failureCount == 0) {
        std::cout << "All IntaRNAnew CLI/IO/output regression tests passed.\n";
    }
    return failureCount == 0 ? 0 : 1;
}
