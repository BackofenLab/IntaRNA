#include "intarnanew/cli.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <set>
#include <string>
#include <string_view>
#include <vector>

namespace {

int failures{};

void check(const bool condition, const std::string_view description) {
    if (!condition) {
        std::cerr << "FAILED: " << description << '\n';
        ++failures;
    }
}

[[nodiscard]] auto lower(std::string_view value) -> std::string {
    std::string result(value);
    std::ranges::transform(result, result.begin(), [](const unsigned char character) {
        return static_cast<char>(std::tolower(character));
    });
    return result;
}

[[nodiscard]] auto fixtureRoot() -> std::filesystem::path {
    if (const auto* configured = std::getenv("INTARNANEW_SOURCE_DIR")) {
        return std::filesystem::path(configured) / "tests" / "fixtures";
    }
    return std::filesystem::path("tests") / "fixtures";
}

[[nodiscard]] auto parse(
    const std::vector<std::string_view>& arguments,
    const std::string_view invocation = {}) {
    return intarnanew::Cli::parse(arguments, invocation);
}

} // namespace

auto main() -> int {
    using namespace intarnanew;

    const auto registry = Cli::optionRegistry();
    check(registry.size() == 92U, "registry exposes the complete 92-option public surface");
    std::set<std::string> longNames;
    std::set<char> shortNames;
    for (const auto& option : registry) {
        check(!option.longName.empty(), "every option has a canonical long name");
        check(longNames.insert(lower(option.longName)).second,
              "long option names are ASCII-case-insensitively unique");
        if (option.shortName != '\0') {
            const auto normalized = static_cast<char>(
                std::tolower(static_cast<unsigned char>(option.shortName)));
            check(shortNames.insert(normalized).second,
                  "short aliases are ASCII-case-insensitively unique");
        }
        check(!option.description.empty(), "every registry entry documents its behavior");
        if (option.valueMode == OptionValueMode::required) {
            check(!option.valueName.empty(), "required values have a help metavariable");
        }
    }
    check(shortNames == std::set<char>{'e', 'h', 'm', 'n', 'q', 't', 'v'},
          "documented short alias set is complete");

    const auto fullHelp = Cli::help(true);
    for (const auto& option : registry) {
        check(fullHelp.find("--" + std::string(option.longName)) != std::string::npos,
              "full help contains every registry option");
    }
    check(fullHelp.find("[compatibility-only]") != std::string::npos,
          "full help identifies accepted compatibility-only controls");
    check(fullHelp.find("qMinE/qSpotProb/qAcc/qPu") != std::string::npos,
          "full help enumerates supported auxiliary output prefixes");
    check(fullHelp.find("IntaRNA/3/1/2/exact/helix/duplex/sTar/seed/ens") !=
              std::string::npos,
          "full help enumerates supported personalities");
    const auto basicHelp = Cli::help(false);
    for (const auto& option : registry) {
        if (option.basic) {
            check(basicHelp.find("--" + std::string(option.longName)) != std::string::npos,
                  "basic help is generated from basic registry entries");
        }
    }

    {
        const auto parameter = (fixtureRoot() / "config" / "base.parameter").string();
        for (const auto& option : registry) {
            std::vector<std::string> owned;
            if (option.longName != "help" && option.longName != "fullhelp" &&
                option.longName != "version") {
                owned.emplace_back("--help");
            }
            std::string spelling = "--" + std::string(option.longName);
            if (option.valueMode == OptionValueMode::required) {
                std::string value(option.defaultValue);
                if (option.longName == "query") value = "G";
                else if (option.longName == "target") value = "C";
                else if (option.longName == "parameterFile") value = parameter;
                else if (value.empty()) value = "x";
                spelling += "=" + value;
            }
            owned.push_back(std::move(spelling));
            std::vector<std::string_view> arguments;
            arguments.reserve(owned.size());
            for (const auto& argument : owned) arguments.push_back(argument);
            check(parse(arguments).has_value(),
                  "every registry entry is routed by the command-line parser");
        }
    }

    {
        const std::vector<std::string_view> arguments{
            "--TaRgEt=C", "--QuErY=G", "--EnErGy=b", "--AcC=n", "--NoSeEd",
            "-M", "h", "-N", "2", "--OuTmOdE=c",
        };
        const auto config = parse(arguments);
        check(config.has_value(), "long/short names and documented choices are case-insensitive");
        if (config) {
            check(config->energy == EnergyKind::basePair, "lower-case energy choice is applied");
            check(config->output.mode == OutputMode::csv, "lower-case output choice is applied");
            check(config->output.number == 2U, "upper-case short alias is applied");
        }
    }

    {
        const auto parameter = (fixtureRoot() / "config" / "base.parameter").string();
        const std::vector<std::string_view> arguments{
            "--personality=IntaRNAexact",
            "--parameterFile", parameter,
            "--mode=S",
            "--seedBP=2",
            "--out=STDERR",
        };
        const auto config = parse(arguments, "/tmp/IntaRNAduplex");
        check(config.has_value(), "personality, nested files, and CLI form one valid layered config");
        if (config) {
            check(config->personality == "IntaRNAexact",
                  "CLI personality overrides the executable personality alias");
            check(config->query.accessibility == AccessibilityKind::compute &&
                  config->target.accessibility == AccessibilityKind::compute,
                  "chosen personality replaces lower-precedence executable preset");
            check(config->query.accessibilityWindow == 0U &&
                  config->target.accessibilityWindow == 0U,
                  "personality defaults apply before parameter and CLI assignments");
            check(config->mode == PredictionMode::seedOnly,
                  "CLI scalar overrides nested and outer parameter-file values");
            check(config->query.regions == "2-3",
                  "nested parameter-file values are loaded relative to the including file");
            check(config->parameterFiles.size() == 2U,
                  "top-level and nested parameter files are retained");
            check(config->output.destinations.size() == 3U &&
                  config->output.destinations[0] == "file-from-base" &&
                  config->output.destinations[1] == "qMinE:tracker.csv" &&
                  config->output.destinations[2] == "STDERR",
                  "repeatable --out values retain deterministic layer/file order");
            const auto modeOrigin = config->provenance.find("mode");
            check(modeOrigin != config->provenance.end() &&
                  modeOrigin->second.source == ConfigSource::commandLine,
                  "effective scalar provenance identifies the winning CLI source");
            check(config->assignmentHistory.size() >= 12U,
                  "assignment history retains overridden sources");
        }
    }

    {
        const std::vector<std::pair<std::string_view, InteractionModel>> personalities{
            {"IntaRNA", InteractionModel::seedExtension},
            {"IntaRNA3", InteractionModel::seedExtension},
            {"IntaRNA1", InteractionModel::singleSite},
            {"IntaRNA2", InteractionModel::singleSite},
            {"IntaRNAexact", InteractionModel::seedExtension},
            {"IntaRNAhelix", InteractionModel::helixBlocks},
            {"IntaRNAduplex", InteractionModel::seedExtension},
            {"IntaRNAsTar", InteractionModel::seedExtension},
            {"IntaRNAseed", InteractionModel::seedExtension},
            {"IntaRNAens", InteractionModel::ensemble},
        };
        for (const auto& [personality, expectedModel] : personalities) {
            const std::vector<std::string_view> arguments{
                "--query=GGGGGGG", "--target=CCCCCCC", "--energy=B",
                "--personality", personality,
            };
            const auto config = parse(arguments);
            check(config.has_value(), "documented personality is accepted");
            if (config) check(config->model == expectedModel, "personality model preset is correct");
        }
        const std::vector<std::string_view> unknown{"--personality=IntaRNAup", "--help"};
        const auto rejected = parse(unknown);
        check(!rejected && rejected.error().find("unknown personality") != std::string::npos,
              "stale/unknown personality names are rejected precisely");
    }

    {
        const std::vector<std::string_view> exactOverride{
            "--query=GGGG", "--target=CCCC", "--personality=IntaRNAexact",
            "--accW=42", "--mode=H",
        };
        const auto config = parse(exactOverride);
        check(config && config->query.accessibilityWindow == 42U &&
              config->target.accessibilityWindow == 42U &&
              config->mode == PredictionMode::heuristic,
              "CLI explicitly overrides personality defaults regardless of argument order");
    }

    {
        const std::vector<std::string_view> aliases{
            "--query=GGGG", "--target=CCCC", "--energy=B",
        };
        const auto duplex = parse(aliases, "/opt/bin/IntaRNAduplex");
        check(duplex && duplex->personality == "IntaRNAduplex" &&
              duplex->query.accessibility == AccessibilityKind::disabled,
              "executable basename selects a personality");
        const auto generic = parse(aliases, "/opt/bin/custom-runner");
        check(generic && generic->personality == "IntaRNA",
              "an embedded non-IntaRNA executable name uses the standard personality");
        check(generic && generic->query.accessibility == AccessibilityKind::compute,
              "generic invocation retains standard computed accessibility");
    }

    {
        const auto conflictFile = (fixtureRoot() / "config" / "conflict.parameter").string();
        const std::vector<std::string_view> arguments{
            "--parameterFile", conflictFile, "--qRegionLenMax=10",
        };
        const auto rejected = parse(arguments);
        check(!rejected &&
              rejected.error().find("parameter file") != std::string::npos &&
              rejected.error().find("command line argument") != std::string::npos,
              "cross-layer conflicts cite both effective origins");
    }

    {
        const auto duplicateFile = (fixtureRoot() / "config" / "duplicate.parameter").string();
        const std::vector<std::string_view> arguments{"--parameterFile", duplicateFile};
        const auto rejected = parse(arguments);
        check(!rejected &&
              rejected.error().find("duplicate --query") != std::string::npos &&
              rejected.error().find("line 1") != std::string::npos &&
              rejected.error().find("line 2") != std::string::npos,
              "same-source scalar duplicates cite both parameter lines");
    }

    {
        const std::vector<std::string_view> supported{
            "--query=GGGG", "--target=CCCC", "--energy=B", "--acc=C",
            "--qPfScale=1.2", "--tPfScale=2", "--accNoLP", "--accNoGUend",
            "--outBestSeedOnly", "--verbose", "--default-log-file=diagnostics.log",
        };
        const auto config = parse(supported);
        check(config.has_value(), "new accessibility and compatibility controls are accepted");
        if (config) {
            check(config->query.partitionScale == 1.2 &&
                  config->target.partitionScale == 2.0,
                  "per-side partition scaling is stored");
            check(config->accessibilityNoLonelyPairs &&
                  config->accessibilityNoGuAtEnds,
                  "accessibility folding flags are stored");
            check(config->output.bestSeedOnly && config->verbose &&
                  config->logFile == "diagnostics.log",
                  "compatibility-only controls are recorded deterministically");
        }
    }

    {
        const std::vector<std::string_view> invalidSeedOnly{
            "--query=G", "--target=C", "--energy=B", "--acc=N",
            "--mode=S", "--noSeed",
        };
        const auto rejected = parse(invalidSeedOnly);
        check(!rejected && rejected.error().find("seed-only") != std::string::npos,
              "seed-only mode rejects a disabled seed without an explicit seed");

        const std::vector<std::string_view> invalidHelix{
            "--query=GG", "--target=CC", "--energy=B", "--acc=N",
            "--model=B", "--mode=M",
        };
        const auto helixRejected = parse(invalidHelix);
        check(!helixRejected && helixRejected.error().find("require heuristic") != std::string::npos,
              "helix-block model rejects a non-heuristic mode");

        const std::vector<std::vector<std::string_view>> invalidEnsembleRequests{
            {"--query=GGGGGGG", "--target=CCCCCCC", "--energy=B", "--model=X", "--outMode=E"},
            {"--query=GGGGGGG", "--target=CCCCCCC", "--energy=B", "--model=X",
             "--outMode=C", "--outCsvCols=id1,P_E"},
            {"--query=GGGGGGG", "--target=CCCCCCC", "--energy=B", "--model=X",
             "--out=spotProb:STDOUT"},
        };
        for (const auto& request : invalidEnsembleRequests) {
            const auto ensembleRejected = parse(request);
            check(!ensembleRejected &&
                  ensembleRejected.error().find("no scientifically valid partition") !=
                      std::string::npos,
                  "seeded model X rejects ensemble-dependent output");
        }
        const std::vector<std::string_view> unseededEnsemble{
            "--query=GGGG", "--target=CCCC", "--energy=B", "--model=X",
            "--noSeed", "--outMode=E",
        };
        check(parse(unseededEnsemble).has_value(),
              "unseeded model X retains scientifically composable ensemble output");
    }

    {
        const auto parameterDirectory = fixtureRoot() / "parameters";
        std::size_t count{};
        for (const auto& entry : std::filesystem::directory_iterator(parameterDirectory)) {
            if (!entry.is_regular_file() || entry.path().extension() != ".parameter") continue;
            const auto path = entry.path().string();
            const std::vector<std::string_view> arguments{"--parameterFile", path};
            const auto config = parse(arguments);
            check(config.has_value(), "normalized public parameter fixture parses");
            ++count;
        }
        check(count == 17U, "all 17 normalized public parameter fixtures are present");
    }

    if (failures != 0) {
        std::cerr << failures << " registry test(s) failed\n";
        return 1;
    }
    std::cout << "all CLI/config registry tests passed\n";
    return 0;
}
