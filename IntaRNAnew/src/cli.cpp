#include "intarnanew/cli.hpp"
#include "intarnanew/folding.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace intarnanew {
namespace {

enum class OptionId {
    query,
    qId,
    qIdxPos0,
    qSet,
    qAcc,
    qAccW,
    qAccL,
    qAccConstr,
    qAccFile,
    qIntLenMax,
    qIntLoopMax,
    qRegion,
    qRegionLenMax,
    qPfScale,
    target,
    tId,
    tIdxPos0,
    tSet,
    tAcc,
    tAccW,
    tAccL,
    tAccConstr,
    tAccFile,
    tIntLenMax,
    tIntLoopMax,
    tRegion,
    tRegionLenMax,
    tPfScale,
    noSeed,
    seedTQ,
    seedBP,
    seedMaxUP,
    seedQMaxUP,
    seedTMaxUP,
    seedMaxE,
    seedMaxEhybrid,
    seedMinPu,
    seedQRange,
    seedTRange,
    seedNoGU,
    seedNoGUend,
    qShape,
    tShape,
    qShapeMethod,
    tShapeMethod,
    qShapeConversion,
    tShapeConversion,
    mode,
    model,
    acc,
    intLenMax,
    intLoopMax,
    accW,
    accL,
    energy,
    energyVRNA,
    energyNoDangles,
    accNoLP,
    accNoGUend,
    energyAdd,
    temperature,
    windowWidth,
    windowOverlap,
    helixMinBP,
    helixMaxBP,
    helixMaxIL,
    helixMinPu,
    helixMaxE,
    helixFullE,
    out,
    outSep,
    outMode,
    outNumber,
    outOverlap,
    outMaxE,
    outMinPu,
    outDeltaE,
    outBestSeedOnly,
    outNoLP,
    outNoGUend,
    outCsvCols,
    outCsvSort,
    outPerRegion,
    outPairwise,
    verbose,
    defaultLogFile,
    threads,
    version,
    personality,
    parameterFile,
    help,
    fullhelp,
};

struct RegisteredOption {
    OptionId id;
    OptionSpec spec;
    double minimum{};
    double maximum{};
};

#define OPTION_REQUIRED(identifier, long_name, short_name, group_name, value_name, default_value, summary, is_basic) \
    RegisteredOption{OptionId::identifier, {long_name, short_name, OptionGroup::group_name, \
        OptionValueMode::required, value_name, default_value, summary, OptionSupport::implemented, false, is_basic}}
#define OPTION_NUMBER(identifier, long_name, short_name, group_name, value_name, default_value, summary, is_basic, minimum_value, maximum_value) \
    RegisteredOption{OptionId::identifier, {long_name, short_name, OptionGroup::group_name, \
        OptionValueMode::required, value_name, default_value, summary, OptionSupport::implemented, false, is_basic}, \
        minimum_value, maximum_value}
#define OPTION_BOOLEAN(identifier, long_name, group_name, default_value, summary, is_basic) \
    RegisteredOption{OptionId::identifier, {long_name, '\0', OptionGroup::group_name, \
        OptionValueMode::optionalBoolean, "BOOL", default_value, summary, OptionSupport::implemented, false, is_basic}}
#define OPTION_FLAG(identifier, long_name, short_name, group_name, summary, is_basic) \
    RegisteredOption{OptionId::identifier, {long_name, short_name, OptionGroup::group_name, \
        OptionValueMode::none, "", "", summary, OptionSupport::implemented, false, is_basic}}

constexpr std::array registeredOptions{
    OPTION_REQUIRED(query, "query", 'q', query, "RNA|FILE|STDIN", "", "Query RNA sequence, FASTA file, or STDIN.", true),
    OPTION_REQUIRED(qId, "qId", '\0', query, "ID", "query", "Query identifier used for literal input.", false),
    OPTION_NUMBER(qIdxPos0, "qIdxPos0", '\0', query, "INDEX", "1", "External coordinate of the query 5' end.", false, -2'000'000'000.0, 2'000'000'000.0),
    OPTION_REQUIRED(qSet, "qSet", '\0', query, "RANGES", "", "One-based query record subset.", false),
    OPTION_REQUIRED(qAcc, "qAcc", '\0', query, "N|C|P|E", "C", "Query accessibility source.", false),
    OPTION_NUMBER(qAccW, "qAccW", '\0', query, "N", "150", "Query accessibility window (0 means global).", false, 0.0, 99'999.0),
    OPTION_NUMBER(qAccL, "qAccL", '\0', query, "N", "100", "Query accessibility base-pair span.", false, 0.0, 99'999.0),
    OPTION_REQUIRED(qAccConstr, "qAccConstr", '\0', query, "CONSTRAINT", "", "Query accessibility structure constraint.", false),
    OPTION_REQUIRED(qAccFile, "qAccFile", '\0', query, "FILE", "", "Query accessibility probability or ED table.", false),
    OPTION_NUMBER(qIntLenMax, "qIntLenMax", '\0', query, "N", "0", "Maximum query interaction-site length.", false, 0.0, 99'999.0),
    OPTION_NUMBER(qIntLoopMax, "qIntLoopMax", '\0', query, "N", "16", "Maximum query-side interior-loop unpaired bases.", false, 0.0, 30.0),
    OPTION_REQUIRED(qRegion, "qRegion", '\0', query, "RANGES", "", "Explicit query prediction regions in external coordinates.", false),
    OPTION_NUMBER(qRegionLenMax, "qRegionLenMax", '\0', query, "N", "0", "Automatic maximum query region length.", false, 0.0, 99'999.0),
    OPTION_NUMBER(qPfScale, "qPfScale", '\0', query, "FACTOR", "1.07", "Query partition-function scaling factor.", false, 1.0, 99'999.0),

    OPTION_REQUIRED(target, "target", 't', target, "RNA|FILE|STDIN", "", "Target RNA sequence, FASTA file, or STDIN.", true),
    OPTION_REQUIRED(tId, "tId", '\0', target, "ID", "target", "Target identifier used for literal input.", false),
    OPTION_NUMBER(tIdxPos0, "tIdxPos0", '\0', target, "INDEX", "1", "External coordinate of the target 5' end.", false, -2'000'000'000.0, 2'000'000'000.0),
    OPTION_REQUIRED(tSet, "tSet", '\0', target, "RANGES", "", "One-based target record subset.", false),
    OPTION_REQUIRED(tAcc, "tAcc", '\0', target, "N|C|P|E", "C", "Target accessibility source.", false),
    OPTION_NUMBER(tAccW, "tAccW", '\0', target, "N", "150", "Target accessibility window (0 means global).", false, 0.0, 99'999.0),
    OPTION_NUMBER(tAccL, "tAccL", '\0', target, "N", "100", "Target accessibility base-pair span.", false, 0.0, 99'999.0),
    OPTION_REQUIRED(tAccConstr, "tAccConstr", '\0', target, "CONSTRAINT", "", "Target accessibility structure constraint.", false),
    OPTION_REQUIRED(tAccFile, "tAccFile", '\0', target, "FILE", "", "Target accessibility probability or ED table.", false),
    OPTION_NUMBER(tIntLenMax, "tIntLenMax", '\0', target, "N", "0", "Maximum target interaction-site length.", false, 0.0, 99'999.0),
    OPTION_NUMBER(tIntLoopMax, "tIntLoopMax", '\0', target, "N", "10", "Maximum target-side interior-loop unpaired bases.", false, 0.0, 30.0),
    OPTION_REQUIRED(tRegion, "tRegion", '\0', target, "RANGES", "", "Explicit target prediction regions in external coordinates.", false),
    OPTION_NUMBER(tRegionLenMax, "tRegionLenMax", '\0', target, "N", "0", "Automatic maximum target region length.", false, 0.0, 99'999.0),
    OPTION_NUMBER(tPfScale, "tPfScale", '\0', target, "FACTOR", "1.07", "Target partition-function scaling factor.", false, 1.0, 99'999.0),

    OPTION_BOOLEAN(noSeed, "noSeed", seed, "false", "Disable the generic seed requirement.", true),
    OPTION_REQUIRED(seedTQ, "seedTQ", '\0', seed, "SEEDS", "", "Explicit target-query seed encodings.", false),
    OPTION_NUMBER(seedBP, "seedBP", '\0', seed, "N", "7", "Required seed base-pair count.", true, 2.0, 20.0),
    OPTION_NUMBER(seedMaxUP, "seedMaxUP", '\0', seed, "N", "0", "Maximum total unpaired seed bases.", false, 0.0, 20.0),
    OPTION_NUMBER(seedQMaxUP, "seedQMaxUP", '\0', seed, "N", "-1", "Maximum unpaired query seed bases (-1 inherits seedMaxUP).", false, -1.0, 20.0),
    OPTION_NUMBER(seedTMaxUP, "seedTMaxUP", '\0', seed, "N", "-1", "Maximum unpaired target seed bases (-1 inherits seedMaxUP).", false, -1.0, 20.0),
    OPTION_NUMBER(seedMaxE, "seedMaxE", '\0', seed, "ENERGY", "0", "Maximum total seed energy.", false, -999.0, 999.0),
    OPTION_NUMBER(seedMaxEhybrid, "seedMaxEhybrid", '\0', seed, "ENERGY", "999", "Maximum seed hybridization energy.", false, -999.0, 999.0),
    OPTION_NUMBER(seedMinPu, "seedMinPu", '\0', seed, "PROBABILITY", "0", "Minimum seed-region unpaired probability.", false, 0.0, 1.0),
    OPTION_REQUIRED(seedQRange, "seedQRange", '\0', seed, "RANGES", "", "Query ranges searched for seeds.", false),
    OPTION_REQUIRED(seedTRange, "seedTRange", '\0', seed, "RANGES", "", "Target ranges searched for seeds.", false),
    OPTION_BOOLEAN(seedNoGU, "seedNoGU", seed, "false", "Disallow GU pairs inside seeds.", false),
    OPTION_BOOLEAN(seedNoGUend, "seedNoGUend", seed, "false", "Disallow GU pairs at seed ends.", false),

    OPTION_REQUIRED(qShape, "qShape", '\0', shape, "FILE", "", "Query SHAPE reactivity input.", false),
    OPTION_REQUIRED(tShape, "tShape", '\0', shape, "FILE", "", "Target SHAPE reactivity input.", false),
    OPTION_REQUIRED(qShapeMethod, "qShapeMethod", '\0', shape, "METHOD", "", "Query SHAPE pseudo-energy method (D, Z, or W encoding).", false),
    OPTION_REQUIRED(tShapeMethod, "tShapeMethod", '\0', shape, "METHOD", "", "Target SHAPE pseudo-energy method (D, Z, or W encoding).", false),
    OPTION_REQUIRED(qShapeConversion, "qShapeConversion", '\0', shape, "CONVERSION", "", "Query SHAPE conversion (M, C, S, L, or O encoding).", false),
    OPTION_REQUIRED(tShapeConversion, "tShapeConversion", '\0', shape, "CONVERSION", "", "Target SHAPE conversion (M, C, S, L, or O encoding).", false),

    OPTION_REQUIRED(mode, "mode", 'm', interaction, "H|M|S", "H", "Prediction mode: heuristic, exact, or seed-only.", true),
    OPTION_REQUIRED(model, "model", '\0', interaction, "S|X|B|P", "X", "Interaction model: single-site, seed extension, helix blocks, or ensemble.", true),
    OPTION_REQUIRED(acc, "acc", '\0', interaction, "N|C", "C", "Set query and target accessibility together.", true),
    OPTION_NUMBER(intLenMax, "intLenMax", '\0', interaction, "N", "0", "Set both maximum interaction-site lengths.", false, 0.0, 99'999.0),
    OPTION_NUMBER(intLoopMax, "intLoopMax", '\0', interaction, "N", "10", "Set both maximum interior-loop sizes.", false, 0.0, 30.0),
    OPTION_NUMBER(accW, "accW", '\0', interaction, "N", "150", "Set both accessibility windows.", false, 0.0, 99'999.0),
    OPTION_NUMBER(accL, "accL", '\0', interaction, "N", "100", "Set both accessibility spans.", false, 0.0, 99'999.0),
    OPTION_REQUIRED(energy, "energy", 'e', interaction, "B|V", "V", "Base-pair or nearest-neighbor energy model.", true),
    OPTION_REQUIRED(energyVRNA, "energyVRNA", '\0', interaction, "MODEL|FILE", "Turner04", "Nearest-neighbor parameter set or file.", false),
    OPTION_BOOLEAN(energyNoDangles, "energyNoDangles", interaction, "false", "Disable dangling-end contributions.", false),
    OPTION_BOOLEAN(accNoLP, "accNoLP", interaction, "false", "Disallow lonely pairs in accessibility folding.", false),
    OPTION_BOOLEAN(accNoGUend, "accNoGUend", interaction, "false", "Disallow GU helix ends in accessibility folding.", false),
    OPTION_NUMBER(energyAdd, "energyAdd", '\0', interaction, "ENERGY", "0", "Additive interaction energy correction.", false, -999.0, 999.0),
    OPTION_NUMBER(temperature, "temperature", '\0', interaction, "CELSIUS", "37", "Folding temperature in Celsius.", false, 0.0, 100.0),
    OPTION_NUMBER(windowWidth, "windowWidth", '\0', interaction, "N", "0", "Prediction-window width (0 disables decomposition).", false, 0.0, 99'999.0),
    OPTION_NUMBER(windowOverlap, "windowOverlap", '\0', interaction, "N", "150", "Overlap between prediction windows.", false, 0.0, 99'999.0),

    OPTION_NUMBER(helixMinBP, "helixMinBP", '\0', helix, "N", "2", "Minimum helix-block base-pair count.", false, 2.0, 4.0),
    OPTION_NUMBER(helixMaxBP, "helixMaxBP", '\0', helix, "N", "10", "Maximum helix-block base-pair count.", false, 2.0, 20.0),
    OPTION_NUMBER(helixMaxIL, "helixMaxIL", '\0', helix, "N", "0", "Maximum internal-loop size inside helix blocks.", false, 0.0, 2.0),
    OPTION_NUMBER(helixMinPu, "helixMinPu", '\0', helix, "PROBABILITY", "0", "Minimum helix-region unpaired probability.", false, 0.0, 1.0),
    OPTION_NUMBER(helixMaxE, "helixMaxE", '\0', helix, "ENERGY", "0", "Maximum accepted helix-block energy.", false, -999.0, 999.0),
    OPTION_BOOLEAN(helixFullE, "helixFullE", helix, "false", "Use full helix energy for helixMaxE.", false),

    RegisteredOption{OptionId::out, {"out", '\0', OptionGroup::output, OptionValueMode::required,
        "DEST", "STDOUT",
        "Primary output or qMinE/qSpotProb/qAcc/qPu/tMinE/tSpotProb/tAcc/tPu/"
        "pMinE/spotProb prefixed auxiliary destination.",
        OptionSupport::implemented, true, true}},
    OPTION_REQUIRED(outSep, "outSep", '\0', output, "CHAR", ";", "Tabular output separator.", false),
    OPTION_REQUIRED(outMode, "outMode", '\0', output, "N|D|C|E", "N", "Normal, detailed, CSV, or ensemble output.", true),
    OPTION_NUMBER(outNumber, "outNumber", 'n', output, "N", "1", "Maximum number of reported interactions.", true, 0.0, 1'000.0),
    OPTION_REQUIRED(outOverlap, "outOverlap", '\0', output, "N|T|Q|B", "B", "Allowed target/query overlap for suboptimal output.", true),
    OPTION_NUMBER(outMaxE, "outMaxE", '\0', output, "ENERGY", "0", "Report interactions below this energy.", false, -999.0, 999.0),
    OPTION_NUMBER(outMinPu, "outMinPu", '\0', output, "PROBABILITY", "0", "Minimum per-position unpaired probability.", false, 0.0, 1.0),
    OPTION_NUMBER(outDeltaE, "outDeltaE", '\0', output, "ENERGY", "100", "Suboptimal energy range above the optimum.", false, 0.0, 100.0),
    RegisteredOption{OptionId::outBestSeedOnly, {"outBestSeedOnly", '\0', OptionGroup::output,
        OptionValueMode::optionalBoolean, "BOOL", "false",
        "Report only the lowest-energy seed of each interaction.",
        OptionSupport::implemented, false, false}},
    OPTION_BOOLEAN(outNoLP, "outNoLP", output, "false", "Disallow lonely intermolecular pairs.", false),
    OPTION_BOOLEAN(outNoGUend, "outNoGUend", output, "false", "Disallow GU pairs at interaction helix ends.", false),
    OPTION_REQUIRED(outCsvCols, "outCsvCols", '\0', output, "COLUMNS", "id1,start1,end1,id2,start2,end2,subseqDP,hybridDP,E", "CSV columns; use '*' or an empty value for every documented column.", true),
    OPTION_REQUIRED(outCsvSort, "outCsvSort", '\0', output, "COLUMN", "", "CSV sort-column identifier.", false),
    OPTION_BOOLEAN(outPerRegion, "outPerRegion", output, "false", "Report each query-target region pair independently.", false),
    OPTION_BOOLEAN(outPairwise, "outPairwise", output, "false", "Predict corresponding query-target records only.", false),
    RegisteredOption{OptionId::verbose, {"verbose", 'v', OptionGroup::output, OptionValueMode::none,
        "", "false", "Record a verbose-logging request; the standalone core remains quiet.",
        OptionSupport::compatibilityOnly, false, false}},
    RegisteredOption{OptionId::defaultLogFile, {"default-log-file", '\0', OptionGroup::output,
        OptionValueMode::required, "FILE", "", "Record the compatibility log destination; no informational log is emitted.",
        OptionSupport::compatibilityOnly, false, false}},

    OPTION_NUMBER(threads, "threads", '\0', general, "N", "1", "Worker count (0 uses available hardware).", true, 0.0, 12.0),
    OPTION_FLAG(version, "version", '\0', general, "Show version information.", true),
    OPTION_REQUIRED(personality, "personality", '\0', general, "NAME", "IntaRNA", "Preset: IntaRNA/3/1/2/exact/helix/duplex/sTar/seed/ens; executable basenames are aliases.", true),
    RegisteredOption{OptionId::parameterFile, {"parameterFile", '\0', OptionGroup::general,
        OptionValueMode::required, "FILE", "", "Read registry options from a key=value file.",
        OptionSupport::implemented, true, true}},
    OPTION_FLAG(help, "help", 'h', general, "Show basic options.", true),
    OPTION_FLAG(fullhelp, "fullhelp", '\0', general, "Show every registered option.", true),
};

#undef OPTION_REQUIRED
#undef OPTION_NUMBER
#undef OPTION_BOOLEAN
#undef OPTION_FLAG

struct Assignment {
    const RegisteredOption* option{};
    std::string value;
    bool hasValue{};
    ConfigOrigin origin;
};

[[nodiscard]] auto equalInsensitive(
    const std::string_view left,
    const std::string_view right) noexcept -> bool {
    if (left.size() != right.size()) return false;
    for (std::size_t index = 0U; index < left.size(); ++index) {
        const auto leftCharacter = static_cast<unsigned char>(left[index]);
        const auto rightCharacter = static_cast<unsigned char>(right[index]);
        if (std::tolower(leftCharacter) != std::tolower(rightCharacter)) return false;
    }
    return true;
}

[[nodiscard]] auto lowerAscii(std::string value) -> std::string {
    std::ranges::transform(value, value.begin(), [](const unsigned char character) {
        return static_cast<char>(std::tolower(character));
    });
    return value;
}

[[nodiscard]] auto trim(std::string value) -> std::string {
    const auto first = std::find_if_not(value.begin(), value.end(), [](const unsigned char character) {
        return std::isspace(character) != 0;
    });
    const auto last = std::find_if_not(value.rbegin(), value.rend(), [](const unsigned char character) {
        return std::isspace(character) != 0;
    }).base();
    return first < last ? std::string(first, last) : std::string{};
}

[[nodiscard]] auto unquote(std::string value) -> std::string {
    if (value.size() >= 2U &&
        ((value.front() == '"' && value.back() == '"') ||
         (value.front() == '\'' && value.back() == '\''))) {
        return value.substr(1U, value.size() - 2U);
    }
    return value;
}

[[nodiscard]] auto originText(const ConfigOrigin& origin) -> std::string {
    switch (origin.source) {
        case ConfigSource::baseline:
            return origin.detail.empty() ? "built-in default" : origin.detail;
        case ConfigSource::personality:
            return origin.detail;
        case ConfigSource::parameterFile:
            return "parameter file '" + origin.detail + "' line " + std::to_string(origin.position);
        case ConfigSource::commandLine:
            return "command line argument " + std::to_string(origin.position);
    }
    return "unknown source";
}

[[nodiscard]] auto findLongOption(const std::string_view name) noexcept
    -> const RegisteredOption* {
    const auto found = std::ranges::find_if(registeredOptions, [name](const RegisteredOption& option) {
        return equalInsensitive(option.spec.longName, name);
    });
    return found == registeredOptions.end() ? nullptr : &*found;
}

[[nodiscard]] auto findShortOption(const char name) noexcept -> const RegisteredOption* {
    const auto normalized = static_cast<char>(
        std::tolower(static_cast<unsigned char>(name)));
    const auto found = std::ranges::find_if(registeredOptions, [normalized](const RegisteredOption& option) {
        return option.spec.shortName != '\0' &&
               std::tolower(static_cast<unsigned char>(option.spec.shortName)) == normalized;
    });
    return found == registeredOptions.end() ? nullptr : &*found;
}

[[nodiscard]] auto looksLikeOption(const std::string_view value) -> bool {
    return value.starts_with("--") ||
           (value.size() == 2U && value.front() == '-' &&
            std::isalpha(static_cast<unsigned char>(value.back())) != 0);
}

[[nodiscard]] auto isBooleanText(const std::string_view value) -> bool {
    return equalInsensitive(value, "1") || equalInsensitive(value, "0") ||
           equalInsensitive(value, "true") || equalInsensitive(value, "false") ||
           equalInsensitive(value, "yes") || equalInsensitive(value, "no") ||
           equalInsensitive(value, "on") || equalInsensitive(value, "off");
}

[[nodiscard]] auto tokenizeCommandLine(std::span<const std::string_view> arguments)
    -> std::expected<std::vector<Assignment>, std::string> {
    std::vector<Assignment> result;
    for (std::size_t index = 0U; index < arguments.size(); ++index) {
        const auto argumentPosition = index + 1U;
        std::string token(arguments[index]);
        if (token.empty()) {
            return std::unexpected("empty positional argument at command line argument " +
                                   std::to_string(argumentPosition));
        }

        const RegisteredOption* option{};
        std::string inlineValue;
        bool hasInlineValue{};
        if (token.starts_with("--")) {
            token.erase(0U, 2U);
            if (const auto equals = token.find('='); equals != std::string::npos) {
                inlineValue = token.substr(equals + 1U);
                token.erase(equals);
                hasInlineValue = true;
            }
            option = findLongOption(token);
        } else if (token.size() == 2U && token.front() == '-') {
            option = findShortOption(token.back());
        } else {
            return std::unexpected("unexpected positional argument '" + token +
                                   "' at command line argument " +
                                   std::to_string(argumentPosition));
        }
        if (option == nullptr) {
            return std::unexpected("unknown option '" + std::string(arguments[index]) +
                                   "' at command line argument " +
                                   std::to_string(argumentPosition));
        }

        Assignment assignment{option, std::move(inlineValue), hasInlineValue,
                              {ConfigSource::commandLine, "command line", argumentPosition}};
        if (option->spec.valueMode == OptionValueMode::none) {
            if (hasInlineValue) {
                return std::unexpected("--" + std::string(option->spec.longName) +
                                       " does not take a value at " +
                                       originText(assignment.origin));
            }
        } else if (option->spec.valueMode == OptionValueMode::required &&
                   !assignment.hasValue) {
            if (index + 1U >= arguments.size() || looksLikeOption(arguments[index + 1U])) {
                return std::unexpected("--" + std::string(option->spec.longName) +
                                       " requires a value at " +
                                       originText(assignment.origin));
            }
            assignment.value = std::string(arguments[++index]);
            assignment.hasValue = true;
        } else if (option->spec.valueMode == OptionValueMode::optionalBoolean &&
                   !assignment.hasValue && index + 1U < arguments.size() &&
                   isBooleanText(arguments[index + 1U])) {
            assignment.value = std::string(arguments[++index]);
            assignment.hasValue = true;
        }
        result.push_back(std::move(assignment));
    }
    return result;
}

[[nodiscard]] auto parseParameterLine(
    std::string line,
    const std::filesystem::path& path,
    const std::size_t lineNumber) -> std::expected<Assignment, std::string> {
    const ConfigOrigin origin{ConfigSource::parameterFile, path.string(), lineNumber};
    if (line.starts_with("--")) line.erase(0U, 2U);

    std::string key;
    std::string value;
    bool hasValue{};
    if (const auto equals = line.find('='); equals != std::string::npos) {
        key = trim(line.substr(0U, equals));
        value = unquote(trim(line.substr(equals + 1U)));
        hasValue = true;
    } else {
        const auto separator = std::find_if(line.begin(), line.end(), [](const unsigned char character) {
            return std::isspace(character) != 0;
        });
        key = std::string(line.begin(), separator);
        if (separator != line.end()) {
            value = unquote(trim(std::string(separator, line.end())));
            hasValue = !value.empty();
        }
    }
    if (key.empty()) {
        return std::unexpected("empty option at " + originText(origin));
    }

    const RegisteredOption* option{};
    if (key.size() == 2U && key.front() == '-') {
        option = findShortOption(key.back());
    } else {
        if (key.starts_with("-")) key.erase(0U, 1U);
        option = findLongOption(key);
        if (option == nullptr && key.size() == 1U) option = findShortOption(key.front());
    }
    if (option == nullptr) {
        return std::unexpected("unknown option '--" + key + "' at " + originText(origin));
    }
    if (option->spec.valueMode == OptionValueMode::none && hasValue) {
        return std::unexpected("--" + std::string(option->spec.longName) +
                               " does not take a value at " + originText(origin));
    }
    if (option->spec.valueMode == OptionValueMode::required && !hasValue) {
        return std::unexpected("--" + std::string(option->spec.longName) +
                               " requires a value at " + originText(origin));
    }
    return Assignment{option, std::move(value), hasValue, origin};
}

[[nodiscard]] auto normalizedPath(const std::filesystem::path& path)
    -> std::filesystem::path {
    std::error_code error;
    auto absolute = std::filesystem::absolute(path, error);
    return (error ? path : absolute).lexically_normal();
}

[[nodiscard]] auto readParameterFile(
    const std::filesystem::path& requestedPath,
    std::vector<std::filesystem::path> includeStack = {})
    -> std::expected<std::vector<Assignment>, std::string> {
    const auto path = normalizedPath(requestedPath);
    if (std::ranges::find(includeStack, path) != includeStack.end()) {
        std::ostringstream message;
        message << "parameter-file include cycle:";
        for (const auto& entry : includeStack) message << " '" << entry.string() << "' ->";
        message << " '" << path.string() << "'";
        return std::unexpected(message.str());
    }
    includeStack.push_back(path);

    std::ifstream input(path);
    if (!input) {
        return std::unexpected("cannot open parameter file '" + path.string() + "'");
    }

    std::vector<Assignment> result;
    std::string line;
    std::size_t lineNumber{};
    while (std::getline(input, line)) {
        ++lineNumber;
        if (const auto comment = line.find('#'); comment != std::string::npos) line.erase(comment);
        line = trim(std::move(line));
        if (line.empty()) continue;

        auto parsed = parseParameterLine(std::move(line), path, lineNumber);
        if (!parsed) return std::unexpected(parsed.error());
        if (parsed->option->id == OptionId::parameterFile) {
            auto nestedPath = std::filesystem::path(parsed->value);
            if (nestedPath.is_relative()) nestedPath = path.parent_path() / nestedPath;
            nestedPath = normalizedPath(nestedPath);
            parsed->value = nestedPath.string();
            result.push_back(*parsed);
            auto nested = readParameterFile(nestedPath, includeStack);
            if (!nested) {
                return std::unexpected(nested.error() + " (included from " +
                                       originText(parsed->origin) + ")");
            }
            result.insert(result.end(),
                          std::make_move_iterator(nested->begin()),
                          std::make_move_iterator(nested->end()));
        } else {
            result.push_back(std::move(*parsed));
        }
    }
    return result;
}

[[nodiscard]] auto validateDuplicates(const std::vector<Assignment>& assignments)
    -> std::expected<void, std::string> {
    using SeenKey = std::tuple<ConfigSource, std::string, std::string_view>;
    std::map<SeenKey, ConfigOrigin> seen;
    for (const auto& assignment : assignments) {
        if (assignment.option->spec.repeatable) continue;
        const auto key = SeenKey{
            assignment.origin.source,
            assignment.origin.detail,
            assignment.option->spec.longName,
        };
        const auto [position, inserted] = seen.emplace(key, assignment.origin);
        if (!inserted) {
            return std::unexpected(
                "duplicate --" + std::string(assignment.option->spec.longName) +
                ": first at " + originText(position->second) +
                ", again at " + originText(assignment.origin));
        }
    }
    return {};
}

[[nodiscard]] auto canonicalPersonality(
    const std::string_view value,
    const ConfigOrigin& origin) -> std::expected<std::string, std::string> {
    struct Name {
        std::string_view input;
        std::string_view canonical;
    };
    constexpr std::array names{
        Name{"default", "IntaRNA"},
        Name{"IntaRNAnew", "IntaRNA"},
        Name{"IntaRNA", "IntaRNA"},
        Name{"IntaRNA3", "IntaRNA3"},
        Name{"IntaRNA1", "IntaRNA1"},
        Name{"IntaRNA2", "IntaRNA2"},
        Name{"IntaRNAexact", "IntaRNAexact"},
        Name{"IntaRNAhelix", "IntaRNAhelix"},
        Name{"IntaRNAduplex", "IntaRNAduplex"},
        Name{"IntaRNAsTar", "IntaRNAsTar"},
        Name{"IntaRNAseed", "IntaRNAseed"},
        Name{"IntaRNAens", "IntaRNAens"},
    };
    const auto match = std::ranges::find_if(names, [value](const Name& name) {
        return equalInsensitive(name.input, value);
    });
    if (match == names.end()) {
        return std::unexpected("unknown personality '" + std::string(value) +
                               "' at " + originText(origin));
    }
    return std::string(match->canonical);
}

[[nodiscard]] auto executablePersonality(const std::string_view invocationName)
    -> std::expected<std::pair<std::string, ConfigOrigin>, std::string> {
    if (invocationName.empty()) {
        return std::pair{std::string{"IntaRNA"}, ConfigOrigin{
            ConfigSource::personality, "default personality", 0U}};
    }
    const auto separator = invocationName.find_last_of("/\\");
    auto basename = std::string(invocationName.substr(
        separator == std::string_view::npos ? 0U : separator + 1U));
    if (basename.size() > 4U &&
        equalInsensitive(std::string_view(basename).substr(basename.size() - 4U), ".exe")) {
        basename.resize(basename.size() - 4U);
    }
    const ConfigOrigin aliasOrigin{
        ConfigSource::personality, "executable alias '" + basename + "'", 0U};
    auto canonical = canonicalPersonality(basename, aliasOrigin);
    if (canonical) return std::pair{std::move(*canonical), aliasOrigin};
    if (!lowerAscii(basename).starts_with("intarna")) {
        return std::pair{std::string{"IntaRNA"}, ConfigOrigin{
            ConfigSource::personality, "default personality", 0U}};
    }
    return std::unexpected(canonical.error());
}

void markPersonalityDefaults(
    Config& config,
    const ConfigOrigin& selector,
    const std::string_view personality,
    const std::initializer_list<std::string_view> options) {
    const ConfigOrigin origin{
        ConfigSource::personality,
        "personality '" + std::string(personality) + "' selected by " + originText(selector),
        0U,
    };
    for (const auto option : options) {
        config.provenance[std::string(option)] = origin;
        config.assignmentHistory.push_back(
            {std::string(option), "<personality default>", origin});
    }
}

void applyPersonality(
    Config& config,
    const std::string_view personality,
    const ConfigOrigin& selector) {
    config.personality = std::string(personality);
    config.provenance["personality"] = selector;
    config.assignmentHistory.push_back(
        {"personality", std::string(personality), selector});

    if (equalInsensitive(personality, "IntaRNA1")) {
        config.model = InteractionModel::singleSite;
        config.mode = PredictionMode::heuristic;
        config.query.accessibilityWindow = 0U;
        config.target.accessibilityWindow = 0U;
        config.query.accessibilitySpan = 0U;
        config.target.accessibilitySpan = 0U;
        config.query.interactionLoopMax = 16U;
        config.target.interactionLoopMax = 16U;
        config.output.overlap = OverlapPolicy::query;
        markPersonalityDefaults(config, selector, personality,
            {"model", "mode", "qAccW", "tAccW", "qAccL", "tAccL",
             "qIntLoopMax", "tIntLoopMax", "outOverlap"});
    } else if (equalInsensitive(personality, "IntaRNA2")) {
        config.model = InteractionModel::singleSite;
        config.mode = PredictionMode::heuristic;
        config.query.interactionLoopMax = 16U;
        config.target.interactionLoopMax = 16U;
        config.output.overlap = OverlapPolicy::query;
        markPersonalityDefaults(config, selector, personality,
            {"model", "mode", "qIntLoopMax", "tIntLoopMax", "outOverlap"});
    } else if (equalInsensitive(personality, "IntaRNAexact")) {
        config.model = InteractionModel::seedExtension;
        config.mode = PredictionMode::exact;
        config.query.accessibilityWindow = 0U;
        config.target.accessibilityWindow = 0U;
        config.query.accessibilitySpan = 0U;
        config.target.accessibilitySpan = 0U;
        config.query.interactionLengthMax = 60U;
        config.target.interactionLengthMax = 60U;
        config.output.overlap = OverlapPolicy::both;
        markPersonalityDefaults(config, selector, personality,
            {"model", "mode", "qAccW", "tAccW", "qAccL", "tAccL",
             "qIntLenMax", "tIntLenMax", "outOverlap"});
    } else if (equalInsensitive(personality, "IntaRNAhelix")) {
        config.model = InteractionModel::helixBlocks;
        config.mode = PredictionMode::heuristic;
        markPersonalityDefaults(config, selector, personality, {"model", "mode"});
    } else if (equalInsensitive(personality, "IntaRNAduplex")) {
        config.query.accessibility = AccessibilityKind::disabled;
        config.target.accessibility = AccessibilityKind::disabled;
        markPersonalityDefaults(config, selector, personality, {"qAcc", "tAcc"});
    } else if (equalInsensitive(personality, "IntaRNAsTar")) {
        config.query.interactionLengthMax = 60U;
        config.target.interactionLengthMax = 60U;
        config.query.interactionLoopMax = 8U;
        config.target.interactionLoopMax = 8U;
        config.seed.noGu = true;
        config.seed.minUnpairedProbability = 0.001;
        config.output.minUnpairedProbability = 0.001;
        config.output.noLonelyPairs = true;
        config.output.noGuAtEnds = true;
        config.output.overlap = OverlapPolicy::query;
        config.output.mode = OutputMode::csv;
        config.output.csvColumns = "id1,id2,start1,end1,start2,end2,E";
        markPersonalityDefaults(config, selector, personality,
            {"qIntLenMax", "tIntLenMax", "qIntLoopMax", "tIntLoopMax",
             "seedNoGU", "seedMinPu", "outMinPu", "outNoLP", "outNoGUend",
             "outOverlap", "outMode", "outCsvCols"});
    } else if (equalInsensitive(personality, "IntaRNAseed")) {
        config.mode = PredictionMode::seedOnly;
        markPersonalityDefaults(config, selector, personality, {"mode"});
    } else if (equalInsensitive(personality, "IntaRNAens")) {
        config.model = InteractionModel::ensemble;
        markPersonalityDefaults(config, selector, personality, {"model"});
    }
}

template <typename Value>
[[nodiscard]] auto numericValue(const Assignment& assignment)
    -> std::expected<Value, std::string> {
    if (!assignment.hasValue) {
        return std::unexpected("--" + std::string(assignment.option->spec.longName) +
                               " requires a value at " + originText(assignment.origin));
    }
    Value value{};
    const auto* begin = assignment.value.data();
    const auto* end = begin + assignment.value.size();
    const auto [position, error] = std::from_chars(begin, end, value);
    if (error != std::errc{} || position != end ||
        (std::is_floating_point_v<Value> && !std::isfinite(value))) {
        return std::unexpected(
            "invalid numeric value '" + assignment.value + "' for --" +
            std::string(assignment.option->spec.longName) + " at " +
            originText(assignment.origin));
    }
    const auto converted = static_cast<long double>(value);
    if (converted < static_cast<long double>(assignment.option->minimum) ||
        converted > static_cast<long double>(assignment.option->maximum)) {
        std::ostringstream message;
        message << "value '" << assignment.value << "' for --"
                << assignment.option->spec.longName << " is outside ["
                << assignment.option->minimum << ',' << assignment.option->maximum
                << "] at " << originText(assignment.origin);
        return std::unexpected(message.str());
    }
    return value;
}

[[nodiscard]] auto booleanValue(const Assignment& assignment)
    -> std::expected<bool, std::string> {
    if (!assignment.hasValue) return true;
    if (equalInsensitive(assignment.value, "1") ||
        equalInsensitive(assignment.value, "true") ||
        equalInsensitive(assignment.value, "yes") ||
        equalInsensitive(assignment.value, "on")) {
        return true;
    }
    if (equalInsensitive(assignment.value, "0") ||
        equalInsensitive(assignment.value, "false") ||
        equalInsensitive(assignment.value, "no") ||
        equalInsensitive(assignment.value, "off")) {
        return false;
    }
    return std::unexpected(
        "invalid Boolean value '" + assignment.value + "' for --" +
        std::string(assignment.option->spec.longName) + " at " +
        originText(assignment.origin));
}

[[nodiscard]] auto accessibilityValue(const Assignment& assignment)
    -> std::expected<AccessibilityKind, std::string> {
    if (!assignment.hasValue || assignment.value.size() != 1U) {
        return std::unexpected("--" + std::string(assignment.option->spec.longName) +
                               " expects N, C, P, or E at " +
                               originText(assignment.origin));
    }
    switch (static_cast<char>(std::toupper(
        static_cast<unsigned char>(assignment.value.front())))) {
        case 'N': return AccessibilityKind::disabled;
        case 'C': return AccessibilityKind::compute;
        case 'P': return AccessibilityKind::probabilitiesFile;
        case 'E': return AccessibilityKind::energiesFile;
        default:
            return std::unexpected("--" + std::string(assignment.option->spec.longName) +
                                   " expects N, C, P, or E at " +
                                   originText(assignment.origin));
    }
}

void recordAssignment(Config& config, const Assignment& assignment) {
    const auto name = std::string(assignment.option->spec.longName);
    config.provenance[name] = assignment.origin;
    config.assignmentHistory.push_back(
        {name, assignment.hasValue ? assignment.value : "true", assignment.origin});
}

[[nodiscard]] auto applyAssignment(
    Config& config,
    const Assignment& assignment,
    bool& outputWasSet) -> std::expected<void, std::string> {
    const auto stringTo = [&](std::string& destination) -> std::expected<void, std::string> {
        if (!assignment.hasValue) {
            return std::unexpected("--" + std::string(assignment.option->spec.longName) +
                                   " requires a value at " + originText(assignment.origin));
        }
        destination = assignment.value;
        return {};
    };
    const auto booleanTo = [&](bool& destination) -> std::expected<void, std::string> {
        auto value = booleanValue(assignment);
        if (!value) return std::unexpected(value.error());
        destination = *value;
        return {};
    };
    const auto sizeTo = [&](std::size_t& destination) -> std::expected<void, std::string> {
        auto value = numericValue<std::size_t>(assignment);
        if (!value) return std::unexpected(value.error());
        destination = *value;
        return {};
    };
    const auto intTo = [&](int& destination) -> std::expected<void, std::string> {
        auto value = numericValue<int>(assignment);
        if (!value) return std::unexpected(value.error());
        destination = *value;
        return {};
    };
    const auto longTo = [&](long long& destination) -> std::expected<void, std::string> {
        auto value = numericValue<long long>(assignment);
        if (!value) return std::unexpected(value.error());
        destination = *value;
        return {};
    };
    const auto doubleTo = [&](double& destination) -> std::expected<void, std::string> {
        auto value = numericValue<double>(assignment);
        if (!value) return std::unexpected(value.error());
        destination = *value;
        return {};
    };
    const auto choice = [&](const std::string_view accepted)
        -> std::expected<char, std::string> {
        if (!assignment.hasValue || assignment.value.size() != 1U) {
            return std::unexpected("--" + std::string(assignment.option->spec.longName) +
                                   " expects one of " + std::string(accepted) + " at " +
                                   originText(assignment.origin));
        }
        const auto value = static_cast<char>(std::toupper(
            static_cast<unsigned char>(assignment.value.front())));
        if (accepted.find(value) == std::string_view::npos) {
            return std::unexpected("--" + std::string(assignment.option->spec.longName) +
                                   " expects one of " + std::string(accepted) + " at " +
                                   originText(assignment.origin));
        }
        return value;
    };

    std::expected<void, std::string> status{};
    switch (assignment.option->id) {
        case OptionId::query: status = stringTo(config.query.input); break;
        case OptionId::qId: status = stringTo(config.query.id); break;
        case OptionId::qIdxPos0: status = longTo(config.query.firstPosition); break;
        case OptionId::qSet: status = stringTo(config.query.subset); break;
        case OptionId::qAcc: {
            auto value = accessibilityValue(assignment);
            if (!value) status = std::unexpected(value.error());
            else config.query.accessibility = *value;
            break;
        }
        case OptionId::qAccW: status = sizeTo(config.query.accessibilityWindow); break;
        case OptionId::qAccL: status = sizeTo(config.query.accessibilitySpan); break;
        case OptionId::qAccConstr: status = stringTo(config.query.accessibilityConstraint); break;
        case OptionId::qAccFile: status = stringTo(config.query.accessibilityFile); break;
        case OptionId::qIntLenMax: status = sizeTo(config.query.interactionLengthMax); break;
        case OptionId::qIntLoopMax: status = sizeTo(config.query.interactionLoopMax); break;
        case OptionId::qRegion: status = stringTo(config.query.regions); break;
        case OptionId::qRegionLenMax: status = sizeTo(config.query.regionLengthMax); break;
        case OptionId::qPfScale: status = doubleTo(config.query.partitionScale); break;

        case OptionId::target: status = stringTo(config.target.input); break;
        case OptionId::tId: status = stringTo(config.target.id); break;
        case OptionId::tIdxPos0: status = longTo(config.target.firstPosition); break;
        case OptionId::tSet: status = stringTo(config.target.subset); break;
        case OptionId::tAcc: {
            auto value = accessibilityValue(assignment);
            if (!value) status = std::unexpected(value.error());
            else config.target.accessibility = *value;
            break;
        }
        case OptionId::tAccW: status = sizeTo(config.target.accessibilityWindow); break;
        case OptionId::tAccL: status = sizeTo(config.target.accessibilitySpan); break;
        case OptionId::tAccConstr: status = stringTo(config.target.accessibilityConstraint); break;
        case OptionId::tAccFile: status = stringTo(config.target.accessibilityFile); break;
        case OptionId::tIntLenMax: status = sizeTo(config.target.interactionLengthMax); break;
        case OptionId::tIntLoopMax: status = sizeTo(config.target.interactionLoopMax); break;
        case OptionId::tRegion: status = stringTo(config.target.regions); break;
        case OptionId::tRegionLenMax: status = sizeTo(config.target.regionLengthMax); break;
        case OptionId::tPfScale: status = doubleTo(config.target.partitionScale); break;

        case OptionId::noSeed: {
            auto value = booleanValue(assignment);
            if (!value) status = std::unexpected(value.error());
            else config.seed.required = !*value;
            break;
        }
        case OptionId::seedTQ: status = stringTo(config.seed.explicitSeeds); break;
        case OptionId::seedBP: status = sizeTo(config.seed.basePairs); break;
        case OptionId::seedMaxUP: status = sizeTo(config.seed.maxUnpaired); break;
        case OptionId::seedQMaxUP: status = intTo(config.seed.queryMaxUnpaired); break;
        case OptionId::seedTMaxUP: status = intTo(config.seed.targetMaxUnpaired); break;
        case OptionId::seedMaxE: status = doubleTo(config.seed.maxEnergy); break;
        case OptionId::seedMaxEhybrid: status = doubleTo(config.seed.maxHybridEnergy); break;
        case OptionId::seedMinPu: status = doubleTo(config.seed.minUnpairedProbability); break;
        case OptionId::seedQRange: status = stringTo(config.seed.queryRanges); break;
        case OptionId::seedTRange: status = stringTo(config.seed.targetRanges); break;
        case OptionId::seedNoGU: status = booleanTo(config.seed.noGu); break;
        case OptionId::seedNoGUend: status = booleanTo(config.seed.noGuAtEnds); break;

        case OptionId::qShape: status = stringTo(config.query.shapeFile); break;
        case OptionId::tShape: status = stringTo(config.target.shapeFile); break;
        case OptionId::qShapeMethod: status = stringTo(config.query.shapeMethod); break;
        case OptionId::tShapeMethod: status = stringTo(config.target.shapeMethod); break;
        case OptionId::qShapeConversion: status = stringTo(config.query.shapeConversion); break;
        case OptionId::tShapeConversion: status = stringTo(config.target.shapeConversion); break;

        case OptionId::mode: {
            auto value = choice("HMS");
            if (!value) status = std::unexpected(value.error());
            else if (*value == 'H') config.mode = PredictionMode::heuristic;
            else if (*value == 'M') config.mode = PredictionMode::exact;
            else config.mode = PredictionMode::seedOnly;
            break;
        }
        case OptionId::model: {
            auto value = choice("SXBP");
            if (!value) status = std::unexpected(value.error());
            else if (*value == 'S') config.model = InteractionModel::singleSite;
            else if (*value == 'X') config.model = InteractionModel::seedExtension;
            else if (*value == 'B') config.model = InteractionModel::helixBlocks;
            else config.model = InteractionModel::ensemble;
            break;
        }
        case OptionId::acc: {
            auto value = accessibilityValue(assignment);
            if (!value) status = std::unexpected(value.error());
            else if (*value == AccessibilityKind::probabilitiesFile ||
                     *value == AccessibilityKind::energiesFile) {
                status = std::unexpected("--acc accepts only N or C at " +
                                         originText(assignment.origin));
            } else {
                config.query.accessibility = *value;
                config.target.accessibility = *value;
            }
            break;
        }
        case OptionId::intLenMax: {
            auto value = numericValue<std::size_t>(assignment);
            if (!value) status = std::unexpected(value.error());
            else {
                config.query.interactionLengthMax = *value;
                config.target.interactionLengthMax = *value;
            }
            break;
        }
        case OptionId::intLoopMax: {
            auto value = numericValue<std::size_t>(assignment);
            if (!value) status = std::unexpected(value.error());
            else {
                config.query.interactionLoopMax = *value;
                config.target.interactionLoopMax = *value;
            }
            break;
        }
        case OptionId::accW: {
            auto value = numericValue<std::size_t>(assignment);
            if (!value) status = std::unexpected(value.error());
            else {
                config.query.accessibilityWindow = *value;
                config.target.accessibilityWindow = *value;
            }
            break;
        }
        case OptionId::accL: {
            auto value = numericValue<std::size_t>(assignment);
            if (!value) status = std::unexpected(value.error());
            else {
                config.query.accessibilitySpan = *value;
                config.target.accessibilitySpan = *value;
            }
            break;
        }
        case OptionId::energy: {
            auto value = choice("BV");
            if (!value) status = std::unexpected(value.error());
            else config.energy = *value == 'B' ? EnergyKind::basePair :
                                                 EnergyKind::nearestNeighbor;
            break;
        }
        case OptionId::energyVRNA: status = stringTo(config.energyParameters); break;
        case OptionId::energyNoDangles: status = booleanTo(config.noDangles); break;
        case OptionId::accNoLP: status = booleanTo(config.accessibilityNoLonelyPairs); break;
        case OptionId::accNoGUend: status = booleanTo(config.accessibilityNoGuAtEnds); break;
        case OptionId::energyAdd: status = doubleTo(config.additiveEnergy); break;
        case OptionId::temperature: status = doubleTo(config.temperatureCelsius); break;
        case OptionId::windowWidth: status = sizeTo(config.windowWidth); break;
        case OptionId::windowOverlap: status = sizeTo(config.windowOverlap); break;

        case OptionId::helixMinBP: status = sizeTo(config.helix.minBasePairs); break;
        case OptionId::helixMaxBP: status = sizeTo(config.helix.maxBasePairs); break;
        case OptionId::helixMaxIL: status = sizeTo(config.helix.maxInternalLoop); break;
        case OptionId::helixMinPu: status = doubleTo(config.helix.minUnpairedProbability); break;
        case OptionId::helixMaxE: status = doubleTo(config.helix.maxEnergy); break;
        case OptionId::helixFullE: status = booleanTo(config.helix.useFullEnergy); break;

        case OptionId::out:
            if (!assignment.hasValue) {
                status = std::unexpected("--out requires a value at " +
                                         originText(assignment.origin));
            } else {
                if (!outputWasSet) {
                    config.output.destinations.clear();
                    outputWasSet = true;
                }
                config.output.destinations.push_back(assignment.value);
            }
            break;
        case OptionId::outSep:
            if (!assignment.hasValue || assignment.value.size() != 1U) {
                status = std::unexpected("--outSep expects one character at " +
                                         originText(assignment.origin));
            } else {
                config.output.separator = assignment.value.front();
            }
            break;
        case OptionId::outMode: {
            auto value = choice("NDCE");
            if (!value) status = std::unexpected(value.error());
            else if (*value == 'N') config.output.mode = OutputMode::normal;
            else if (*value == 'D') config.output.mode = OutputMode::detailed;
            else if (*value == 'C') config.output.mode = OutputMode::csv;
            else config.output.mode = OutputMode::ensemble;
            break;
        }
        case OptionId::outNumber: status = sizeTo(config.output.number); break;
        case OptionId::outOverlap: {
            auto value = choice("NTQB");
            if (!value) status = std::unexpected(value.error());
            else if (*value == 'N') config.output.overlap = OverlapPolicy::neither;
            else if (*value == 'T') config.output.overlap = OverlapPolicy::target;
            else if (*value == 'Q') config.output.overlap = OverlapPolicy::query;
            else config.output.overlap = OverlapPolicy::both;
            break;
        }
        case OptionId::outMaxE: status = doubleTo(config.output.maxEnergy); break;
        case OptionId::outMinPu: status = doubleTo(config.output.minUnpairedProbability); break;
        case OptionId::outDeltaE: status = doubleTo(config.output.deltaEnergy); break;
        case OptionId::outBestSeedOnly: status = booleanTo(config.output.bestSeedOnly); break;
        case OptionId::outNoLP: status = booleanTo(config.output.noLonelyPairs); break;
        case OptionId::outNoGUend: status = booleanTo(config.output.noGuAtEnds); break;
        case OptionId::outCsvCols: status = stringTo(config.output.csvColumns); break;
        case OptionId::outCsvSort: status = stringTo(config.output.csvSort); break;
        case OptionId::outPerRegion: status = booleanTo(config.output.perRegion); break;
        case OptionId::outPairwise: status = booleanTo(config.output.pairwise); break;
        case OptionId::verbose: config.verbose = true; break;
        case OptionId::defaultLogFile: status = stringTo(config.logFile); break;

        case OptionId::threads: {
            auto value = numericValue<unsigned int>(assignment);
            if (!value) status = std::unexpected(value.error());
            else config.threads = *value;
            break;
        }
        case OptionId::version: config.action = RunAction::version; break;
        case OptionId::personality: {
            auto value = canonicalPersonality(assignment.value, assignment.origin);
            if (!value) status = std::unexpected(value.error());
            else config.personality = std::move(*value);
            break;
        }
        case OptionId::parameterFile:
            if (!assignment.hasValue) {
                status = std::unexpected("--parameterFile requires a value at " +
                                         originText(assignment.origin));
            } else {
                config.parameterFile = assignment.value;
                config.parameterFiles.push_back(assignment.value);
            }
            break;
        case OptionId::help: config.action = RunAction::help; break;
        case OptionId::fullhelp: config.action = RunAction::fullHelp; break;
    }

    if (!status) return std::unexpected(status.error());
    recordAssignment(config, assignment);
    return {};
}

[[nodiscard]] auto effectiveOrigin(
    const Config& config,
    const std::string_view option) -> ConfigOrigin {
    const auto found = config.provenance.find(option);
    return found == config.provenance.end() ? ConfigOrigin{} : found->second;
}

[[nodiscard]] auto citedOption(
    const Config& config,
    const std::string_view option) -> std::string {
    return "--" + std::string(option) + " (" +
           originText(effectiveOrigin(config, option)) + ")";
}

[[nodiscard]] auto validateShapeSide(
    const Config& config,
    const SideConfig& side,
    const std::string_view prefix) -> std::expected<void, std::string> {
    const auto shape = std::string(prefix) + "Shape";
    const auto method = std::string(prefix) + "ShapeMethod";
    const auto conversion = std::string(prefix) + "ShapeConversion";
    if (side.shapeFile.empty() &&
        (!side.shapeMethod.empty() || !side.shapeConversion.empty())) {
        const auto supplied = !side.shapeMethod.empty() ? method : conversion;
        return std::unexpected(citedOption(config, supplied) + " requires --" + shape);
    }
    if (!side.shapeFile.empty() && side.accessibility != AccessibilityKind::compute) {
        return std::unexpected(citedOption(config, shape) +
                               " requires computed accessibility (--" +
                               std::string(prefix) + "Acc=C)");
    }
    if (!side.shapeFile.empty()) {
        try {
            validateShapeEncoding(side.shapeMethod, side.shapeConversion);
        } catch (const std::invalid_argument& error) {
            const auto supplied = !side.shapeMethod.empty() ? method : conversion;
            return std::unexpected(citedOption(config, supplied) + ": " + error.what());
        }
    }
    return {};
}

[[nodiscard]] auto effectiveInteractionLength(const SideConfig& side) -> std::size_t {
    const auto interaction = side.interactionLengthMax;
    if (side.accessibility != AccessibilityKind::compute) return interaction;
    const auto accessibility = side.accessibilityWindow;
    if (interaction == 0U) return accessibility;
    if (accessibility == 0U) return interaction;
    return std::min(interaction, accessibility);
}

[[nodiscard]] auto hasUnsafeWindowTracker(const Config& config) -> bool {
    for (const auto& destination : config.output.destinations) {
        const auto separator = destination.find(':');
        const auto prefix = lowerAscii(destination.substr(0U, separator));
        if (prefix == "qspotprob" || prefix == "tspotprob" || prefix == "spotprob") {
            return true;
        }
    }
    return false;
}

[[nodiscard]] auto requestsEnsembleCsvColumn(const std::string_view columns) -> bool {
    std::size_t begin{};
    while (begin <= columns.size()) {
        const auto comma = columns.find(',', begin);
        const auto end = comma == std::string_view::npos ? columns.size() : comma;
        const auto column = trim(std::string(columns.substr(begin, end - begin)));
        if (equalInsensitive(column, "w") ||
            equalInsensitive(column, "Eall") ||
            equalInsensitive(column, "EallTotal") ||
            equalInsensitive(column, "Zall") ||
            equalInsensitive(column, "P_E")) {
            return true;
        }
        if (comma == std::string_view::npos) break;
        begin = comma + 1U;
    }
    return false;
}

[[nodiscard]] auto requestsSpotProbability(const Config& config) -> bool {
    for (const auto& destination : config.output.destinations) {
        const auto separator = destination.find(':');
        const auto prefix = destination.substr(0U, separator);
        if (equalInsensitive(prefix, "qSpotProb") ||
            equalInsensitive(prefix, "tSpotProb") ||
            equalInsensitive(prefix, "spotProb")) {
            return true;
        }
    }
    return false;
}

[[nodiscard]] auto validate(const Config& config) -> std::expected<void, std::string> {
    if (config.action != RunAction::predict) return {};
    if (config.query.input.empty()) return std::unexpected("missing required --query input");
    if (config.target.input.empty()) return std::unexpected("missing required --target input");
    if ((equalInsensitive(config.query.input, "STDIN") || config.query.input == "-") &&
        (equalInsensitive(config.target.input, "STDIN") || config.target.input == "-")) {
        return std::unexpected("query and target cannot both be read from standard input");
    }
    if (config.seed.required && config.seed.basePairs < 2U) {
        return std::unexpected("a required seed needs at least two base pairs");
    }
    if (config.mode == PredictionMode::seedOnly &&
        !config.seed.required && config.seed.explicitSeeds.empty()) {
        return std::unexpected(
            citedOption(config, "mode") +
            " selects seed-only prediction but " + citedOption(config, "noSeed") +
            " disables every seed and no --seedTQ was supplied");
    }
    if (config.model == InteractionModel::helixBlocks &&
        config.mode != PredictionMode::heuristic) {
        return std::unexpected(
            citedOption(config, "model") +
            " selects helix blocks, which require heuristic " +
            citedOption(config, "mode"));
    }
    if (config.model == InteractionModel::seedExtension &&
        (config.seed.required || !config.seed.explicitSeeds.empty()) &&
        (config.output.mode == OutputMode::ensemble ||
         requestsEnsembleCsvColumn(config.output.csvColumns) ||
         requestsSpotProbability(config))) {
        return std::unexpected(
            "seeded model X has no scientifically valid partition function; "
            "use model S/P or remove ensemble fields");
    }
    if (config.helix.minBasePairs > config.helix.maxBasePairs) {
        return std::unexpected(
            citedOption(config, "helixMinBP") + " cannot exceed " +
            citedOption(config, "helixMaxBP"));
    }
    if (!config.query.regions.empty() && config.query.regionLengthMax != 0U) {
        return std::unexpected(
            citedOption(config, "qRegion") + " conflicts with " +
            citedOption(config, "qRegionLenMax"));
    }
    if (!config.target.regions.empty() && config.target.regionLengthMax != 0U) {
        return std::unexpected(
            citedOption(config, "tRegion") + " conflicts with " +
            citedOption(config, "tRegionLenMax"));
    }
    if (auto status = validateShapeSide(config, config.query, "q"); !status) {
        return std::unexpected(status.error());
    }
    if (auto status = validateShapeSide(config, config.target, "t"); !status) {
        return std::unexpected(status.error());
    }
    if ((config.query.accessibility == AccessibilityKind::probabilitiesFile ||
         config.query.accessibility == AccessibilityKind::energiesFile) &&
        config.query.accessibilityFile.empty()) {
        return std::unexpected(citedOption(config, "qAcc") +
                               " requires --qAccFile");
    }
    if ((config.target.accessibility == AccessibilityKind::probabilitiesFile ||
         config.target.accessibility == AccessibilityKind::energiesFile) &&
        config.target.accessibilityFile.empty()) {
        return std::unexpected(citedOption(config, "tAcc") +
                               " requires --tAccFile");
    }
    if (config.windowWidth != 0U && config.windowWidth < 10U) {
        return std::unexpected(citedOption(config, "windowWidth") +
                               " has to be either zero or at least 10");
    }
    if (config.windowWidth != 0U && config.windowOverlap >= config.windowWidth) {
        return std::unexpected(citedOption(config, "windowOverlap") +
                               " has to be smaller than " +
                               citedOption(config, "windowWidth"));
    }
    if (config.windowWidth != 0U) {
        const auto queryLength = effectiveInteractionLength(config.query);
        const auto targetLength = effectiveInteractionLength(config.target);
        if (queryLength == 0U) {
            return std::unexpected(
                "windowed computation needs a finite query interaction length");
        }
        if (targetLength == 0U) {
            return std::unexpected(
                "windowed computation needs a finite target interaction length");
        }
        const auto requiredOverlap = std::max(queryLength, targetLength);
        if (config.windowOverlap < requiredOverlap) {
            return std::unexpected(
                citedOption(config, "windowOverlap") +
                " has to cover the maximal interaction length (at least " +
                std::to_string(requiredOverlap) + ")");
        }
        if (config.model == InteractionModel::ensemble ||
            config.output.mode == OutputMode::ensemble) {
            return std::unexpected(
                "ensemble/partition output is not scientifically composable across overlapping windows");
        }
        if (hasUnsafeWindowTracker(config)) {
            return std::unexpected(
                "spot-probability trackers are not scientifically composable across overlapping windows");
        }
    }
    return {};
}

[[nodiscard]] auto groupName(const OptionGroup group) -> std::string_view {
    switch (group) {
        case OptionGroup::query: return "Query";
        case OptionGroup::target: return "Target";
        case OptionGroup::seed: return "Seed";
        case OptionGroup::shape: return "SHAPE";
        case OptionGroup::interaction: return "Interaction";
        case OptionGroup::helix: return "Helix";
        case OptionGroup::output: return "Output";
        case OptionGroup::general: return "General";
    }
    return "Options";
}

[[nodiscard]] auto supportLabel(const OptionSupport support) -> std::string_view {
    switch (support) {
        case OptionSupport::implemented: return "";
        case OptionSupport::compatibilityOnly: return " [compatibility-only]";
        case OptionSupport::unavailable: return " [unavailable]";
    }
    return "";
}

} // namespace

auto Cli::parse(const std::span<const std::string_view> arguments)
    -> std::expected<Config, std::string> {
    return parse(arguments, {});
}

auto Cli::parse(
    const std::span<const std::string_view> arguments,
    const std::string_view invocationName) -> std::expected<Config, std::string> {
    auto commandLine = tokenizeCommandLine(arguments);
    if (!commandLine) return std::unexpected(commandLine.error());
    if (auto status = validateDuplicates(*commandLine); !status) {
        return std::unexpected(status.error());
    }

    std::vector<Assignment> parameterAssignments;
    for (const auto& assignment : *commandLine) {
        if (assignment.option->id != OptionId::parameterFile) continue;
        auto fromFile = readParameterFile(assignment.value);
        if (!fromFile) {
            return std::unexpected(fromFile.error() + " (requested at " +
                                   originText(assignment.origin) + ")");
        }
        parameterAssignments.insert(
            parameterAssignments.end(),
            std::make_move_iterator(fromFile->begin()),
            std::make_move_iterator(fromFile->end()));
    }
    if (auto status = validateDuplicates(parameterAssignments); !status) {
        return std::unexpected(status.error());
    }

    auto invocation = executablePersonality(invocationName);
    if (!invocation) return std::unexpected(invocation.error());
    auto selectedPersonality = invocation->first;
    auto personalityOrigin = invocation->second;
    for (const auto& assignment : parameterAssignments) {
        if (assignment.option->id != OptionId::personality) continue;
        auto selected = canonicalPersonality(assignment.value, assignment.origin);
        if (!selected) return std::unexpected(selected.error());
        selectedPersonality = std::move(*selected);
        personalityOrigin = assignment.origin;
    }
    for (const auto& assignment : *commandLine) {
        if (assignment.option->id != OptionId::personality) continue;
        auto selected = canonicalPersonality(assignment.value, assignment.origin);
        if (!selected) return std::unexpected(selected.error());
        selectedPersonality = std::move(*selected);
        personalityOrigin = assignment.origin;
    }

    Config config;
    for (const auto& option : registeredOptions) {
        config.provenance.emplace(
            std::string(option.spec.longName),
            ConfigOrigin{ConfigSource::baseline, "built-in default", 0U});
    }
    applyPersonality(config, selectedPersonality, personalityOrigin);

    bool outputWasSet{};
    for (const auto& assignment : parameterAssignments) {
        if (auto status = applyAssignment(config, assignment, outputWasSet); !status) {
            return std::unexpected(status.error());
        }
    }
    for (const auto& assignment : *commandLine) {
        if (auto status = applyAssignment(config, assignment, outputWasSet); !status) {
            return std::unexpected(status.error());
        }
    }

    if (config.threads == 0U) {
        config.threads = std::max(1U, std::thread::hardware_concurrency());
    }
    if (auto status = validate(config); !status) return std::unexpected(status.error());
    return config;
}

auto Cli::optionRegistry() noexcept -> std::span<const OptionSpec> {
    static const std::vector<OptionSpec> publicOptions = [] {
        std::vector<OptionSpec> result;
        result.reserve(registeredOptions.size());
        for (const auto& option : registeredOptions) result.push_back(option.spec);
        return result;
    }();
    return publicOptions;
}

auto Cli::version() -> std::string {
    return "IntaRNAnew 4.0.0\nstandalone C++23 clean-room implementation\n";
}

auto Cli::help(const bool full) -> std::string {
    std::ostringstream output;
    output << "IntaRNAnew predicts RNA-RNA interactions using a standalone C++23 engine.\n"
           << "Option names and documented choice values are ASCII case-insensitive.\n"
           << "Precedence: built-in defaults < personality < parameter file(s) < command line.\n\n";

    constexpr std::array groups{
        OptionGroup::query,
        OptionGroup::target,
        OptionGroup::seed,
        OptionGroup::shape,
        OptionGroup::interaction,
        OptionGroup::helix,
        OptionGroup::output,
        OptionGroup::general,
    };
    for (const auto group : groups) {
        const auto hasEntries = std::ranges::any_of(
            registeredOptions, [group, full](const RegisteredOption& option) {
                return option.spec.group == group && (full || option.spec.basic);
            });
        if (!hasEntries) continue;
        output << groupName(group) << ":\n";
        for (const auto& option : registeredOptions) {
            if (option.spec.group != group || (!full && !option.spec.basic)) continue;
            std::ostringstream spelling;
            if (option.spec.shortName != '\0') {
                spelling << '-' << option.spec.shortName << ", ";
            } else {
                spelling << "    ";
            }
            spelling << "--" << option.spec.longName;
            if (option.spec.valueMode == OptionValueMode::required) {
                spelling << " <" << option.spec.valueName << '>';
            } else if (option.spec.valueMode == OptionValueMode::optionalBoolean) {
                spelling << "[=BOOL]";
            }
            output << "  " << std::left << std::setw(37) << spelling.str()
                   << option.spec.description << supportLabel(option.spec.support);
            if (!option.spec.defaultValue.empty()) {
                output << " (default: " << option.spec.defaultValue << ')';
            }
            if (option.spec.repeatable) output << " (repeatable)";
            output << '\n';
        }
        output << '\n';
    }
    if (!full) output << "Use --fullhelp to list all registered compatibility options.\n";
    if (full) {
        output << "Compatibility-only options are parsed, range-checked, recorded in Config, "
                  "and otherwise have no standalone runtime side effect.\n";
    }
    return output.str();
}

} // namespace intarnanew
