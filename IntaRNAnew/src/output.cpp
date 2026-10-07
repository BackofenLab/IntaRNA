#include "intarnanew/output.hpp"

#include "intarnanew/output_plan.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <map>
#include <optional>
#include <ranges>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

namespace intarnanew {
namespace {

const std::vector<std::string> allColumns{
    "id1", "id2", "seq1", "seq2", "subseq1", "subseq2", "subseqDP", "subseqDB",
    "start1", "end1", "start2", "end2", "hybridDP", "hybridDB", "hybridDPfull",
    "hybridDBfull", "bpList", "E", "Etotal", "ED1", "ED2", "Pu1", "Pu2", "E_init",
    "E_loops", "E_dangleL", "E_dangleR", "E_endL", "E_endR", "E_hybrid", "E_norm",
    "E_hybridNorm", "E_add", "seedStart1", "seedEnd1", "seedStart2", "seedEnd2", "seedE",
    "seedED1", "seedED2", "seedPu1", "seedPu2", "w", "Eall", "Eall1", "Eall2", "Zall",
    "Zall1", "Zall2", "EallTotal", "P_E", "RT",
};

[[nodiscard]] auto number(const double value) -> std::string {
    if (!std::isfinite(value)) return "NA";
    if (std::abs(value) < 0.0000005) return "0";
    std::array<char, 64U> buffer{};
    const auto [end, error] = std::to_chars(
        buffer.data(), buffer.data() + buffer.size(), value, std::chars_format::general, 6);
    if (error == std::errc{}) return std::string(buffer.data(), end);
    std::ostringstream fallback;
    fallback.imbue(std::locale::classic());
    fallback << std::setprecision(6) << std::defaultfloat << value;
    return fallback.str();
}

[[nodiscard]] auto centikcal(const Energy value) noexcept -> Energy {
    if (!std::isfinite(value)) return value;
    const auto scaled = value * 100.0;
    const auto tolerance = 64.0 * std::numeric_limits<double>::epsilon() *
                           std::max(1.0, std::abs(scaled));
    return std::trunc(scaled + std::copysign(tolerance, scaled)) / 100.0;
}

[[nodiscard]] auto energyNumber(const Energy value) -> std::string {
    return number(centikcal(value));
}

[[nodiscard]] auto fixedEnergyNumber(const Energy value) -> std::string {
    if (!std::isfinite(value)) return "NA";
    std::ostringstream output;
    output.imbue(std::locale::classic());
    const auto rounded = centikcal(value);
    output << std::fixed << std::setprecision(2) << (rounded == 0.0 ? 0.0 : rounded);
    return output.str();
}

[[nodiscard]] auto scientificNumber(const double value) -> std::string {
    if (!std::isfinite(value)) return "NA";
    std::ostringstream output;
    output.imbue(std::locale::classic());
    output << std::scientific << std::setprecision(6) << value;
    return output.str();
}

[[nodiscard]] auto combinedEnergy(
    const Energy interactionEnergy,
    const PredictionResult& result) noexcept -> Energy {
    if (!std::isfinite(interactionEnergy) ||
        !std::isfinite(result.targetEnsembleFreeEnergy) ||
        !std::isfinite(result.queryEnsembleFreeEnergy)) {
        return infinity;
    }
    return interactionEnergy + result.targetEnsembleFreeEnergy +
           result.queryEnsembleFreeEnergy;
}

[[nodiscard]] auto totalEnsembleEnergy(const PredictionResult& result) noexcept -> Energy {
    return combinedEnergy(result.ensembleFreeEnergy, result);
}

[[nodiscard]] auto partitionWeight(const double logPartition) noexcept -> double {
    return std::isfinite(logPartition) ? std::exp(logPartition) : infinity;
}

[[nodiscard]] auto lowerAscii(std::string_view value) -> std::string {
    std::string result(value);
    std::ranges::transform(result, result.begin(), [](const unsigned char character) {
        return static_cast<char>(std::tolower(character));
    });
    return result;
}

[[nodiscard]] auto split(const std::string_view value, const char delimiter) -> std::vector<std::string> {
    std::vector<std::string> result;
    std::size_t start{};
    for (std::size_t index = 0; index <= value.size(); ++index) {
        if (index == value.size() || value[index] == delimiter) {
            result.emplace_back(value.substr(start, index - start));
            start = index + 1U;
        }
    }
    return result;
}

[[nodiscard]] auto patterns(
    const Interaction& interaction) -> std::pair<std::string, std::string> {
    const auto targetRange = interaction.targetRange();
    const auto queryRange = interaction.queryRange();
    std::string targetPattern(targetRange.size(), '.');
    std::string queryPattern(queryRange.size(), '.');
    for (const auto pair : interaction.pairs) {
        targetPattern[pair.target - targetRange.begin] = '(';
        queryPattern[pair.query - queryRange.begin] = ')';
    }
    return {std::move(targetPattern), std::move(queryPattern)};
}

[[nodiscard]] auto barPatterns(
    const Interaction& interaction,
    const Sequence& target,
    const Sequence& query,
    const bool full) -> std::pair<std::string, std::string> {
    const auto targetRange = full ? Interval{0U, target.size() - 1U} : interaction.targetRange();
    const auto queryRange = full ? Interval{0U, query.size() - 1U} : interaction.queryRange();
    std::string targetPattern(targetRange.size(), '.');
    std::string queryPattern(queryRange.size(), '.');
    for (const auto pair : interaction.pairs) {
        targetPattern[pair.target - targetRange.begin] = '|';
        queryPattern[pair.query - queryRange.begin] = '|';
    }
    return {std::to_string(target.externalIndex(targetRange.begin)) + targetPattern,
            std::to_string(query.externalIndex(queryRange.begin)) + queryPattern};
}

[[nodiscard]] auto selectedSeeds(
    const Interaction& interaction,
    const bool bestSeedOnly) -> std::vector<const SeedMatch*> {
    std::vector<const SeedMatch*> result;
    result.reserve(interaction.seeds.size());
    for (const auto& seed : interaction.seeds) result.push_back(&seed);
    std::ranges::sort(result, [](const SeedMatch* left, const SeedMatch* right) {
        return seedMatchLess(*left, *right);
    });
    if (bestSeedOnly && result.size() > 1U) result.resize(1U);
    return result;
}

template <typename Formatter>
[[nodiscard]] auto joinedSeedValues(
    const std::vector<const SeedMatch*>& seeds,
    const std::string_view delimiter,
    Formatter formatter) -> std::string {
    std::string result;
    for (const auto* seed : seeds) {
        if (!result.empty()) result += delimiter;
        result += formatter(*seed);
    }
    return result;
}

[[nodiscard]] auto csvColumns(const std::string& specification) -> std::expected<std::vector<std::string>, std::string> {
    auto result = specification.empty() || specification == "*" ? allColumns : split(specification, ',');
    const std::unordered_set<std::string> available(allColumns.begin(), allColumns.end());
    for (const auto& column : result) {
        if (!available.contains(column)) return std::unexpected("unknown CSV column '" + column + "'");
    }
    return result;
}

[[nodiscard]] auto columnValue(
    const std::string& column,
    const Interaction& interaction,
    const Sequence& target,
    const Sequence& query,
    const PredictionResult& result,
    const bool bestSeedOnly) -> std::string {
    const auto targetRange = interaction.targetRange();
    const auto queryRange = interaction.queryRange();
    const auto [targetPattern, queryPattern] = patterns(interaction);
    if (column == "id1") return interaction.targetId;
    if (column == "id2") return interaction.queryId;
    if (column == "seq1") return target.str();
    if (column == "seq2") return query.str();
    if (column == "subseq1") return target.str().substr(targetRange.begin, targetRange.size());
    if (column == "subseq2") return query.str().substr(queryRange.begin, queryRange.size());
    if (column == "subseqDP") return target.str().substr(targetRange.begin, targetRange.size()) + "&" +
                                      query.str().substr(queryRange.begin, queryRange.size());
    if (column == "subseqDB") {
        return std::to_string(target.externalIndex(targetRange.begin)) +
               target.str().substr(targetRange.begin, targetRange.size()) + "&" +
               std::to_string(query.externalIndex(queryRange.begin)) +
               query.str().substr(queryRange.begin, queryRange.size());
    }
    if (column == "start1") return std::to_string(target.externalIndex(targetRange.begin));
    if (column == "end1") return std::to_string(target.externalIndex(targetRange.end));
    if (column == "start2") return std::to_string(query.externalIndex(queryRange.begin));
    if (column == "end2") return std::to_string(query.externalIndex(queryRange.end));
    if (column == "hybridDP") return targetPattern + "&" + queryPattern;
    if (column == "hybridDB") {
        const auto bars = barPatterns(interaction, target, query, false);
        return bars.first + "&" + bars.second;
    }
    if (column == "hybridDPfull") {
        std::string targetFull(target.size(), '.');
        std::string queryFull(query.size(), '.');
        for (const auto pair : interaction.pairs) {
            targetFull[pair.target] = '(';
            queryFull[pair.query] = ')';
        }
        return targetFull + "&" + queryFull;
    }
    if (column == "hybridDBfull") {
        const auto bars = barPatterns(interaction, target, query, true);
        return bars.first + "&" + bars.second;
    }
    if (column == "bpList") {
        std::string value;
        for (const auto pair : interaction.pairs) {
            if (!value.empty()) value.push_back(':');
            value += "(" + std::to_string(target.externalIndex(pair.target)) + "," +
                     std::to_string(query.externalIndex(pair.query)) + ")";
        }
        return value;
    }
    if (column == "E") return energyNumber(interaction.energy.total());
    if (column == "Etotal") return energyNumber(combinedEnergy(interaction.energy.total(), result));
    if (column == "ED1") return energyNumber(interaction.energy.openingTarget);
    if (column == "ED2") return energyNumber(interaction.energy.openingQuery);
    if (column == "Pu1") return number(interaction.unpairedTarget);
    if (column == "Pu2") return number(interaction.unpairedQuery);
    if (column == "E_init") return energyNumber(interaction.energy.initiation);
    if (column == "E_loops") return energyNumber(interaction.energy.loops);
    if (column == "E_dangleL") return energyNumber(interaction.energy.dangleLeft);
    if (column == "E_dangleR") return energyNumber(interaction.energy.dangleRight);
    if (column == "E_endL") return energyNumber(interaction.energy.endLeft);
    if (column == "E_endR") return energyNumber(interaction.energy.endRight);
    if (column == "E_hybrid") return energyNumber(interaction.energy.hybrid());
    if (column == "E_norm") {
        const auto denominator = std::log(static_cast<double>(target.size())) +
                                 std::log(static_cast<double>(query.size()));
        return number(denominator > 0.0 ? interaction.energy.total() / denominator :
                      std::numeric_limits<double>::quiet_NaN());
    }
    if (column == "E_hybridNorm") {
        const auto denominator = std::log(static_cast<double>(target.size())) +
                                 std::log(static_cast<double>(query.size()));
        return number(denominator > 0.0 ? interaction.energy.hybrid() / denominator :
                      std::numeric_limits<double>::quiet_NaN());
    }
    if (column == "E_add") return energyNumber(interaction.energy.additive);
    if (column.starts_with("seed")) {
        const auto seeds = selectedSeeds(interaction, bestSeedOnly);
        if (seeds.empty()) return "NAN";
        return joinedSeedValues(seeds, ":", [&](const SeedMatch& seed) {
            const auto& first = interaction.pairs.at(seed.firstPair);
            const auto& last = interaction.pairs.at(seed.lastPair);
            if (column == "seedStart1") return std::to_string(target.externalIndex(first.target));
            if (column == "seedEnd1") return std::to_string(target.externalIndex(last.target));
            if (column == "seedStart2") return std::to_string(query.externalIndex(last.query));
            if (column == "seedEnd2") return std::to_string(query.externalIndex(first.query));
            if (column == "seedE") return energyNumber(seed.energy);
            if (column == "seedED1") return energyNumber(seed.openingTarget);
            if (column == "seedED2") return energyNumber(seed.openingQuery);
            if (column == "seedPu1") return number(seed.unpairedTarget);
            if (column == "seedPu2") return number(seed.unpairedQuery);
            return std::string{};
        });
    }
    if (column == "w") return number(std::exp(-centikcal(interaction.energy.total()) / result.rt));
    if (column == "Eall") return energyNumber(result.ensembleFreeEnergy);
    if (column == "Eall1") return energyNumber(result.targetEnsembleFreeEnergy);
    if (column == "Eall2") return energyNumber(result.queryEnsembleFreeEnergy);
    if (column == "Zall") return number(partitionWeight(result.logPartition));
    if (column == "Zall1") return number(partitionWeight(result.targetLogPartition));
    if (column == "Zall2") return number(partitionWeight(result.queryLogPartition));
    if (column == "EallTotal") return energyNumber(totalEnsembleEnergy(result));
    if (column == "P_E") {
        return std::isfinite(result.logPartition)
            ? number(interaction.probability) : "NA";
    }
    if (column == "RT") return number(result.rt);
    return "";
}

enum class CsvSortKind { text, integer, floating };

[[nodiscard]] auto csvSortKind(const std::string_view column) -> CsvSortKind {
    static const std::unordered_set<std::string_view> integerColumns{
        "start1", "end1", "start2", "end2",
        "seedStart1", "seedEnd1", "seedStart2", "seedEnd2",
    };
    static const std::unordered_set<std::string_view> floatingColumns{
        "E", "Etotal", "ED1", "ED2", "Pu1", "Pu2", "E_init", "E_loops",
        "E_dangleL", "E_dangleR", "E_endL", "E_endR", "E_hybrid", "E_norm",
        "E_hybridNorm", "E_add", "seedE", "seedED1", "seedED2", "seedPu1",
        "seedPu2", "w", "Eall", "Eall1", "Eall2", "Zall", "Zall1", "Zall2",
        "EallTotal", "P_E", "RT",
    };
    if (integerColumns.contains(column)) return CsvSortKind::integer;
    if (floatingColumns.contains(column)) return CsvSortKind::floating;
    return CsvSortKind::text;
}

[[nodiscard]] auto integerColumnValue(
    const std::string_view column,
    const Interaction& interaction,
    const Sequence& target,
    const Sequence& query) -> std::optional<long long> {
    const auto targetRange = interaction.targetRange();
    const auto queryRange = interaction.queryRange();
    if (column == "start1") return target.externalIndex(targetRange.begin);
    if (column == "end1") return target.externalIndex(targetRange.end);
    if (column == "start2") return query.externalIndex(queryRange.begin);
    if (column == "end2") return query.externalIndex(queryRange.end);
    const auto* seed = interaction.bestSeed();
    if (seed == nullptr) return std::nullopt;
    const auto& first = interaction.pairs[seed->firstPair];
    const auto& last = interaction.pairs[seed->lastPair];
    if (column == "seedStart1") return target.externalIndex(first.target);
    if (column == "seedEnd1") return target.externalIndex(last.target);
    if (column == "seedStart2") return query.externalIndex(last.query);
    if (column == "seedEnd2") return query.externalIndex(first.query);
    return std::nullopt;
}

[[nodiscard]] auto finite(const double value) -> std::optional<double> {
    return std::isfinite(value) ? std::optional<double>{value} : std::nullopt;
}

[[nodiscard]] auto floatingColumnValue(
    const std::string_view column,
    const Interaction& interaction,
    const Sequence& target,
    const Sequence& query,
    const PredictionResult& result) -> std::optional<double> {
    if (column == "E") return finite(interaction.energy.total());
    if (column == "Etotal") return finite(combinedEnergy(interaction.energy.total(), result));
    if (column == "ED1") return finite(interaction.energy.openingTarget);
    if (column == "ED2") return finite(interaction.energy.openingQuery);
    if (column == "Pu1") return finite(interaction.unpairedTarget);
    if (column == "Pu2") return finite(interaction.unpairedQuery);
    if (column == "E_init") return finite(interaction.energy.initiation);
    if (column == "E_loops") return finite(interaction.energy.loops);
    if (column == "E_dangleL") return finite(interaction.energy.dangleLeft);
    if (column == "E_dangleR") return finite(interaction.energy.dangleRight);
    if (column == "E_endL") return finite(interaction.energy.endLeft);
    if (column == "E_endR") return finite(interaction.energy.endRight);
    if (column == "E_hybrid") return finite(interaction.energy.hybrid());
    if (column == "E_norm" || column == "E_hybridNorm") {
        const auto denominator = std::log(static_cast<double>(target.size())) +
                                 std::log(static_cast<double>(query.size()));
        if (!(denominator > 0.0)) return std::nullopt;
        return finite((column == "E_norm" ? interaction.energy.total() : interaction.energy.hybrid()) /
                      denominator);
    }
    if (column == "E_add") return finite(interaction.energy.additive);
    if (column == "seedE") {
        const auto* seed = interaction.bestSeed();
        return seed != nullptr ? finite(seed->energy) : std::nullopt;
    }
    if (column == "seedED1") {
        const auto* seed = interaction.bestSeed();
        return seed != nullptr ? finite(seed->openingTarget) : std::nullopt;
    }
    if (column == "seedED2") {
        const auto* seed = interaction.bestSeed();
        return seed != nullptr ? finite(seed->openingQuery) : std::nullopt;
    }
    if (column == "seedPu1") {
        const auto* seed = interaction.bestSeed();
        return seed != nullptr ? finite(seed->unpairedTarget) : std::nullopt;
    }
    if (column == "seedPu2") {
        const auto* seed = interaction.bestSeed();
        return seed != nullptr ? finite(seed->unpairedQuery) : std::nullopt;
    }
    if (column == "w") {
        return finite(std::exp(-centikcal(interaction.energy.total()) / result.rt));
    }
    if (column == "Eall") return finite(result.ensembleFreeEnergy);
    if (column == "Eall1") return finite(result.targetEnsembleFreeEnergy);
    if (column == "Eall2") return finite(result.queryEnsembleFreeEnergy);
    if (column == "Zall") return finite(partitionWeight(result.logPartition));
    if (column == "Zall1") return finite(partitionWeight(result.targetLogPartition));
    if (column == "Zall2") return finite(partitionWeight(result.queryLogPartition));
    if (column == "EallTotal") return finite(totalEnsembleEnergy(result));
    if (column == "P_E") {
        return std::isfinite(result.logPartition)
            ? finite(interaction.probability) : std::nullopt;
    }
    if (column == "RT") return finite(result.rt);
    return std::nullopt;
}

template <typename Value>
[[nodiscard]] auto optionalLess(
    const std::optional<Value>& left,
    const std::optional<Value>& right) -> bool {
    if (!left) return false;
    if (!right) return true;
    return *left < *right;
}

[[nodiscard]] auto csvField(const std::string& value, const char separator) -> std::string {
    if (value.find_first_of(std::string{separator} + "\"\r\n") == std::string::npos) return value;
    std::string escaped;
    escaped.reserve(value.size() + 2U);
    escaped.push_back('"');
    for (const char character : value) {
        if (character == '"') escaped.push_back('"');
        escaped.push_back(character);
    }
    escaped.push_back('"');
    return escaped;
}

[[nodiscard]] auto csv(
    const Config& config,
    const Sequence& target,
    const Sequence& query,
    const PredictionResult& result,
    const bool includeHeader) -> std::expected<std::string, std::string> {
    auto columns = csvColumns(config.output.csvColumns);
    if (!columns) return std::unexpected(columns.error());
    auto rows = result.interactions;
    if (!config.output.csvSort.empty()) {
        if (std::ranges::find(*columns, config.output.csvSort) == columns->end()) {
            return std::unexpected("--outCsvSort column is not selected in --outCsvCols");
        }
        const auto kind = csvSortKind(config.output.csvSort);
        std::ranges::stable_sort(rows, [&](const Interaction& left, const Interaction& right) {
            if (kind == CsvSortKind::integer) {
                return optionalLess(
                    integerColumnValue(config.output.csvSort, left, target, query),
                    integerColumnValue(config.output.csvSort, right, target, query));
            }
            if (kind == CsvSortKind::floating) {
                return optionalLess(
                    floatingColumnValue(config.output.csvSort, left, target, query, result),
                    floatingColumnValue(config.output.csvSort, right, target, query, result));
            }
            return columnValue(config.output.csvSort, left, target, query, result,
                               config.output.bestSeedOnly) <
                   columnValue(config.output.csvSort, right, target, query, result,
                               config.output.bestSeedOnly);
        });
    }
    std::ostringstream output;
    output.imbue(std::locale::classic());
    if (includeHeader) {
        for (Index index = 0; index < columns->size(); ++index) {
            if (index != 0U) output << config.output.separator;
            output << csvField((*columns)[index], config.output.separator);
        }
        output << '\n';
    }
    for (const auto& interaction : rows) {
        for (Index index = 0; index < columns->size(); ++index) {
            if (index != 0U) output << config.output.separator;
            output << csvField(
                columnValue((*columns)[index], interaction, target, query, result,
                            config.output.bestSeedOnly),
                config.output.separator);
        }
        output << '\n';
    }
    return output.str();
}

[[nodiscard]] auto reverseString(std::string value) -> std::string {
    std::ranges::reverse(value);
    return value;
}

[[nodiscard]] auto leftFlank(std::string value) -> std::string {
    constexpr std::size_t flankWidth{10U};
    if (value.size() <= flankWidth) return value;
    return value.substr(0U, 3U) + "..." + value.substr(value.size() - 4U);
}

[[nodiscard]] auto rightFlank(std::string value) -> std::string {
    constexpr std::size_t flankWidth{10U};
    if (value.size() <= flankWidth) return value;
    return value.substr(0U, 4U) + "..." + value.substr(value.size() - 3U);
}

[[nodiscard]] auto annotatedIndex(const Sequence& sequence, const Index index) -> std::string {
    const auto external = sequence.externalIndex(index);
    if (sequence.firstPosition() < 0 && external > 0) {
        return "+" + std::to_string(external);
    }
    return std::to_string(external);
}

[[nodiscard]] auto coordinateLine(
    const std::size_t firstColumn,
    const std::size_t lastColumn,
    const std::string& first,
    const std::string& last) -> std::string {
    std::string line(firstColumn + 1U - std::min(firstColumn + 1U, first.size()), ' ');
    line += first;
    if (lastColumn != firstColumn) {
        if (line.size() < lastColumn) line.append(lastColumn - line.size(), ' ');
        line += last;
    }
    return line;
}

[[nodiscard]] auto markerLine(
    const std::size_t firstColumn,
    const std::size_t lastColumn) -> std::string {
    std::string line(firstColumn, ' ');
    line.push_back('|');
    if (lastColumn != firstColumn) {
        line.append(lastColumn - firstColumn - 1U, ' ');
        line.push_back('|');
    }
    return line;
}

struct DisplayPair {
    BasePair pair;
    Index originalIndex{};
};

struct InteractionChart {
    std::string targetBackbone;
    std::string targetPairs;
    std::string pairSymbols;
    std::string queryPairs;
    std::string queryBackbone;
    std::vector<DisplayPair> pairs;
};

[[nodiscard]] auto interactionChart(
    const Sequence& target,
    const Sequence& query,
    const Interaction& interaction,
    const bool bestSeedOnly) -> InteractionChart {
    InteractionChart chart;
    chart.pairs.reserve(interaction.pairs.size());
    for (Index index = 0U; index < interaction.pairs.size(); ++index) {
        chart.pairs.push_back({interaction.pairs[index], index});
    }
    std::ranges::sort(chart.pairs, {}, [](const DisplayPair& value) {
        return value.pair.target;
    });

    std::vector<std::size_t> columns(chart.pairs.size(), 0U);
    for (Index index = 1U; index < chart.pairs.size(); ++index) {
        const auto& previous = chart.pairs[index - 1U].pair;
        const auto& current = chart.pairs[index].pair;
        const auto targetGap = current.target - previous.target - 1U;
        const auto queryGap = previous.query - current.query - 1U;
        columns[index] = columns[index - 1U] + 1U + std::max(targetGap, queryGap);
    }

    const auto width = columns.back() + 1U;
    chart.targetBackbone.assign(width, ' ');
    chart.targetPairs.assign(width, ' ');
    chart.pairSymbols.assign(width, ' ');
    chart.queryPairs.assign(width, ' ');
    chart.queryBackbone.assign(width, ' ');

    const auto seeds = selectedSeeds(interaction, bestSeedOnly);
    for (Index index = 0U; index < chart.pairs.size(); ++index) {
        const auto& display = chart.pairs[index];
        const auto column = columns[index];
        chart.targetPairs[column] = target[display.pair.target];
        chart.queryPairs[column] = query[display.pair.query];
        const auto inSeed = std::ranges::any_of(seeds, [&](const SeedMatch* seed) {
            const auto first = std::min(seed->firstPair, seed->lastPair);
            const auto last = std::max(seed->firstPair, seed->lastPair);
            return display.originalIndex >= first && display.originalIndex <= last;
        });
        chart.pairSymbols[column] = inSeed ? '+' :
            (isGuPair(target[display.pair.target], query[display.pair.query]) ? ':' : '|');

        if (index + 1U == chart.pairs.size()) continue;
        const auto& next = chart.pairs[index + 1U].pair;
        const auto targetGap = next.target - display.pair.target - 1U;
        const auto queryGap = display.pair.query - next.query - 1U;
        const auto gapWidth = std::max(targetGap, queryGap);
        const auto gapStart = column + 1U;
        for (Index offset = 0U; offset < gapWidth; ++offset) {
            chart.targetBackbone[gapStart + offset] = offset < targetGap
                ? target[display.pair.target + 1U + offset] : '-';
            chart.queryBackbone[gapStart + offset] = offset < queryGap
                ? query[display.pair.query - 1U - offset] : '-';
        }
    }
    return chart;
}

[[nodiscard]] auto asciiInteraction(
    const Sequence& target,
    const Sequence& query,
    const Interaction& interaction,
    const bool detailed,
    const bool bestSeedOnly) -> std::string {
    constexpr std::size_t firstPairColumn{13U};
    const auto chart = interactionChart(target, query, interaction, bestSeedOnly);
    const auto& first = chart.pairs.front().pair;
    const auto& last = chart.pairs.back().pair;
    const auto lastPairColumn = firstPairColumn + chart.targetPairs.size() - 1U;

    const auto targetLeft = leftFlank(target.str().substr(0U, first.target));
    const auto targetRight = rightFlank(target.str().substr(last.target + 1U));
    const auto queryLeft = leftFlank(reverseString(query.str().substr(first.query + 1U)));
    const auto queryRight = rightFlank(reverseString(query.str().substr(0U, last.query)));
    const auto outerLine = [](const std::string& left, const std::string_view orientation,
                              const std::string& backbone, const std::string& right,
                              const std::string_view terminus) {
        return std::string(10U - left.size(), ' ') + std::string(orientation) + left +
               backbone + right + std::string(terminus);
    };

    std::ostringstream output;
    output.imbue(std::locale::classic());
    output << '\n' << target.id() << '\n'
           << coordinateLine(firstPairColumn, lastPairColumn,
                  annotatedIndex(target, first.target), annotatedIndex(target, last.target)) << '\n'
           << markerLine(firstPairColumn, lastPairColumn) << '\n'
           << outerLine(targetLeft, "5'-", chart.targetBackbone, targetRight, "-3'") << '\n'
           << std::string(firstPairColumn, ' ') << chart.targetPairs << '\n'
           << std::string(firstPairColumn, ' ') << chart.pairSymbols << '\n'
           << std::string(firstPairColumn, ' ') << chart.queryPairs << '\n'
           << outerLine(queryLeft, "3'-", chart.queryBackbone, queryRight, "-5'") << '\n'
           << markerLine(firstPairColumn, lastPairColumn) << '\n'
           << coordinateLine(firstPairColumn, lastPairColumn,
                  annotatedIndex(query, first.query), annotatedIndex(query, last.query)) << '\n'
           << query.id() << "\n\n";

    if (detailed) {
        output << "interaction seq1   = " << target.externalIndex(first.target) << ".."
               << target.externalIndex(last.target) << '\n'
               << "interaction seq2   = " << query.externalIndex(last.query) << ".."
               << query.externalIndex(first.query) << "\n\n";
    }

    output << "interaction energy = " << energyNumber(interaction.energy.total())
           << " kcal/mol\n";
    if (!detailed) return output.str();

    output << "  = E(init)        = " << energyNumber(interaction.energy.initiation) << '\n'
           << "  + E(loops)       = " << energyNumber(interaction.energy.loops) << '\n'
           << "  + E(dangleLeft)  = " << energyNumber(interaction.energy.dangleLeft) << '\n'
           << "  + E(dangleRight) = " << energyNumber(interaction.energy.dangleRight) << '\n'
           << "  + E(endLeft)     = " << energyNumber(interaction.energy.endLeft) << '\n'
           << "  + E(endRight)    = " << energyNumber(interaction.energy.endRight) << '\n'
           << "    : E(hybrid)    = " << energyNumber(interaction.energy.hybrid()) << '\n'
           << "  + ED(seq1)       = " << energyNumber(interaction.energy.openingTarget) << '\n'
           << "    : Pu(seq1)     = " << number(interaction.unpairedTarget) << '\n'
           << "  + ED(seq2)       = " << energyNumber(interaction.energy.openingQuery) << '\n'
           << "    : Pu(seq2)     = " << number(interaction.unpairedQuery) << '\n';
    if (centikcal(interaction.energy.additive) != 0.0) {
        output << "  + E(add)         = " << energyNumber(interaction.energy.additive) << '\n';
    }
    const auto seeds = selectedSeeds(interaction, bestSeedOnly);
    if (!seeds.empty()) {
        const auto seedCoordinate = [&](const SeedMatch& seed, const bool targetSide) {
            const auto seedFirst = std::min(seed.firstPair, seed.lastPair);
            const auto seedLast = std::max(seed.firstPair, seed.lastPair);
            const auto& seedLeft = interaction.pairs.at(seedFirst);
            const auto& seedRight = interaction.pairs.at(seedLast);
            return targetSide
                ? std::to_string(target.externalIndex(seedLeft.target)) + ".." +
                      std::to_string(target.externalIndex(seedRight.target))
                : std::to_string(query.externalIndex(seedRight.query)) + ".." +
                      std::to_string(query.externalIndex(seedLeft.query));
        };
        output << '\n'
               << "seed seq1   = " << joinedSeedValues(seeds, " | ", [&](const SeedMatch& seed) {
                      return seedCoordinate(seed, true);
                  }) << '\n'
               << "seed seq2   = " << joinedSeedValues(seeds, " | ", [&](const SeedMatch& seed) {
                      return seedCoordinate(seed, false);
                  }) << '\n'
               << "seed energy = " << joinedSeedValues(seeds, " | ", [](const SeedMatch& seed) {
                      return energyNumber(seed.energy);
                  }) << '\n'
               << "seed ED1    = " << joinedSeedValues(seeds, " | ", [](const SeedMatch& seed) {
                      return energyNumber(seed.openingTarget);
                  }) << '\n'
               << "seed ED2    = " << joinedSeedValues(seeds, " | ", [](const SeedMatch& seed) {
                      return energyNumber(seed.openingQuery);
                  }) << '\n'
               << "seed Pu1    = " << joinedSeedValues(seeds, " | ", [&](const SeedMatch& seed) {
                      return number(seed.unpairedTarget);
                  }) << '\n'
               << "seed Pu2    = " << joinedSeedValues(seeds, " | ", [&](const SeedMatch& seed) {
                      return number(seed.unpairedQuery);
                  }) << '\n';
    }
    return output.str();
}

[[nodiscard]] auto profile(
    const Sequence& sequence,
    const std::vector<Interaction>& interactions,
    const bool querySide,
    const bool probability,
    const char separator) -> std::string {
    std::ostringstream output;
    output << "idx" << separator << sequence.id() << separator
           << (probability ? "spotProb" : "minE") << '\n';
    for (Index position = 0U; position < sequence.size(); ++position) {
        double value = probability ? 0.0 : infinity;
        for (const auto& interaction : interactions) {
            const auto interval = querySide ? interaction.queryRange() : interaction.targetRange();
            if (!interval.contains(position)) continue;
            if (probability) value += interaction.probability;
            else value = std::min(value, interaction.energy.total());
        }
        output << sequence.externalIndex(position) << separator << sequence[position]
               << separator << number(value) << '\n';
    }
    return output.str();
}

[[nodiscard]] auto effectiveAccessibilityLength(
    const Sequence& sequence,
    const SideConfig& config) noexcept -> Index {
    auto length = sequence.size();
    if (config.interactionLengthMax != 0U) {
        length = std::min(length, config.interactionLengthMax);
    }
    if (config.accessibility == AccessibilityKind::compute &&
        config.accessibilityWindow != 0U) {
        length = std::min(length, config.accessibilityWindow);
    }
    return length;
}

[[nodiscard]] auto accessibilityTable(
    const Sequence& sequence,
    const SideConfig& config,
    const AccessibilityProvider& accessibility,
    const bool probabilities) -> std::string {
    std::ostringstream output;
    output.imbue(std::locale::classic());
    output << (probabilities ? "#unpaired probabilities\n"
                             : "#ensemble delta energy to unpair a region ED\n");
    const auto maxLength = effectiveAccessibilityLength(sequence, config);
    output << " #i$\tl=1";
    for (Index length = 2U; length <= maxLength; ++length) output << '\t' << length;
    output << "\t\n";
    for (Index end = 0U; end < sequence.size(); ++end) {
        output << sequence.externalIndex(end) << '\t';
        for (Index length = 1U; length <= maxLength; ++length) {
            if (length > end + 1U) {
                output << "NA\t";
                continue;
            }
            const Interval interval{end + 1U - length, end};
            const auto value = probabilities
                ? accessibility.unpairedProbability(interval)
                : accessibility.openingEnergy(interval);
            output << scientificNumber(value) << '\t';
        }
        output << '\n';
    }
    return output.str();
}

[[nodiscard]] auto coveringProbability(
    const std::vector<Interaction>& sites,
    const Index targetPosition,
    const Index queryPosition) noexcept -> double {
    double probability{};
    for (const auto& interaction : sites) {
        if (interaction.targetRange().contains(targetPosition) &&
            interaction.queryRange().contains(queryPosition)) {
            probability += interaction.probability;
        }
    }
    return std::clamp(probability, 0.0, 1.0);
}

} // namespace

auto OutputFormatter::primary(
    const Config& config,
    const Sequence& target,
    const Sequence& query,
    const PredictionResult& result,
    const bool includeCsvHeader) -> std::expected<std::string, std::string> {
    if (config.output.mode == OutputMode::csv) {
        return csv(config, target, query, result, includeCsvHeader);
    }
    if (config.output.mode == OutputMode::ensemble) {
        const auto interactionEnsemble = result.ensembleSites.empty()
            ? 0.0 : result.ensembleFreeEnergy;
        // EallTotal describes formation of an interaction plus both monomer
        // ensembles. With no favorable interaction there is no formation
        // event, so the public contract reports zero rather than the sum of
        // unrelated monomer folding free energies.
        const auto totalEnsemble = result.ensembleSites.empty()
            ? 0.0 : totalEnsembleEnergy(result);
        std::ostringstream output;
        output << "id1 " << target.id() << '\n'
               << "id2 " << query.id() << '\n'
               << "RT " << number(result.rt) << '\n'
               << "Eall " << fixedEnergyNumber(interactionEnsemble) << '\n'
               << "Eall1 " << fixedEnergyNumber(result.targetEnsembleFreeEnergy) << '\n'
               << "Eall2 " << fixedEnergyNumber(result.queryEnsembleFreeEnergy) << '\n'
               << "EallTotal " << fixedEnergyNumber(totalEnsemble) << '\n';
        return output.str();
    }
    if (result.interactions.empty()) {
        return std::string("\nno favorable interaction found\n");
    }
    std::string output;
    for (const auto& interaction : result.interactions) {
        output += asciiInteraction(target, query, interaction,
            config.output.mode == OutputMode::detailed,
            config.output.bestSeedOnly);
    }
    return output;
}

auto OutputFormatter::auxiliary(
    const std::string_view descriptor,
    const Config& config,
    const Sequence& target,
    const Sequence& query,
    const AccessibilityProvider& targetAccessibility,
    const AccessibilityProvider& queryAccessibility,
    const PredictionResult& result) -> std::expected<std::string, std::string> {
    const auto separator = descriptor.find(':');
    const auto originalKind = descriptor.substr(0U, separator);
    const auto kind = lowerAscii(originalKind);
    if (kind == "qmine") return profile(query, result.ensembleSites, true, false, config.output.separator);
    if (kind == "tmine") return profile(target, result.ensembleSites, false, false, config.output.separator);
    if (kind == "qspotprob") return profile(query, result.ensembleSites, true, true, config.output.separator);
    if (kind == "tspotprob") return profile(target, result.ensembleSites, false, true, config.output.separator);
    if (kind == "qacc") return accessibilityTable(query, config.query, queryAccessibility, false);
    if (kind == "qpu") return accessibilityTable(query, config.query, queryAccessibility, true);
    if (kind == "tacc") return accessibilityTable(target, config.target, targetAccessibility, false);
    if (kind == "tpu") return accessibilityTable(target, config.target, targetAccessibility, true);
    if (kind == "pmine") {
        std::ostringstream output;
        output << "minE";
        for (Index queryPosition = 0U; queryPosition < query.size(); ++queryPosition) {
            output << config.output.separator << query[queryPosition] << '_'
                   << query.externalIndex(queryPosition);
        }
        output << '\n';
        for (Index targetPosition = 0U; targetPosition < target.size(); ++targetPosition) {
            output << target[targetPosition] << '_' << target.externalIndex(targetPosition);
            for (Index queryPosition = 0U; queryPosition < query.size(); ++queryPosition) {
                auto minimum = infinity;
                for (const auto& interaction : result.ensembleSites) {
                    if (interaction.targetRange().contains(targetPosition) &&
                        interaction.queryRange().contains(queryPosition)) {
                        minimum = std::min(minimum, interaction.energy.total());
                    }
                }
                output << config.output.separator << number(minimum);
            }
            output << '\n';
        }
        return output.str();
    }
    if (kind == "spotprob") {
        std::string spots;
        if (separator != std::string_view::npos) {
            const auto second = descriptor.find(':', separator + 1U);
            if (second != std::string_view::npos) {
                spots = std::string(descriptor.substr(
                    separator + 1U, second - separator - 1U));
            }
        }
        if (spots.empty()) {
            std::ostringstream output;
            output << "spotProb";
            for (Index queryPosition = 0U; queryPosition < query.size(); ++queryPosition) {
                output << config.output.separator << query[queryPosition] << '_'
                       << query.externalIndex(queryPosition);
            }
            output << '\n';
            for (Index targetPosition = 0U; targetPosition < target.size(); ++targetPosition) {
                output << target[targetPosition] << '_' << target.externalIndex(targetPosition);
                for (Index queryPosition = 0U; queryPosition < query.size(); ++queryPosition) {
                    output << config.output.separator << number(coveringProbability(
                        result.ensembleSites, targetPosition, queryPosition));
                }
                output << '\n';
            }
            return output.str();
        }

        struct RequestedSpot {
            std::string label;
            Index target{};
            Index query{};
        };
        std::vector<RequestedSpot> requested;
        for (const auto& token : split(spots, ',')) {
            if (token.empty()) continue;
            const auto ampersand = token.find('&');
            if (ampersand == std::string::npos) return std::unexpected("spot must use target&query coordinates");
            long long targetExternal{};
            long long queryExternal{};
            const auto targetResult = std::from_chars(token.data(), token.data() + ampersand, targetExternal);
            const auto queryResult = std::from_chars(token.data() + ampersand + 1U,
                token.data() + token.size(), queryExternal);
            if (targetResult.ec != std::errc{} ||
                targetResult.ptr != token.data() + ampersand ||
                queryResult.ec != std::errc{} ||
                queryResult.ptr != token.data() + token.size()) {
                return std::unexpected("invalid spot coordinate");
            }
            auto targetPosition = target.internalIndex(targetExternal);
            auto queryPosition = query.internalIndex(queryExternal);
            if (!targetPosition || !queryPosition) return std::unexpected("spot coordinate is out of range");
            requested.push_back({token, *targetPosition, *queryPosition});
        }
        double covered{};
        for (const auto& interaction : result.ensembleSites) {
            if (std::ranges::any_of(requested, [&](const RequestedSpot& spot) {
                    return interaction.targetRange().contains(spot.target) &&
                           interaction.queryRange().contains(spot.query);
                })) {
                covered += interaction.probability;
            }
        }
        const auto targetBefore = target.externalIndex(0U);
        const auto queryBefore = query.externalIndex(0U);
        if (targetBefore == std::numeric_limits<long long>::min() ||
            queryBefore == std::numeric_limits<long long>::min()) {
            return std::unexpected("cannot encode the coordinate preceding the sequence origin");
        }
        std::ostringstream output;
        output << "spot" << config.output.separator << "probability\n"
               << targetBefore - 1LL << '&' << queryBefore - 1LL << config.output.separator
               << number(std::clamp(1.0 - covered, 0.0, 1.0)) << '\n';
        for (const auto& spot : requested) {
            output << spot.label << config.output.separator
                   << number(coveringProbability(result.ensembleSites, spot.target, spot.query))
                   << '\n';
        }
        return output.str();
    }
    return std::unexpected("unknown auxiliary output kind '" + std::string(originalKind) + "'");
}

auto writeOutput(const std::string_view destination, const std::string_view content)
    -> std::expected<void, std::string> {
    const std::array artifact{OutputArtifact{std::string(destination), std::string(content)}};
    return publishOutputs(artifact);
}

} // namespace intarnanew
