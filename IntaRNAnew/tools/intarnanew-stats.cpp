#include "intarnanew/tools/csv.hpp"
#include "intarnanew/tools/statistics.hpp"

#include <array>
#include <charconv>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <span>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace {

using intarnanew::tools::AdjustmentMethod;
using intarnanew::tools::CsvTable;
using intarnanew::tools::DistributionFit;
using intarnanew::tools::DistributionKind;

struct Arguments {
    std::string command;
    std::string input{"-"};
    std::string output{"-"};
    std::string column{"E"};
    std::string outputColumn;
    std::optional<char> separator;
    DistributionKind distribution{DistributionKind::gev};
    AdjustmentMethod adjustment{AdjustmentMethod::benjaminiHochberg};
    bool empirical{};
    bool columnSpecified{};
    bool distributionSpecified{};
    bool adjustmentSpecified{};
    bool outputColumnSpecified{};
};

[[nodiscard]] auto help() -> std::string {
    return R"(intarnanew-stats - deterministic interaction-score statistics

Usage:
  intarnanew-stats fit [options]
  intarnanew-stats pvalue [options]
  intarnanew-stats adjust [options]

Commands:
  fit      fit a distribution and print name, location, scale, shape, and NLL
  pvalue   append lower-tail p-values for interaction energies
  adjust   append multiple-testing adjusted p-values

Options:
  -i, --input FILE          CSV/TSV input, or '-' for standard input
  -o, --output FILE         output file, or '-' for standard output
  -c, --column NAME         numeric input column (E; p-value for adjust)
      --output-column NAME  appended column (default: p-value or p-adjusted)
  -d, --distribution KIND  gev, gumbel, or gauss (default: gev)
      --empirical           pvalue: exact plus-one empirical lower tail
  -m, --method METHOD       adjust: none, bonferroni, holm, hochberg, bh, by
      --separator CHAR      force input/output separator; use 'tab' for TSV
  -h, --help                show this help

The fitted p-value is P(X <= E), because more-negative interaction energy is
more extreme. The empirical estimate is (#{sample <= E}+1)/(n+1). Missing,
infinite, and malformed numeric cells are rejected rather than imputed.
)";
}

[[nodiscard]] auto optionValue(
    const std::span<const std::string_view> values,
    std::size_t& index,
    const std::string_view option) -> std::expected<std::string_view, std::string> {
    if (index + 1U >= values.size()) {
        return std::unexpected("missing value for " + std::string{option});
    }
    return values[++index];
}

[[nodiscard]] auto parse(std::span<const std::string_view> values)
    -> std::expected<std::optional<Arguments>, std::string> {
    if (values.empty()) {
        return std::unexpected("missing command; use --help");
    }
    if (values.front() == "-h" || values.front() == "--help") {
        return std::optional<Arguments>{};
    }
    Arguments result;
    result.command = values.front();
    if (result.command != "fit" && result.command != "pvalue" && result.command != "adjust") {
        return std::unexpected("unknown command '" + result.command + "'");
    }
    for (std::size_t index = 1U; index < values.size(); ++index) {
        const auto option = values[index];
        if (option == "-h" || option == "--help") {
            return std::optional<Arguments>{};
        }
        auto get = [&] { return optionValue(values, index, option); };
        if (option == "-i" || option == "--input") {
            auto value = get(); if (!value) return std::unexpected(value.error());
            result.input = *value;
        } else if (option == "-o" || option == "--output") {
            auto value = get(); if (!value) return std::unexpected(value.error());
            result.output = *value;
        } else if (option == "-c" || option == "--column") {
            auto value = get(); if (!value) return std::unexpected(value.error());
            result.column = *value;
            result.columnSpecified = true;
        } else if (option == "--output-column") {
            auto value = get(); if (!value) return std::unexpected(value.error());
            result.outputColumn = *value;
            result.outputColumnSpecified = true;
        } else if (option == "-d" || option == "--distribution") {
            auto value = get(); if (!value) return std::unexpected(value.error());
            auto kind = intarnanew::tools::parseDistribution(*value);
            if (!kind) return std::unexpected(kind.error());
            result.distribution = *kind;
            result.distributionSpecified = true;
        } else if (option == "-m" || option == "--method") {
            auto value = get(); if (!value) return std::unexpected(value.error());
            auto method = intarnanew::tools::parseAdjustment(*value);
            if (!method) return std::unexpected(method.error());
            result.adjustment = *method;
            result.adjustmentSpecified = true;
        } else if (option == "--separator") {
            auto value = get(); if (!value) return std::unexpected(value.error());
            if (*value == "tab") result.separator = '\t';
            else if (value->size() == 1U) result.separator = value->front();
            else return std::unexpected("separator must be one character or 'tab'");
        } else if (option == "--empirical") {
            result.empirical = true;
        } else {
            return std::unexpected("unknown option '" + std::string{option} + "'");
        }
    }
    if (result.command == "adjust" && !result.columnSpecified) {
        result.column = "p-value";
    }
    if (result.command != "adjust" && result.adjustmentSpecified) {
        return std::unexpected("--method is only valid for the adjust command");
    }
    if (result.command == "adjust" && result.distributionSpecified) {
        return std::unexpected("--distribution is not valid for the adjust command");
    }
    if (result.command != "pvalue" && result.empirical) {
        return std::unexpected("--empirical is only valid for the pvalue command");
    }
    if (result.empirical && result.distributionSpecified) {
        return std::unexpected("--empirical and --distribution are mutually exclusive");
    }
    if (result.command == "fit" && result.outputColumnSpecified) {
        return std::unexpected("--output-column is not valid for the fit command");
    }
    if (result.command != "fit" && result.outputColumn.empty()) {
        result.outputColumn = result.command == "adjust" ? "p-adjusted" : "p-value";
    }
    return std::optional<Arguments>{std::move(result)};
}

[[nodiscard]] auto readTable(const Arguments& arguments) -> std::expected<CsvTable, std::string> {
    if (arguments.input == "-") {
        return intarnanew::tools::readCsv(
            std::cin, {.separator = arguments.separator});
    }
    std::ifstream input(arguments.input, std::ios::binary);
    if (!input) {
        return std::unexpected("cannot open input file '" + arguments.input + "'");
    }
    return intarnanew::tools::readCsv(input, {.separator = arguments.separator});
}

[[nodiscard]] auto extract(
    const CsvTable& table,
    const std::string_view column) -> std::expected<std::vector<double>, std::string> {
    const auto index = table.column(column);
    if (!index) {
        return std::unexpected("CSV has no column '" + std::string{column} + "'");
    }
    std::vector<double> values;
    values.reserve(table.rows.size());
    for (std::size_t row = 0U; row < table.rows.size(); ++row) {
        auto value = intarnanew::tools::parseFiniteNumber(table.rows[row][*index]);
        if (!value) {
            return std::unexpected(
                "record " + std::to_string(row + 2U) + ", column '" +
                std::string{column} + "': " + value.error());
        }
        values.push_back(*value);
    }
    return values;
}

[[nodiscard]] auto render(const double value) -> std::string {
    std::array<char, 64U> buffer{};
    const auto [end, error] = std::to_chars(
        buffer.data(), buffer.data() + buffer.size(), value);
    if (error != std::errc{}) {
        return {};
    }
    return std::string(buffer.data(), end);
}

[[nodiscard]] auto writeText(
    const std::string_view destination,
    const std::string_view text) -> std::expected<void, std::string> {
    if (destination == "-") {
        std::cout << text;
        if (!std::cout) return std::unexpected("failed to write standard output");
        return {};
    }
    const std::filesystem::path finalPath{destination};
    const auto temporary = finalPath.string() + ".tmp";
    {
        std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
        if (!output) return std::unexpected("cannot open temporary output file '" + temporary + "'");
        output << text;
        if (!output) return std::unexpected("failed while writing temporary output file");
    }
    std::error_code error;
    std::filesystem::rename(temporary, finalPath, error);
    if (error) {
        std::filesystem::remove(temporary);
        return std::unexpected("cannot install output file '" + finalPath.string() + "': " + error.message());
    }
    return {};
}

[[nodiscard]] auto execute(const Arguments& arguments) -> std::expected<void, std::string> {
    auto table = readTable(arguments);
    if (!table) return std::unexpected(table.error());
    auto values = extract(*table, arguments.column);
    if (!values) return std::unexpected(values.error());
    if (arguments.command == "fit") {
        auto fit = intarnanew::tools::fitDistribution(*values, arguments.distribution);
        if (!fit) return std::unexpected(fit.error());
        std::ostringstream output;
        output << "distribution\tlocation\tscale\tshape\tnll\tobservations\tconverged\titerations\n"
               << intarnanew::tools::distributionName(fit->kind) << '\t'
               << std::setprecision(17) << fit->location << '\t' << fit->scale << '\t'
               << fit->shape << '\t' << fit->negativeLogLikelihood << '\t'
               << fit->observations << '\t' << (fit->converged ? "true" : "false")
               << '\t' << fit->iterations << '\n';
        return writeText(arguments.output, output.str());
    }
    if (table->column(arguments.outputColumn)) {
        return std::unexpected("output column '" + arguments.outputColumn + "' already exists");
    }
    std::vector<double> probabilities;
    probabilities.reserve(values->size());
    if (arguments.command == "adjust") {
        auto adjusted = intarnanew::tools::adjustPValues(*values, arguments.adjustment);
        if (!adjusted) return std::unexpected(adjusted.error());
        probabilities = std::move(*adjusted);
    } else if (arguments.empirical) {
        for (const double value : *values) {
            auto probability = intarnanew::tools::empiricalInteractionPValue(value, *values);
            if (!probability) return std::unexpected(probability.error());
            probabilities.push_back(*probability);
        }
    } else {
        auto fit = intarnanew::tools::fitDistribution(*values, arguments.distribution);
        if (!fit) return std::unexpected(fit.error());
        for (const double value : *values) {
            auto probability = intarnanew::tools::interactionEnergyPValue(value, *fit);
            if (!probability) return std::unexpected(probability.error());
            probabilities.push_back(*probability);
        }
    }
    table->header.push_back(arguments.outputColumn);
    for (std::size_t row = 0U; row < table->rows.size(); ++row) {
        table->rows[row].push_back(render(probabilities[row]));
    }
    auto output = intarnanew::tools::csvText(*table);
    if (!output) return std::unexpected(output.error());
    return writeText(arguments.output, *output);
}

} // namespace

int main(const int argc, const char* const* argv) {
    std::vector<std::string_view> values;
    for (int index = 1; index < argc; ++index) values.emplace_back(argv[index]);
    auto arguments = parse(values);
    if (!arguments) {
        std::cerr << "intarnanew-stats: " << arguments.error() << '\n';
        return 2;
    }
    if (!*arguments) {
        std::cout << help();
        return 0;
    }
    auto result = execute(**arguments);
    if (!result) {
        std::cerr << "intarnanew-stats: " << result.error() << '\n';
        return 1;
    }
    return 0;
}
