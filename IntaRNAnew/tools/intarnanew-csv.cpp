#include "intarnanew/tools/csv.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace {

[[nodiscard]] auto help() -> std::string {
    return R"(intarnanew-csv - lossless deterministic interaction CSV fusion

Usage: intarnanew-csv [options] INPUT [INPUT ...]

Options:
  -o, --output FILE       output file, or '-' for standard output (default: -)
      --separator CHAR    force input and output separator; 'tab' selects TSV
      --source-column N   append source file name in column N
      --deduplicate       remove exact duplicate rows, preserving the first
  -h, --help              show this help

Quoted fields, embedded separators/newlines, and escaped quotes are supported.
The output schema is the stable union of input schemas; absent fields are empty.
)";
}

struct Arguments {
    std::vector<std::string> inputs;
    std::string output{"-"};
    std::optional<char> separator;
    std::optional<std::string> sourceColumn;
    bool deduplicate{};
};

[[nodiscard]] auto parse(const int argc, const char* const* argv)
    -> std::expected<std::optional<Arguments>, std::string> {
    Arguments result;
    for (int index = 1; index < argc; ++index) {
        const std::string_view option{argv[index]};
        const auto value = [&]() -> std::expected<std::string_view, std::string> {
            if (++index >= argc) return std::unexpected("missing value for " + std::string{option});
            return argv[index];
        };
        if (option == "-h" || option == "--help") return std::optional<Arguments>{};
        if (option == "-o" || option == "--output") {
            auto item = value(); if (!item) return std::unexpected(item.error());
            result.output = *item;
        } else if (option == "--separator") {
            auto item = value(); if (!item) return std::unexpected(item.error());
            if (*item == "tab") result.separator = '\t';
            else if (item->size() == 1U) result.separator = item->front();
            else return std::unexpected("separator must be one character or 'tab'");
        } else if (option == "--source-column") {
            auto item = value(); if (!item) return std::unexpected(item.error());
            result.sourceColumn = *item;
        } else if (option == "--deduplicate") {
            result.deduplicate = true;
        } else if (option.starts_with('-') && option != "-") {
            return std::unexpected("unknown option '" + std::string{option} + "'");
        } else {
            result.inputs.emplace_back(option);
        }
    }
    if (result.inputs.empty()) return std::unexpected("at least one input is required");
    const auto stdinCount = std::ranges::count(result.inputs, "-");
    if (stdinCount > 1) return std::unexpected("standard input can appear only once");
    return std::optional<Arguments>{std::move(result)};
}

[[nodiscard]] auto writeOutput(const std::string& destination, const std::string& content)
    -> std::expected<void, std::string> {
    if (destination == "-") {
        std::cout << content;
        if (!std::cout) return std::unexpected("failed to write standard output");
        return {};
    }
    std::ofstream output(destination, std::ios::binary | std::ios::trunc);
    if (!output) return std::unexpected("cannot open output file '" + destination + "'");
    output << content;
    if (!output) return std::unexpected("failed while writing output file '" + destination + "'");
    return {};
}

[[nodiscard]] auto readInput(
    const std::string& inputName,
    const std::optional<char> separator) -> std::expected<intarnanew::tools::CsvTable, std::string> {
    if (inputName == "-") {
        return intarnanew::tools::readCsv(std::cin, {.separator = separator});
    }
    std::ifstream input(inputName, std::ios::binary);
    if (!input) {
        return std::unexpected("cannot open '" + inputName + "'");
    }
    return intarnanew::tools::readCsv(input, {.separator = separator});
}

} // namespace

int main(const int argc, const char* const* argv) {
    auto arguments = parse(argc, argv);
    if (!arguments) { std::cerr << "intarnanew-csv: " << arguments.error() << '\n'; return 2; }
    if (!*arguments) { std::cout << help(); return 0; }
    std::vector<intarnanew::tools::CsvTable> tables;
    std::vector<std::string> labels;
    for (const auto& inputName : (**arguments).inputs) {
        auto table = readInput(inputName, (**arguments).separator);
        if (inputName == "-") {
            labels.emplace_back("STDIN");
        } else {
            labels.push_back(std::filesystem::path(inputName).filename().string());
        }
        if (!table) { std::cerr << "intarnanew-csv: " << inputName << ": " << table.error() << '\n'; return 1; }
        tables.push_back(std::move(*table));
    }
    auto fused = intarnanew::tools::fuseCsv(
        tables,
        labels,
        {.sourceColumn = (**arguments).sourceColumn, .deduplicate = (**arguments).deduplicate});
    if (!fused) { std::cerr << "intarnanew-csv: " << fused.error() << '\n'; return 1; }
    if ((**arguments).separator) fused->separator = *(**arguments).separator;
    auto output = intarnanew::tools::csvText(*fused);
    if (!output) { std::cerr << "intarnanew-csv: " << output.error() << '\n'; return 1; }
    auto written = writeOutput((**arguments).output, *output);
    if (!written) { std::cerr << "intarnanew-csv: " << written.error() << '\n'; return 1; }
    return 0;
}
