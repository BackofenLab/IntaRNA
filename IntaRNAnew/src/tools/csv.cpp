#include "intarnanew/tools/csv.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <limits>
#include <ostream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

namespace intarnanew::tools {
namespace {

[[nodiscard]] auto detectSeparator(const std::string_view input) -> char {
    std::size_t semicolons{};
    std::size_t tabs{};
    std::size_t commas{};
    bool quoted{};
    for (std::size_t index = 0; index < input.size(); ++index) {
        const char symbol = input[index];
        if (symbol == '"') {
            if (quoted && index + 1U < input.size() && input[index + 1U] == '"') {
                ++index;
            } else {
                quoted = !quoted;
            }
        } else if (!quoted && (symbol == '\n' || symbol == '\r')) {
            break;
        } else if (!quoted) {
            semicolons += symbol == ';';
            tabs += symbol == '\t';
            commas += symbol == ',';
        }
    }
    if (tabs > semicolons && tabs >= commas) {
        return '\t';
    }
    if (commas > semicolons) {
        return ',';
    }
    return ';';
}

[[nodiscard]] auto parseRecords(
    const std::string_view input,
    const char separator,
    const bool allowEmptyLines)
    -> std::expected<std::vector<std::vector<std::string>>, std::string> {
    std::vector<std::vector<std::string>> records;
    std::vector<std::string> record;
    std::string field;
    bool quoted{};
    bool afterQuote{};

    const auto finishField = [&] {
        record.push_back(std::move(field));
        field.clear();
        afterQuote = false;
    };
    const auto finishRecord = [&] {
        const bool syntacticallyBlank = record.empty() && field.empty() && !afterQuote;
        finishField();
        if (!allowEmptyLines || !syntacticallyBlank) {
            records.push_back(std::move(record));
        }
        record.clear();
    };

    for (std::size_t index = 0; index < input.size(); ++index) {
        const char symbol = input[index];
        if (quoted) {
            if (symbol == '"') {
                if (index + 1U < input.size() && input[index + 1U] == '"') {
                    field.push_back('"');
                    ++index;
                } else {
                    quoted = false;
                    afterQuote = true;
                }
            } else {
                field.push_back(symbol);
            }
            continue;
        }

        if (afterQuote) {
            if (symbol == separator) {
                finishField();
            } else if (symbol == '\n') {
                finishRecord();
            } else if (symbol == '\r') {
                if (index + 1U < input.size() && input[index + 1U] == '\n') {
                    ++index;
                }
                finishRecord();
            } else if (symbol != ' ' && symbol != '\t') {
                return std::unexpected(
                    "unexpected character after a closing CSV quote at byte " +
                    std::to_string(index));
            }
            continue;
        }

        if (symbol == '"') {
            if (!field.empty()) {
                return std::unexpected(
                    "CSV quote starts inside an unquoted field at byte " +
                    std::to_string(index));
            }
            quoted = true;
        } else if (symbol == separator) {
            finishField();
        } else if (symbol == '\n') {
            finishRecord();
        } else if (symbol == '\r') {
            if (index + 1U < input.size() && input[index + 1U] == '\n') {
                ++index;
            }
            finishRecord();
        } else {
            field.push_back(symbol);
        }
    }

    if (quoted) {
        return std::unexpected("unterminated quoted CSV field");
    }
    if (!field.empty() || !record.empty() || afterQuote) {
        finishRecord();
    }
    return records;
}

[[nodiscard]] auto escaped(
    const std::string_view value,
    const char separator,
    const bool forceQuote = false) -> std::string {
    const bool quote = forceQuote || value.find(separator) != std::string_view::npos ||
                       value.find('"') != std::string_view::npos ||
                       value.find('\n') != std::string_view::npos ||
                       value.find('\r') != std::string_view::npos;
    if (!quote) {
        return std::string{value};
    }
    std::string result;
    result.reserve(value.size() + 2U);
    result.push_back('"');
    for (const char symbol : value) {
        if (symbol == '"') {
            result.push_back('"');
        }
        result.push_back(symbol);
    }
    result.push_back('"');
    return result;
}

[[nodiscard]] auto rowIdentity(const std::vector<std::string>& row) -> std::string {
    std::string identity;
    for (const auto& value : row) {
        identity += std::to_string(value.size());
        identity.push_back(':');
        identity += value;
        identity.push_back('|');
    }
    return identity;
}

} // namespace

auto CsvTable::column(const std::string_view name) const noexcept -> std::optional<std::size_t> {
    const auto found = std::find(header.begin(), header.end(), name);
    if (found == header.end()) {
        return std::nullopt;
    }
    return static_cast<std::size_t>(std::distance(header.begin(), found));
}

auto readCsv(std::istream& input, const CsvReadOptions options)
    -> std::expected<CsvTable, std::string> {
    std::ostringstream buffer;
    buffer << input.rdbuf();
    if (input.bad()) {
        return std::unexpected("failed while reading CSV input");
    }
    auto text = buffer.str();
    if (text.starts_with("\xEF\xBB\xBF")) {
        text.erase(0U, 3U);
    }
    const auto separator = options.separator.value_or(detectSeparator(text));
    if (separator == '"' || separator == '\r' || separator == '\n' || separator == '\0') {
        return std::unexpected("invalid CSV separator");
    }
    auto parsed = parseRecords(text, separator, options.allowEmptyLines);
    if (!parsed) {
        return std::unexpected(parsed.error());
    }
    if (parsed->empty()) {
        return std::unexpected("CSV input is empty");
    }

    CsvTable table;
    table.separator = separator;
    if (options.requireHeader) {
        table.header = std::move(parsed->front());
        parsed->erase(parsed->begin());
    } else {
        const auto width = parsed->front().size();
        table.header.reserve(width);
        for (std::size_t index = 0; index < width; ++index) {
            table.header.push_back("column" + std::to_string(index + 1U));
        }
    }
    if (table.header.empty()) {
        return std::unexpected("CSV header is empty");
    }
    std::unordered_set<std::string> headerNames;
    for (const auto& name : table.header) {
        if (name.empty()) {
            return std::unexpected("CSV contains an empty column name");
        }
        if (!headerNames.insert(name).second) {
            return std::unexpected("duplicate CSV column name '" + name + "'");
        }
    }
    for (std::size_t index = 0; index < parsed->size(); ++index) {
        if ((*parsed)[index].size() != table.header.size()) {
            return std::unexpected(
                "CSV record " + std::to_string(index + 2U) + " has " +
                std::to_string((*parsed)[index].size()) + " fields; expected " +
                std::to_string(table.header.size()));
        }
    }
    table.rows = std::move(*parsed);
    return table;
}

auto writeCsv(const CsvTable& table, std::ostream& output) -> std::expected<void, std::string> {
    if (table.header.empty()) {
        return std::unexpected("cannot write a CSV table without columns");
    }
    if (table.separator == '"' || table.separator == '\r' ||
        table.separator == '\n' || table.separator == '\0') {
        return std::unexpected("cannot write CSV with an invalid separator");
    }
    std::unordered_set<std::string> names;
    for (const auto& name : table.header) {
        if (name.empty()) {
            return std::unexpected("cannot write CSV with an empty column name");
        }
        if (!names.insert(name).second) {
            return std::unexpected("cannot write CSV with duplicate column '" + name + "'");
        }
    }
    const auto writeRow = [&](const std::vector<std::string>& row) {
        for (std::size_t index = 0; index < row.size(); ++index) {
            if (index != 0U) {
                output.put(table.separator);
            }
            // An unquoted empty value in a one-column table is indistinguishable
            // from an ignorable blank record when read back.
            output << escaped(
                row[index], table.separator, row.size() == 1U && row[index].empty());
        }
        output.put('\n');
    };
    writeRow(table.header);
    for (std::size_t index = 0; index < table.rows.size(); ++index) {
        if (table.rows[index].size() != table.header.size()) {
            return std::unexpected(
                "cannot write malformed CSV record " + std::to_string(index + 1U));
        }
        writeRow(table.rows[index]);
    }
    if (!output) {
        return std::unexpected("failed while writing CSV output");
    }
    return {};
}

auto csvText(const CsvTable& table) -> std::expected<std::string, std::string> {
    std::ostringstream output;
    auto written = writeCsv(table, output);
    if (!written) {
        return std::unexpected(written.error());
    }
    return output.str();
}

auto fuseCsv(
    const std::span<const CsvTable> tables,
    const std::span<const std::string> sourceLabels,
    const CsvFusionOptions options) -> std::expected<CsvTable, std::string> {
    if (tables.empty()) {
        return std::unexpected("at least one CSV table is required for fusion");
    }
    if (!sourceLabels.empty() && sourceLabels.size() != tables.size()) {
        return std::unexpected("source label count must equal table count");
    }

    CsvTable result;
    result.separator = tables.front().separator;
    std::unordered_map<std::string, std::size_t> outputIndices;
    for (std::size_t tableIndex = 0U; tableIndex < tables.size(); ++tableIndex) {
        const auto& table = tables[tableIndex];
        if (table.header.empty()) {
            return std::unexpected(
                "input table " + std::to_string(tableIndex + 1U) + " has no columns");
        }
        std::unordered_set<std::string> names;
        for (const auto& name : table.header) {
            if (name.empty()) {
                return std::unexpected(
                    "input table " + std::to_string(tableIndex + 1U) +
                    " contains an empty column name");
            }
            if (!names.insert(name).second) {
                return std::unexpected(
                    "input table " + std::to_string(tableIndex + 1U) +
                    " contains duplicate column '" + name + "'");
            }
            if (!outputIndices.contains(name)) {
                outputIndices.emplace(name, result.header.size());
                result.header.push_back(name);
            }
        }
        for (const auto& row : table.rows) {
            if (row.size() != table.header.size()) {
                return std::unexpected(
                    "input table " + std::to_string(tableIndex + 1U) + " is malformed");
            }
        }
    }
    if (options.sourceColumn) {
        if (options.sourceColumn->empty()) {
            return std::unexpected("source column name cannot be empty");
        }
        if (outputIndices.contains(*options.sourceColumn)) {
            return std::unexpected(
                "source column '" + *options.sourceColumn + "' already exists");
        }
        outputIndices.emplace(*options.sourceColumn, result.header.size());
        result.header.push_back(*options.sourceColumn);
    }

    std::unordered_set<std::string> seen;
    for (std::size_t tableIndex = 0; tableIndex < tables.size(); ++tableIndex) {
        const auto& table = tables[tableIndex];
        std::vector<std::size_t> projection;
        projection.reserve(table.header.size());
        for (const auto& name : table.header) {
            projection.push_back(outputIndices.at(name));
        }
        for (const auto& inputRow : table.rows) {
            if (inputRow.size() != table.header.size()) {
                return std::unexpected(
                    "input table " + std::to_string(tableIndex + 1U) + " is malformed");
            }
            std::vector<std::string> outputRow(result.header.size());
            for (std::size_t columnIndex = 0; columnIndex < inputRow.size(); ++columnIndex) {
                outputRow[projection[columnIndex]] = inputRow[columnIndex];
            }
            if (options.sourceColumn) {
                outputRow.back() = sourceLabels.empty()
                    ? std::to_string(tableIndex + 1U)
                    : sourceLabels[tableIndex];
            }
            if (!options.deduplicate || seen.insert(rowIdentity(outputRow)).second) {
                result.rows.push_back(std::move(outputRow));
            }
        }
    }
    return result;
}

auto parseFiniteNumber(const std::string_view text) -> std::expected<double, std::string> {
    if (text.empty()) {
        return std::unexpected("numeric value is empty");
    }
    double value{};
    const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
    if (error != std::errc{} || end != text.data() + text.size()) {
        return std::unexpected("invalid numeric value '" + std::string{text} + "'");
    }
    if (!std::isfinite(value)) {
        return std::unexpected("numeric value must be finite");
    }
    return value;
}

} // namespace intarnanew::tools
