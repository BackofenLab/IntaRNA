#pragma once

#include <expected>
#include <iosfwd>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace intarnanew::tools {

struct CsvTable {
    std::vector<std::string> header;
    std::vector<std::vector<std::string>> rows;
    char separator{';'};

    [[nodiscard]] auto column(std::string_view name) const noexcept -> std::optional<std::size_t>;
};

struct CsvReadOptions {
    // A null separator auto-detects ';', tab, or ',' from the first record.
    std::optional<char> separator;
    bool requireHeader{true};
    bool allowEmptyLines{true};
};

[[nodiscard]] auto readCsv(std::istream& input, CsvReadOptions options = {})
    -> std::expected<CsvTable, std::string>;

[[nodiscard]] auto writeCsv(const CsvTable& table, std::ostream& output)
    -> std::expected<void, std::string>;

[[nodiscard]] auto csvText(const CsvTable& table) -> std::expected<std::string, std::string>;

struct CsvFusionOptions {
    // If set, this column is appended and populated from sourceLabels.
    std::optional<std::string> sourceColumn;
    // Exact duplicate output rows are removed while preserving first occurrence.
    bool deduplicate{};
};

// Produces the stable schema union of all inputs. Columns are ordered by first
// occurrence, and rows retain input-table and within-table order.
[[nodiscard]] auto fuseCsv(
    std::span<const CsvTable> tables,
    std::span<const std::string> sourceLabels = {},
    CsvFusionOptions options = {}) -> std::expected<CsvTable, std::string>;

[[nodiscard]] auto parseFiniteNumber(std::string_view text)
    -> std::expected<double, std::string>;

} // namespace intarnanew::tools
