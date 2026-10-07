#include "intarnanew/tools/csv.hpp"
#include "intarnanew/tools/svg.hpp"

#include <charconv>
#include <fstream>
#include <iostream>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace {

[[nodiscard]] auto help() -> std::string {
    return R"(intarnanew-svg - direct deterministic SVG profiles and heatmaps

Usage:
  intarnanew-svg profile INPUT OUTPUT [options]
  intarnanew-svg heatmap INPUT OUTPUT [options]
  intarnanew-svg regions INPUT OUTPUT [options]

Options:
  --x-column NAME       profile x column (default: first column)
  --y-column NAME       profile value column (default: last column)
  --title TEXT          plot title
  --x-label TEXT        horizontal axis label
  --y-label TEXT        vertical axis label
  --width PIXELS        SVG width (default: 960)
  --height PIXELS       SVG height (profile 480, heatmap 720)
  --separator CHAR      force CSV separator; 'tab' selects TSV
  --keep-positive       heatmap: do not clamp positive energies to zero
  --sequence SIDE       regions: use id/start/end columns suffixed 1 or 2
  -h, --help            show this help

Profile input is a headered table; NA and empty values split line segments.
Heatmap input follows IntaRNA pMinE layout: the first column supplies row labels
and all remaining column names are x labels. NA and empty cells remain missing.
Regions input follows IntaRNA CSV output and defaults to id1/start1/end1.
)";
}

struct Arguments {
    std::string mode;
    std::string input;
    std::string output;
    std::optional<std::string> xColumn;
    std::optional<std::string> yColumn;
    std::optional<char> separator;
    std::string title;
    std::string xLabel;
    std::string yLabel;
    std::size_t width{960U};
    std::size_t height{};
    bool keepPositive{};
    unsigned sequenceSide{1U};
};

[[nodiscard]] auto parseSize(const std::string_view text) -> std::expected<std::size_t, std::string> {
    std::size_t value{};
    const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
    if (error != std::errc{} || end != text.data() + text.size()) {
        return std::unexpected("invalid pixel size '" + std::string{text} + "'");
    }
    return value;
}

[[nodiscard]] auto parseCoordinate(const std::string_view text)
    -> std::expected<long long, std::string> {
    long long value{};
    const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
    if (error != std::errc{} || end != text.data() + text.size()) {
        return std::unexpected("invalid integral coordinate '" + std::string{text} + "'");
    }
    return value;
}

[[nodiscard]] auto parse(const int argc, const char* const* argv)
    -> std::expected<std::optional<Arguments>, std::string> {
    if (argc == 2 &&
        (std::string_view{argv[1]} == "--help" || std::string_view{argv[1]} == "-h")) {
        return std::optional<Arguments>{};
    }
    if (argc < 4) return std::unexpected("expected MODE INPUT OUTPUT; use --help");
    Arguments result;
    result.mode = argv[1];
    result.input = argv[2];
    result.output = argv[3];
    if (result.mode != "profile" && result.mode != "heatmap" && result.mode != "regions") {
        return std::unexpected("mode must be 'profile', 'heatmap', or 'regions'");
    }
    result.height = result.mode == "profile" ? 480U : 720U;
    for (int index = 4; index < argc; ++index) {
        const std::string_view option{argv[index]};
        const auto value = [&]() -> std::expected<std::string_view, std::string> {
            if (++index >= argc) return std::unexpected("missing value for " + std::string{option});
            return argv[index];
        };
        if (option == "-h" || option == "--help") return std::optional<Arguments>{};
        if (option == "--x-column") {
            auto item = value(); if (!item) return std::unexpected(item.error()); result.xColumn = *item;
        } else if (option == "--y-column") {
            auto item = value(); if (!item) return std::unexpected(item.error()); result.yColumn = *item;
        } else if (option == "--title") {
            auto item = value(); if (!item) return std::unexpected(item.error()); result.title = *item;
        } else if (option == "--x-label") {
            auto item = value(); if (!item) return std::unexpected(item.error()); result.xLabel = *item;
        } else if (option == "--y-label") {
            auto item = value(); if (!item) return std::unexpected(item.error()); result.yLabel = *item;
        } else if (option == "--width" || option == "--height") {
            auto item = value(); if (!item) return std::unexpected(item.error());
            auto size = parseSize(*item); if (!size) return std::unexpected(size.error());
            if (option == "--width") result.width = *size; else result.height = *size;
        } else if (option == "--separator") {
            auto item = value(); if (!item) return std::unexpected(item.error());
            if (*item == "tab") result.separator = '\t';
            else if (item->size() == 1U) result.separator = item->front();
            else return std::unexpected("separator must be one character or 'tab'");
        } else if (option == "--sequence") {
            auto item = value(); if (!item) return std::unexpected(item.error());
            if (*item == "1") result.sequenceSide = 1U;
            else if (*item == "2") result.sequenceSide = 2U;
            else return std::unexpected("--sequence must be 1 or 2");
        } else if (option == "--keep-positive") result.keepPositive = true;
        else return std::unexpected("unknown option '" + std::string{option} + "'");
    }
    return std::optional<Arguments>{std::move(result)};
}

[[nodiscard]] auto missing(const std::string_view text) noexcept -> bool {
    return text.empty() || text == "NA" || text == "NaN" || text == "nan";
}

[[nodiscard]] auto makeSvg(
    const intarnanew::tools::CsvTable& table,
    const Arguments& arguments) -> std::expected<std::string, std::string> {
    if (arguments.mode == "profile") {
        const auto x = arguments.xColumn
            ? table.column(*arguments.xColumn)
            : std::optional<std::size_t>{0U};
        const auto y = arguments.yColumn ? table.column(*arguments.yColumn)
                                         : std::optional{table.header.size() - 1U};
        if (!x) return std::unexpected("profile x column not found");
        if (!y) return std::unexpected("profile y column not found");
        std::vector<intarnanew::tools::ProfilePoint> points;
        for (std::size_t row = 0U; row < table.rows.size(); ++row) {
            auto position = intarnanew::tools::parseFiniteNumber(table.rows[row][*x]);
            if (!position) return std::unexpected("row " + std::to_string(row + 2U) + ": " + position.error());
            std::optional<double> value;
            if (!missing(table.rows[row][*y])) {
                auto parsed = intarnanew::tools::parseFiniteNumber(table.rows[row][*y]);
                if (!parsed) return std::unexpected("row " + std::to_string(row + 2U) + ": " + parsed.error());
                value = *parsed;
            }
            points.push_back(intarnanew::tools::ProfilePoint{*position, value});
        }
        intarnanew::tools::ProfileSvgOptions options;
        options.width = arguments.width; options.height = arguments.height;
        if (!arguments.title.empty()) options.title = arguments.title;
        options.xLabel = arguments.xLabel.empty() ? table.header[*x] : arguments.xLabel;
        options.yLabel = arguments.yLabel.empty() ? table.header[*y] : arguments.yLabel;
        return intarnanew::tools::profileSvg(points, options);
    }
    if (arguments.mode == "regions") {
        const auto suffix = std::to_string(arguments.sequenceSide);
        const auto idColumn = table.column("id" + suffix);
        const auto startColumn = table.column("start" + suffix);
        const auto endColumn = table.column("end" + suffix);
        if (!idColumn || !startColumn || !endColumn) {
            return std::unexpected(
                "regions input requires id" + suffix + ", start" + suffix +
                ", and end" + suffix + " columns");
        }
        std::vector<intarnanew::tools::RegionSpan> regions;
        regions.reserve(table.rows.size());
        for (std::size_t row = 0U; row < table.rows.size(); ++row) {
            auto start = parseCoordinate(table.rows[row][*startColumn]);
            auto end = parseCoordinate(table.rows[row][*endColumn]);
            if (!start || !end) {
                return std::unexpected(
                    "row " + std::to_string(row + 2U) +
                    " has a non-integral or out-of-range region coordinate");
            }
            regions.push_back({
                table.rows[row][*idColumn], *start, *end});
        }
        intarnanew::tools::RegionSvgOptions options;
        options.width = arguments.width; options.height = arguments.height;
        if (!arguments.title.empty()) options.title = arguments.title;
        if (!arguments.xLabel.empty()) options.xLabel = arguments.xLabel;
        return intarnanew::tools::regionsSvg(regions, options);
    }
    if (table.header.size() < 2U) return std::unexpected("heatmap needs a row-label column and at least one value column");
    intarnanew::tools::HeatmapData data;
    data.xLabels.assign(table.header.begin() + 1, table.header.end());
    for (std::size_t row = 0U; row < table.rows.size(); ++row) {
        data.yLabels.push_back(table.rows[row][0]);
        for (std::size_t column = 1U; column < table.header.size(); ++column) {
            if (missing(table.rows[row][column])) data.values.emplace_back(std::nullopt);
            else {
                auto parsed = intarnanew::tools::parseFiniteNumber(table.rows[row][column]);
                if (!parsed) return std::unexpected("row " + std::to_string(row + 2U) + ": " + parsed.error());
                data.values.emplace_back(*parsed);
            }
        }
    }
    intarnanew::tools::HeatmapSvgOptions options;
    options.width = arguments.width; options.height = arguments.height;
    options.clampPositiveToZero = !arguments.keepPositive;
    if (!arguments.title.empty()) options.title = arguments.title;
    if (!arguments.xLabel.empty()) options.xLabel = arguments.xLabel;
    if (!arguments.yLabel.empty()) options.yLabel = arguments.yLabel;
    return intarnanew::tools::heatmapSvg(data, options);
}

} // namespace

int main(const int argc, const char* const* argv) {
    auto arguments = parse(argc, argv);
    if (!arguments) { std::cerr << "intarnanew-svg: " << arguments.error() << '\n'; return 2; }
    if (!*arguments) { std::cout << help(); return 0; }
    std::ifstream input((**arguments).input, std::ios::binary);
    if (!input) { std::cerr << "intarnanew-svg: cannot open input file\n"; return 1; }
    auto table = intarnanew::tools::readCsv(input, {.separator = (**arguments).separator});
    if (!table) { std::cerr << "intarnanew-svg: " << table.error() << '\n'; return 1; }
    auto svg = makeSvg(*table, **arguments);
    if (!svg) { std::cerr << "intarnanew-svg: " << svg.error() << '\n'; return 1; }
    std::ofstream output((**arguments).output, std::ios::binary | std::ios::trunc);
    if (!output) { std::cerr << "intarnanew-svg: cannot open output file\n"; return 1; }
    output << *svg;
    if (!output) { std::cerr << "intarnanew-svg: failed while writing output\n"; return 1; }
    return 0;
}
