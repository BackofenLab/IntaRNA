#pragma once

#include <cstddef>
#include <expected>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace intarnanew::tools {

struct ProfilePoint {
    double position{};
    std::optional<double> value;
};

struct ProfileSvgOptions {
    std::size_t width{960U};
    std::size_t height{480U};
    std::string title{"IntaRNA profile"};
    std::string xLabel{"sequence position"};
    std::string yLabel{"value"};
    std::string stroke{"#2166ac"};
    bool zeroLine{true};
};

// Generates a self-contained, deterministic SVG. Missing profile points split
// the line into independent segments instead of being silently interpolated.
[[nodiscard]] auto profileSvg(
    std::span<const ProfilePoint> points,
    const ProfileSvgOptions& options = {}) -> std::expected<std::string, std::string>;

struct HeatmapData {
    std::vector<std::string> xLabels;
    std::vector<std::string> yLabels;
    // Row-major: values[y * xLabels.size() + x]. Null values are missing.
    std::vector<std::optional<double>> values;
};

struct HeatmapSvgOptions {
    std::size_t width{960U};
    std::size_t height{720U};
    std::string title{"IntaRNA interaction-energy heatmap"};
    std::string xLabel{"query position"};
    std::string yLabel{"target position"};
    std::string missingColor{"#e6e6e6"};
    // Matches the documented pMinE visualization convention: positive
    // energies are treated as zero unless this is false.
    bool clampPositiveToZero{true};
};

[[nodiscard]] auto heatmapSvg(
    const HeatmapData& data,
    const HeatmapSvgOptions& options = {}) -> std::expected<std::string, std::string>;

struct RegionSpan {
    std::string id;
    long long start{};
    long long end{};
};

struct RegionSvgOptions {
    std::size_t width{960U};
    std::size_t height{720U};
    std::string title{"IntaRNA interaction-covered regions"};
    std::string xLabel{"sequence position"};
    std::string fill{"#4393c3"};
};

// Draws one labeled horizontal track per input span, preserving row order.
// Coordinates are inclusive and may be negative, but start must not exceed end.
[[nodiscard]] auto regionsSvg(
    std::span<const RegionSpan> regions,
    const RegionSvgOptions& options = {}) -> std::expected<std::string, std::string>;

[[nodiscard]] auto xmlEscape(std::string_view text) -> std::string;

} // namespace intarnanew::tools
