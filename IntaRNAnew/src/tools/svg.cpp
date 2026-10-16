#include "intarnanew/tools/svg.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace intarnanew::tools {
namespace {

struct PlotBox {
    double left{};
    double top{};
    double width{};
    double height{};
};

[[nodiscard]] auto number(const double value) -> std::string {
    std::ostringstream output;
    output << std::fixed << std::setprecision(3) << value;
    auto text = output.str();
    while (text.size() > 1U && text.back() == '0') {
        text.pop_back();
    }
    if (!text.empty() && text.back() == '.') {
        text.pop_back();
    }
    if (text == "-0") {
        text = "0";
    }
    return text;
}

[[nodiscard]] auto validDimensions(
    const std::size_t width,
    const std::size_t height) -> std::expected<void, std::string> {
    if (width < 240U || height < 180U) {
        return std::unexpected("SVG dimensions must be at least 240 by 180 pixels");
    }
    if (width > 32768U || height > 32768U) {
        return std::unexpected("SVG dimensions exceed the 32768-pixel limit");
    }
    return {};
}

[[nodiscard]] auto documentStart(
    const std::size_t width,
    const std::size_t height,
    const std::string_view title) -> std::string {
    std::ostringstream output;
    output << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
           << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width
           << "\" height=\"" << height << "\" viewBox=\"0 0 " << width << ' '
           << height << "\" role=\"img\" aria-labelledby=\"title\">\n"
           << "<title id=\"title\">" << xmlEscape(title) << "</title>\n"
           << "<rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>\n"
           << "<style>text{font-family:system-ui,sans-serif;fill:#222}"
              ".axis{stroke:#222;stroke-width:1}.grid{stroke:#ddd;stroke-width:1}"
              ".tick{font-size:11px}.label{font-size:13px}.heading{font-size:17px;font-weight:600}"
              "</style>\n";
    return output.str();
}

[[nodiscard]] auto scale(
    const double value,
    const double minimum,
    const double maximum,
    const double outputMinimum,
    const double outputMaximum) noexcept -> double {
    if (minimum == maximum) {
        return 0.5 * (outputMinimum + outputMaximum);
    }
    return outputMinimum + (value - minimum) * (outputMaximum - outputMinimum) /
                               (maximum - minimum);
}

[[nodiscard]] auto tickValue(const double minimum, const double maximum, const int tick) noexcept
    -> double {
    return minimum + (maximum - minimum) * static_cast<double>(tick) / 5.0;
}

struct Rgb { int red{}; int green{}; int blue{}; };

[[nodiscard]] auto interpolate(const Rgb low, const Rgb high, const double fraction) noexcept -> Rgb {
    const auto channel = [fraction](const int start, const int end) {
        return static_cast<int>(std::lround(start + fraction * static_cast<double>(end - start)));
    };
    return {channel(low.red, high.red), channel(low.green, high.green), channel(low.blue, high.blue)};
}

[[nodiscard]] auto color(const double value, const double minimum, const double maximum) -> std::string {
    const double fraction = maximum == minimum
        ? 0.5
        : std::clamp((value - minimum) / (maximum - minimum), 0.0, 1.0);
    // ColorBrewer RdBu reversed: strongest (most negative) values are blue,
    // values at the upper bound are red/near white depending on the data.
    constexpr Rgb blue{33, 102, 172};
    constexpr Rgb pale{247, 247, 247};
    constexpr Rgb red{178, 24, 43};
    const Rgb result = fraction <= 0.5
        ? interpolate(blue, pale, fraction * 2.0)
        : interpolate(pale, red, (fraction - 0.5) * 2.0);
    std::ostringstream output;
    output << '#' << std::hex << std::setfill('0') << std::setw(2) << result.red
           << std::setw(2) << result.green << std::setw(2) << result.blue;
    return output.str();
}

} // namespace

auto xmlEscape(const std::string_view text) -> std::string {
    std::string result;
    result.reserve(text.size());
    for (const char symbol : text) {
        switch (symbol) {
            case '&': result += "&amp;"; break;
            case '<': result += "&lt;"; break;
            case '>': result += "&gt;"; break;
            case '"': result += "&quot;"; break;
            case '\'': result += "&apos;"; break;
            default: result.push_back(symbol); break;
        }
    }
    return result;
}

auto profileSvg(
    const std::span<const ProfilePoint> points,
    const ProfileSvgOptions& options) -> std::expected<std::string, std::string> {
    auto dimensions = validDimensions(options.width, options.height);
    if (!dimensions) {
        return std::unexpected(dimensions.error());
    }
    if (points.empty()) {
        return std::unexpected("profile requires at least one point");
    }
    double xMinimum = std::numeric_limits<double>::infinity();
    double xMaximum = -std::numeric_limits<double>::infinity();
    double yMinimum = std::numeric_limits<double>::infinity();
    double yMaximum = -std::numeric_limits<double>::infinity();
    std::size_t present{};
    for (const auto& point : points) {
        if (!std::isfinite(point.position)) {
            return std::unexpected("profile positions must be finite");
        }
        xMinimum = std::min(xMinimum, point.position);
        xMaximum = std::max(xMaximum, point.position);
        if (point.value) {
            if (!std::isfinite(*point.value)) {
                return std::unexpected("profile values must be finite or missing");
            }
            yMinimum = std::min(yMinimum, *point.value);
            yMaximum = std::max(yMaximum, *point.value);
            ++present;
        }
    }
    if (present == 0U) {
        return std::unexpected("profile contains no finite values");
    }
    if (yMinimum == yMaximum) {
        const double padding = std::max(1.0, std::abs(yMinimum) * 0.1);
        yMinimum -= padding;
        yMaximum += padding;
    }
    if (options.zeroLine) {
        yMinimum = std::min(yMinimum, 0.0);
        yMaximum = std::max(yMaximum, 0.0);
    }

    const PlotBox plot{88.0, 54.0, static_cast<double>(options.width) - 120.0,
                       static_cast<double>(options.height) - 132.0};
    std::ostringstream output;
    output << documentStart(options.width, options.height, options.title)
           << "<text class=\"heading\" x=\"" << options.width / 2U
           << "\" y=\"28\" text-anchor=\"middle\">" << xmlEscape(options.title)
           << "</text>\n";
    for (int tick = 0; tick <= 5; ++tick) {
        const double x = plot.left + plot.width * static_cast<double>(tick) / 5.0;
        const double y = plot.top + plot.height * static_cast<double>(tick) / 5.0;
        const double xValue = tickValue(xMinimum, xMaximum, tick);
        const double yValue = tickValue(yMaximum, yMinimum, tick);
        output << "<line class=\"grid\" x1=\"" << number(x) << "\" y1=\""
               << number(plot.top) << "\" x2=\"" << number(x) << "\" y2=\""
               << number(plot.top + plot.height) << "\"/>\n"
               << "<line class=\"grid\" x1=\"" << number(plot.left) << "\" y1=\""
               << number(y) << "\" x2=\"" << number(plot.left + plot.width)
               << "\" y2=\"" << number(y) << "\"/>\n"
               << "<text class=\"tick\" x=\"" << number(x) << "\" y=\""
               << number(plot.top + plot.height + 18.0) << "\" text-anchor=\"middle\">"
               << number(xValue) << "</text>\n"
               << "<text class=\"tick\" x=\"" << number(plot.left - 8.0) << "\" y=\""
               << number(y + 4.0) << "\" text-anchor=\"end\">" << number(yValue)
               << "</text>\n";
    }
    output << "<line class=\"axis\" x1=\"" << number(plot.left) << "\" y1=\""
           << number(plot.top + plot.height) << "\" x2=\""
           << number(plot.left + plot.width) << "\" y2=\""
           << number(plot.top + plot.height) << "\"/>\n"
           << "<line class=\"axis\" x1=\"" << number(plot.left) << "\" y1=\""
           << number(plot.top) << "\" x2=\"" << number(plot.left) << "\" y2=\""
           << number(plot.top + plot.height) << "\"/>\n";
    if (options.zeroLine && yMinimum <= 0.0 && yMaximum >= 0.0) {
        const double zero = scale(0.0, yMinimum, yMaximum, plot.top + plot.height, plot.top);
        output << "<line x1=\"" << number(plot.left) << "\" y1=\"" << number(zero)
               << "\" x2=\"" << number(plot.left + plot.width) << "\" y2=\""
               << number(zero) << "\" stroke=\"#b2182b\" stroke-width=\"1\""
                  " stroke-dasharray=\"5 4\"/>\n";
    }

    bool inSegment{};
    for (const auto& point : points) {
        if (!point.value) {
            if (inSegment) {
                output << "\"/>\n";
                inSegment = false;
            }
            continue;
        }
        const double x = scale(point.position, xMinimum, xMaximum, plot.left, plot.left + plot.width);
        const double y = scale(*point.value, yMinimum, yMaximum, plot.top + plot.height, plot.top);
        if (!inSegment) {
            output << "<polyline fill=\"none\" stroke=\"" << xmlEscape(options.stroke)
                   << "\" stroke-width=\"2\" stroke-linejoin=\"round\" points=\"";
            inSegment = true;
        } else {
            output << ' ';
        }
        output << number(x) << ',' << number(y);
    }
    if (inSegment) {
        output << "\"/>\n";
    }
    // A single-point segment has no visible stroke in SVG, so draw every
    // present observation as a small point as well as part of the polyline.
    for (const auto& point : points) {
        if (!point.value) {
            continue;
        }
        const double x = scale(point.position, xMinimum, xMaximum, plot.left, plot.left + plot.width);
        const double y = scale(*point.value, yMinimum, yMaximum, plot.top + plot.height, plot.top);
        output << "<circle cx=\"" << number(x) << "\" cy=\"" << number(y)
               << "\" r=\"2.2\" fill=\"" << xmlEscape(options.stroke) << "\"><title>"
               << number(point.position) << ": " << number(*point.value)
               << "</title></circle>\n";
    }
    output << "<text class=\"label\" x=\"" << number(plot.left + plot.width / 2.0)
           << "\" y=\"" << options.height - 18U << "\" text-anchor=\"middle\">"
           << xmlEscape(options.xLabel) << "</text>\n"
           << "<text class=\"label\" transform=\"translate(18 "
           << number(plot.top + plot.height / 2.0) << ") rotate(-90)\" text-anchor=\"middle\">"
           << xmlEscape(options.yLabel) << "</text>\n</svg>\n";
    return output.str();
}

auto heatmapSvg(
    const HeatmapData& data,
    const HeatmapSvgOptions& options) -> std::expected<std::string, std::string> {
    auto dimensions = validDimensions(options.width, options.height);
    if (!dimensions) {
        return std::unexpected(dimensions.error());
    }
    if (data.xLabels.empty() || data.yLabels.empty()) {
        return std::unexpected("heatmap requires non-empty x and y labels");
    }
    if (data.values.size() != data.xLabels.size() * data.yLabels.size()) {
        return std::unexpected("heatmap value count does not match its dimensions");
    }
    double minimum = std::numeric_limits<double>::infinity();
    double maximum = -std::numeric_limits<double>::infinity();
    std::vector<std::optional<double>> values = data.values;
    for (auto& value : values) {
        if (!value) {
            continue;
        }
        if (!std::isfinite(*value)) {
            return std::unexpected("heatmap values must be finite or missing");
        }
        if (options.clampPositiveToZero && *value > 0.0) {
            *value = 0.0;
        }
        minimum = std::min(minimum, *value);
        maximum = std::max(maximum, *value);
    }
    if (!std::isfinite(minimum)) {
        return std::unexpected("heatmap contains no finite values");
    }
    const PlotBox plot{108.0, 56.0, static_cast<double>(options.width) - 172.0,
                       static_cast<double>(options.height) - 144.0};
    const double cellWidth = plot.width / static_cast<double>(data.xLabels.size());
    const double cellHeight = plot.height / static_cast<double>(data.yLabels.size());

    std::ostringstream output;
    output << documentStart(options.width, options.height, options.title)
           << "<text class=\"heading\" x=\"" << options.width / 2U
           << "\" y=\"28\" text-anchor=\"middle\">" << xmlEscape(options.title)
           << "</text>\n";
    for (std::size_t y = 0; y < data.yLabels.size(); ++y) {
        for (std::size_t x = 0; x < data.xLabels.size(); ++x) {
            const auto& value = values[y * data.xLabels.size() + x];
            output << "<rect x=\"" << number(plot.left + cellWidth * static_cast<double>(x))
                   << "\" y=\"" << number(plot.top + cellHeight * static_cast<double>(y))
                   << "\" width=\"" << number(cellWidth + 0.01) << "\" height=\""
                   << number(cellHeight + 0.01) << "\" fill=\""
                   << (value ? color(*value, minimum, maximum) : xmlEscape(options.missingColor))
                   << "\"><title>" << xmlEscape(data.xLabels[x]) << " × "
                   << xmlEscape(data.yLabels[y]) << ": "
                   << (value ? number(*value) : std::string{"NA"}) << "</title></rect>\n";
        }
    }
    output << "<rect x=\"" << number(plot.left) << "\" y=\"" << number(plot.top)
           << "\" width=\"" << number(plot.width) << "\" height=\"" << number(plot.height)
           << "\" fill=\"none\" class=\"axis\"/>\n";

    const auto labelStride = [](const std::size_t count) {
        return std::max<std::size_t>(1U, (count + 11U) / 12U);
    };
    for (std::size_t x = 0; x < data.xLabels.size(); x += labelStride(data.xLabels.size())) {
        output << "<text class=\"tick\" x=\""
               << number(plot.left + (static_cast<double>(x) + 0.5) * cellWidth)
               << "\" y=\"" << number(plot.top + plot.height + 17.0)
               << "\" text-anchor=\"middle\">" << xmlEscape(data.xLabels[x]) << "</text>\n";
    }
    for (std::size_t y = 0; y < data.yLabels.size(); y += labelStride(data.yLabels.size())) {
        output << "<text class=\"tick\" x=\"" << number(plot.left - 7.0) << "\" y=\""
               << number(plot.top + (static_cast<double>(y) + 0.5) * cellHeight + 4.0)
               << "\" text-anchor=\"end\">" << xmlEscape(data.yLabels[y]) << "</text>\n";
    }
    output << "<text class=\"label\" x=\"" << number(plot.left + plot.width / 2.0)
           << "\" y=\"" << options.height - 17U << "\" text-anchor=\"middle\">"
           << xmlEscape(options.xLabel) << "</text>\n"
           << "<text class=\"label\" transform=\"translate(19 "
           << number(plot.top + plot.height / 2.0) << ") rotate(-90)\" text-anchor=\"middle\">"
           << xmlEscape(options.yLabel) << "</text>\n"
           << "<text class=\"tick\" x=\"" << options.width - 58U << "\" y=\""
           << number(plot.top + 12.0) << "\">min " << number(minimum) << "</text>\n"
           << "<text class=\"tick\" x=\"" << options.width - 58U << "\" y=\""
           << number(plot.top + 28.0) << "\">max " << number(maximum) << "</text>\n"
           << "</svg>\n";
    return output.str();
}

auto regionsSvg(
    const std::span<const RegionSpan> regions,
    const RegionSvgOptions& options) -> std::expected<std::string, std::string> {
    auto dimensions = validDimensions(options.width, options.height);
    if (!dimensions) {
        return std::unexpected(dimensions.error());
    }
    if (regions.empty()) {
        return std::unexpected("region plot requires at least one span");
    }
    long long minimum = std::numeric_limits<long long>::max();
    long long maximum = std::numeric_limits<long long>::min();
    for (const auto& region : regions) {
        if (region.id.empty()) {
            return std::unexpected("region identifiers must not be empty");
        }
        if (region.start > region.end) {
            return std::unexpected(
                "region '" + region.id + "' starts after it ends");
        }
        minimum = std::min(minimum, region.start);
        maximum = std::max(maximum, region.end);
    }
    const PlotBox plot{150.0, 55.0, static_cast<double>(options.width) - 185.0,
                       static_cast<double>(options.height) - 128.0};
    const double laneHeight = plot.height / static_cast<double>(regions.size());
    const double barHeight = std::clamp(laneHeight * 0.58, 2.0, 18.0);
    const auto x = [&](const long long position) {
        return scale(
            static_cast<double>(position), static_cast<double>(minimum),
            static_cast<double>(maximum), plot.left, plot.left + plot.width);
    };

    std::ostringstream output;
    output << documentStart(options.width, options.height, options.title)
           << "<text class=\"heading\" x=\"" << options.width / 2U
           << "\" y=\"28\" text-anchor=\"middle\">" << xmlEscape(options.title)
           << "</text>\n";
    for (int tick = 0; tick <= 5; ++tick) {
        const double coordinate = plot.left + plot.width * static_cast<double>(tick) / 5.0;
        const double value = tickValue(
            static_cast<double>(minimum), static_cast<double>(maximum), tick);
        output << "<line class=\"grid\" x1=\"" << number(coordinate) << "\" y1=\""
               << number(plot.top) << "\" x2=\"" << number(coordinate) << "\" y2=\""
               << number(plot.top + plot.height) << "\"/>\n"
               << "<text class=\"tick\" x=\"" << number(coordinate) << "\" y=\""
               << number(plot.top + plot.height + 18.0) << "\" text-anchor=\"middle\">"
               << number(value) << "</text>\n";
    }
    for (std::size_t index = 0U; index < regions.size(); ++index) {
        const auto& region = regions[index];
        const double center = plot.top + (static_cast<double>(index) + 0.5) * laneHeight;
        const double left = x(region.start);
        const double right = x(region.end);
        output << "<line x1=\"" << number(plot.left) << "\" y1=\"" << number(center)
               << "\" x2=\"" << number(plot.left + plot.width) << "\" y2=\""
               << number(center) << "\" stroke=\"#eeeeee\"/>\n"
               << "<text class=\"tick\" x=\"" << number(plot.left - 8.0) << "\" y=\""
               << number(center + 4.0) << "\" text-anchor=\"end\">"
               << xmlEscape(region.id) << "</text>\n"
               << "<rect x=\"" << number(left) << "\" y=\""
               << number(center - barHeight / 2.0) << "\" width=\""
               << number(std::max(2.0, right - left)) << "\" height=\""
               << number(barHeight) << "\" rx=\"1\" fill=\"" << xmlEscape(options.fill)
               << "\"><title>" << xmlEscape(region.id) << ": " << region.start << "–"
               << region.end << "</title></rect>\n";
    }
    output << "<rect x=\"" << number(plot.left) << "\" y=\"" << number(plot.top)
           << "\" width=\"" << number(plot.width) << "\" height=\"" << number(plot.height)
           << "\" fill=\"none\" class=\"axis\"/>\n"
           << "<text class=\"label\" x=\"" << number(plot.left + plot.width / 2.0)
           << "\" y=\"" << options.height - 17U << "\" text-anchor=\"middle\">"
           << xmlEscape(options.xLabel) << "</text>\n</svg>\n";
    return output.str();
}

} // namespace intarnanew::tools
