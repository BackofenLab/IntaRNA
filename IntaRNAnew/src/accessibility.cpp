#include "intarnanew/accessibility.hpp"

#include "intarnanew/compression.hpp"

#include "noncrossing_partition.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cctype>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string_view>

namespace intarnanew {
namespace {

inline constexpr double logZero = -std::numeric_limits<double>::infinity();
inline constexpr std::size_t maximumAccessibilityFileBytes = 256U * 1024U * 1024U;

[[nodiscard]] auto parseConstraint(
    const Sequence& sequence,
    const std::string_view encoding) -> std::vector<char> {
    std::vector<char> result(sequence.size(), '.');
    if (encoding.empty()) return result;
    if (encoding.size() == sequence.size() && encoding.find(':') == std::string_view::npos) {
        for (Index index = 0; index < sequence.size(); ++index) {
            const auto symbol = static_cast<char>(std::tolower(static_cast<unsigned char>(encoding[index])));
            if (symbol != '.' && symbol != 'x' && symbol != 'p' && symbol != 'b') {
                throw std::invalid_argument("invalid accessibility constraint symbol");
            }
            result[index] = symbol;
        }
        return result;
    }

    char kind = '.';
    std::size_t cursor{};
    while (cursor < encoding.size()) {
        while (cursor < encoding.size() && (encoding[cursor] == ',' || std::isspace(
                   static_cast<unsigned char>(encoding[cursor])) != 0)) ++cursor;
        if (cursor >= encoding.size()) break;
        if (cursor + 1U < encoding.size() && encoding[cursor + 1U] == ':') {
            kind = static_cast<char>(std::tolower(static_cast<unsigned char>(encoding[cursor])));
            if (kind != 'x' && kind != 'p' && kind != 'b') {
                throw std::invalid_argument("invalid range accessibility constraint type");
            }
            cursor += 2U;
        }
        const auto start = cursor;
        while (cursor < encoding.size() && encoding[cursor] != ',' && encoding[cursor] != ':') ++cursor;
        auto token = encoding.substr(start, cursor - start);
        while (!token.empty() && std::isspace(static_cast<unsigned char>(token.back())) != 0) {
            token.remove_suffix(1);
        }
        // Start after a possible sign of the first coordinate. This also
        // handles encodings such as -5--3.
        const auto dash = token.find('-', 1U);
        if (dash == std::string_view::npos) {
            throw std::invalid_argument("constraint range must use FROM-TO encoding");
        }
        long long first{};
        long long last{};
        const auto firstText = token.substr(0, dash);
        const auto lastText = token.substr(dash + 1U);
        const auto firstResult = std::from_chars(firstText.data(), firstText.data() + firstText.size(), first);
        const auto lastResult = std::from_chars(lastText.data(), lastText.data() + lastText.size(), last);
        if (firstResult.ec != std::errc{} || firstResult.ptr != firstText.data() + firstText.size() ||
            lastResult.ec != std::errc{} || lastResult.ptr != lastText.data() + lastText.size() ||
            first > last) {
            throw std::invalid_argument("invalid accessibility constraint range");
        }
        const auto internalFirst = sequence.internalIndex(first);
        const auto internalLast = sequence.internalIndex(last);
        if (!internalFirst || !internalLast) {
            throw std::invalid_argument("accessibility constraint range is outside the sequence");
        }
        std::fill(result.begin() + static_cast<std::ptrdiff_t>(*internalFirst),
                  result.begin() + static_cast<std::ptrdiff_t>(*internalLast + 1U), kind);
    }
    return result;
}

[[nodiscard]] auto pairEnergy(const char left, const char right) noexcept -> double {
    const auto a = static_cast<char>(std::toupper(static_cast<unsigned char>(left)));
    const auto b = static_cast<char>(std::toupper(static_cast<unsigned char>(right)));
    if ((a == 'G' && b == 'C') || (a == 'C' && b == 'G')) return -2.4;
    if ((a == 'A' && (b == 'U' || b == 'T')) ||
        ((a == 'U' || a == 'T') && b == 'A')) return -1.1;
    if ((a == 'G' && (b == 'U' || b == 'T')) ||
        ((a == 'U' || a == 'T') && b == 'G')) return -0.5;
    return -0.8; // conservative value for an ambiguous IUPAC-compatible pair
}

[[nodiscard]] auto intervalValid(const Interval interval, const Index length) noexcept -> bool {
    return interval.begin <= interval.end && interval.end < length;
}

[[nodiscard]] auto parseNumericRow(const std::string& line) -> std::vector<double> {
    std::vector<double> values;
    std::istringstream input(line);
    double value{};
    while (input >> value) values.push_back(value);
    return values;
}

[[nodiscard]] auto hasGzipSuffix(const std::string_view path) noexcept -> bool {
    if (path.size() < 3U) return false;
    const auto offset = path.size() - 3U;
    return path[offset] == '.' &&
        static_cast<char>(std::tolower(static_cast<unsigned char>(path[offset + 1U]))) == 'g' &&
        static_cast<char>(std::tolower(static_cast<unsigned char>(path[offset + 2U]))) == 'z';
}

[[nodiscard]] auto readBoundedAccessibilityFile(const std::string& path) -> std::string {
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::invalid_argument("cannot open accessibility file '" + path + "'");

    std::string bytes;
    std::array<char, 65536U> buffer{};
    while (input) {
        input.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
        const auto count = input.gcount();
        if (count < 0) throw std::invalid_argument("failed to read accessibility file '" + path + "'");
        const auto amount = static_cast<std::size_t>(count);
        if (amount > maximumAccessibilityFileBytes - bytes.size()) {
            throw std::invalid_argument("accessibility file exceeds the 256 MiB input-byte limit: '" +
                                        path + "'");
        }
        bytes.append(buffer.data(), amount);
    }
    if (!input.eof()) throw std::invalid_argument("failed to read accessibility file '" + path + "'");

    if (hasGzipMagic(bytes)) {
        GzipLimits limits;
        limits.maxCompressedBytes = maximumAccessibilityFileBytes;
        limits.maxDecompressedBytes = maximumAccessibilityFileBytes;
        auto decoded = gzipDecompress(bytes, limits);
        if (!decoded) {
            throw std::invalid_argument("cannot decode accessibility file '" + path + "': " +
                                        decoded.error());
        }
        return std::move(*decoded);
    }
    if (hasGzipSuffix(path)) {
        throw std::invalid_argument("accessibility file '" + path +
                                    "' has a .gz suffix but no gzip signature");
    }
    return bytes;
}

[[nodiscard]] auto makeConfiguredFoldingEnsemble(
    const Sequence& sequence,
    const EnergyKind energy,
    const FoldingOptions& options) -> std::unique_ptr<FoldingEnsemble> {
    return energy == EnergyKind::basePair
        ? makeBasePairFoldingEnsemble(sequence, options)
        : makeTurnerFoldingEnsemble(sequence, options);
}

// Eall/Zall describe the free monomer, not the accessibility calculation.
// In particular, they are independent of acc=N/C/P/E, accW/accL, interval
// constraints, and accessibility-only noLP/noGUend/SHAPE settings.
[[nodiscard]] auto makeMonomerSummary(
    const Sequence& sequence,
    const Config& globalConfig) -> std::unique_ptr<FoldingEnsemble> {
    FoldingOptions options;
    options.temperatureCelsius = globalConfig.temperatureCelsius;
    options.parameterSet = globalConfig.energyParameters;
    options.maximumPairSpan = sequence.size();
    options.includeDangles = !globalConfig.noDangles;
    return makeConfiguredFoldingEnsemble(sequence, globalConfig.energy, options);
}

} // namespace

TurnerAccessibility::TurnerAccessibility(
    const Sequence& sequence,
    const SideConfig& config,
    const Config& globalConfig)
    : length_(sequence.size()),
      stride_(sequence.size() + 1U),
      windowLength_(config.accessibilityWindow == 0U
          ? sequence.size() : std::min(sequence.size(), config.accessibilityWindow)),
      rt_(gasConstantKcal * (globalConfig.temperatureCelsius + 273.15)),
      probabilityExponent_(globalConfig.energy == EnergyKind::basePair ? 1.0 / rt_ : 1.0),
      blocked_(length_, false),
      intervalProbabilities_(length_ * stride_, std::numeric_limits<double>::quiet_NaN()) {
    if (!std::isfinite(rt_) || rt_ <= 0.0) {
        throw std::invalid_argument("Turner accessibility temperature must be above absolute zero");
    }
    if (!std::isfinite(config.partitionScale) || config.partitionScale < 1.0) {
        throw std::invalid_argument("partition scale must be finite and at least 1");
    }
    if (config.accessibilitySpan != 0U && config.accessibilityWindow != 0U &&
        config.accessibilitySpan > config.accessibilityWindow) {
        throw std::invalid_argument("accessibility pair span must not exceed its folding window");
    }
    const auto parsed = parseConstraint(sequence, config.accessibilityConstraint);
    for (Index position{}; position < length_; ++position) blocked_[position] = parsed[position] == 'b';

    auto summary = makeMonomerSummary(sequence, globalConfig);
    summary_ = std::shared_ptr<FoldingEnsemble>(std::move(summary));

    // The public energy=B contract is one global Nussinov ensemble. accW
    // limits the largest reported/usable interval, but does not localize the
    // fold; accL, accNoLP, and accNoGUend are unsupported for this model.
    if (globalConfig.energy == EnergyKind::basePair) {
        FoldingOptions options;
        options.temperatureCelsius = globalConfig.temperatureCelsius;
        options.maximumPairSpan = length_;
        options.partitionScale = config.partitionScale;
        options.constraint = config.accessibilityConstraint;
        std::shared_ptr<FoldingEnsemble> ensemble;
        if (config.accessibilityConstraint.empty()) {
            ensemble = summary_;
        } else {
            auto owned = makeBasePairFoldingEnsemble(sequence, options);
            ensemble = std::shared_ptr<FoldingEnsemble>(std::move(owned));
        }
        windows_.push_back(Window{0U, length_, std::move(ensemble)});
        return;
    }

    const auto windowCount = length_ <= windowLength_ ? 1U : length_ - windowLength_ + 1U;
    windows_.reserve(windowCount);
    for (Index begin{}; begin < windowCount; ++begin) {
        const auto end = begin + windowLength_;
        std::string constraint(parsed.begin() + static_cast<std::ptrdiff_t>(begin),
                               parsed.begin() + static_cast<std::ptrdiff_t>(end));
        FoldingOptions options;
        options.temperatureCelsius = globalConfig.temperatureCelsius;
        options.parameterSet = globalConfig.energyParameters;
        options.maximumPairSpan = config.accessibilitySpan == 0U
            ? windowLength_ : std::min(windowLength_, config.accessibilitySpan);
        options.noLonelyPairs = globalConfig.accessibilityNoLonelyPairs;
        options.noGuHelixEnds = globalConfig.accessibilityNoGuAtEnds;
        options.includeDangles = !globalConfig.noDangles;
        options.partitionScale = config.partitionScale;
        options.constraint = std::move(constraint);
        options.shapeFile = config.shapeFile;
        options.shapeMethod = config.shapeMethod;
        options.shapeConversion = config.shapeConversion;
        const Sequence subsequence(
            sequence.id(), sequence.str().substr(begin, windowLength_), sequence.externalIndex(begin));
        const auto canReuseSummary = begin == 0U && windowLength_ == length_ &&
            (config.accessibilitySpan == 0U || config.accessibilitySpan >= length_) &&
            !globalConfig.accessibilityNoLonelyPairs && !globalConfig.accessibilityNoGuAtEnds &&
            config.accessibilityConstraint.empty() && config.shapeFile.empty() &&
            config.shapeMethod.empty() && config.shapeConversion.empty();
        std::shared_ptr<FoldingEnsemble> ensemble;
        if (canReuseSummary) {
            ensemble = summary_;
        } else {
            auto owned = makeConfiguredFoldingEnsemble(subsequence, globalConfig.energy, options);
            ensemble = std::shared_ptr<FoldingEnsemble>(std::move(owned));
        }
        windows_.push_back(Window{begin, end, std::move(ensemble)});
    }
}

auto TurnerAccessibility::openingEnergy(const Interval interval) const -> Energy {
    const auto probability = unpairedProbability(interval);
    return probability <= 0.0 ? infinity : -rt_ * std::log(probability);
}

auto TurnerAccessibility::unpairedProbability(const Interval interval) const -> double {
    if (!intervalValid(interval, length_)) throw std::out_of_range("accessibility interval is out of range");
    if (std::any_of(blocked_.begin() + static_cast<std::ptrdiff_t>(interval.begin),
                    blocked_.begin() + static_cast<std::ptrdiff_t>(interval.end + 1U),
                    [](const bool value) { return value; })) return 0.0;
    if (interval.size() > windowLength_) return 0.0;
    const auto offset = interval.begin * stride_ + interval.size();
    {
        const std::scoped_lock lock(cacheMutex_);
        if (!std::isnan(intervalProbabilities_[offset])) return intervalProbabilities_[offset];
    }

    double sum{};
    Index count{};
    for (const auto& window : windows_) {
        if (window.begin > interval.begin || window.end <= interval.end) continue;
        const Interval local{interval.begin - window.begin, interval.end - window.begin};
        sum += window.ensemble->jointUnpairedProbability(local);
        ++count;
    }
    if (count == 0U) throw std::logic_error("no accessibility folding window contains a valid interval");
    const auto ensembleProbability = std::clamp(sum / static_cast<double>(count), 0.0, 1.0);
    // energy=B defines its monomer ensemble with RT=1, while accessibility ED
    // remains an energy that is converted to Pu using the physical RT. Thus
    // Pu = exp(-ED/RT_phys) = (Z_constrained/Z)^(1/RT_phys).
    double probability = ensembleProbability;
    if (ensembleProbability > 0.0 && probabilityExponent_ != 1.0) {
        // The legacy public energy interface stores ED in integer centikcal
        // units. Preserve that observable boundary before converting ED to Pu.
        const auto opening = std::trunc(-std::log(ensembleProbability) * 100.0) / 100.0;
        probability = std::exp(-opening * probabilityExponent_);
    }
    probability = std::clamp(probability, 0.0, 1.0);
    {
        const std::scoped_lock lock(cacheMutex_);
        intervalProbabilities_[offset] = probability;
    }
    return probability;
}

auto TurnerAccessibility::positionUnpairedProbability(const Index position) const -> double {
    if (position >= length_) throw std::out_of_range("accessibility position is out of range");
    return unpairedProbability({position, position});
}

auto TurnerAccessibility::blocked(const Index position) const -> bool {
    if (position >= length_) throw std::out_of_range("accessibility position is out of range");
    return blocked_[position];
}

auto TurnerAccessibility::ensembleLogPartition() const noexcept -> std::optional<double> {
    return summary_ ? std::optional<double>{summary_->logPartition()} : std::nullopt;
}

auto TurnerAccessibility::ensembleFreeEnergy() const noexcept -> std::optional<Energy> {
    return summary_ ? std::optional<Energy>{summary_->ensembleFreeEnergy()} : std::nullopt;
}

DisabledAccessibility::DisabledAccessibility(const Sequence& sequence, const std::string_view constraint)
    : length_(sequence.size()), blocked_(length_, false) {
    const auto parsed = parseConstraint(sequence, constraint);
    for (Index index = 0; index < length_; ++index) blocked_[index] = parsed[index] == 'b';
}

DisabledAccessibility::DisabledAccessibility(
    const Sequence& sequence,
    const SideConfig& side,
    const Config& globalConfig)
    : DisabledAccessibility(sequence, side.accessibilityConstraint) {
    if (!side.shapeFile.empty() || !side.shapeMethod.empty() || !side.shapeConversion.empty()) {
        throw std::invalid_argument("SHAPE data requires computed accessibility");
    }
    if (globalConfig.predictionRequirements.computeMonomerPartition) {
        ensemble_ = makeMonomerSummary(sequence, globalConfig);
    }
}

auto DisabledAccessibility::openingEnergy(const Interval interval) const -> Energy {
    if (!intervalValid(interval, length_)) throw std::out_of_range("accessibility interval is out of range");
    return std::any_of(blocked_.begin() + static_cast<std::ptrdiff_t>(interval.begin),
                       blocked_.begin() + static_cast<std::ptrdiff_t>(interval.end + 1U),
                       [](const bool value) { return value; }) ? infinity : 0.0;
}

auto DisabledAccessibility::unpairedProbability(const Interval interval) const -> double {
    return std::isfinite(openingEnergy(interval)) ? 1.0 : 0.0;
}

auto DisabledAccessibility::positionUnpairedProbability(const Index position) const -> double {
    if (position >= length_) throw std::out_of_range("accessibility position is out of range");
    return blocked_[position] ? 0.0 : 1.0;
}

auto DisabledAccessibility::blocked(const Index position) const -> bool {
    if (position >= length_) throw std::out_of_range("accessibility position is out of range");
    return blocked_[position];
}

auto DisabledAccessibility::unconstrained(const Interval interval) const noexcept -> bool {
    return intervalValid(interval, length_) &&
        std::none_of(blocked_.begin() + static_cast<std::ptrdiff_t>(interval.begin),
                     blocked_.begin() + static_cast<std::ptrdiff_t>(interval.end + 1U),
                     [](const bool value) { return value; });
}

auto DisabledAccessibility::ensembleLogPartition() const noexcept -> std::optional<double> {
    return ensemble_ ? std::optional<double>{ensemble_->logPartition()} : std::nullopt;
}

auto DisabledAccessibility::ensembleFreeEnergy() const noexcept -> std::optional<Energy> {
    return ensemble_ ? std::optional<Energy>{ensemble_->ensembleFreeEnergy()} : std::nullopt;
}

NativeAccessibility::NativeAccessibility(
    const Sequence& sequence,
    const SideConfig& config,
    const double temperatureCelsius)
    : length_(sequence.size()),
      stride_(sequence.size() + 1U),
      maxPairSpan_(config.accessibilitySpan == 0U
          ? sequence.size()
          : std::min(sequence.size(), config.accessibilitySpan)),
      rt_(gasConstantKcal * (temperatureCelsius + 273.15)),
      sequence_(sequence.str()),
      constraints_(parseConstraint(sequence, config.accessibilityConstraint)),
      blocked_(length_, false),
      intervalProbabilities_(length_ * stride_, std::numeric_limits<double>::quiet_NaN()) {
    if (!std::isfinite(rt_) || rt_ <= 0.0) {
        throw std::invalid_argument("native accessibility temperature must be above absolute zero");
    }
    for (Index index = 0; index < length_; ++index) {
        blocked_[index] = constraints_[index] == 'b';
    }
    constexpr Index minimumHairpinDistance = 4U;
    auto partitions = detail::computeNoncrossingIntervalPartitions(
        length_,
        [this](const Index position) { return constraints_[position] != 'p'; },
        [this](const Index left, const Index right) {
            if (right < left + minimumHairpinDistance || right - left > maxPairSpan_ ||
                constraints_[left] == 'x' || constraints_[left] == 'b' ||
                constraints_[right] == 'x' || constraints_[right] == 'b' ||
                !canPair(sequence_[left], sequence_[right])) {
                return logZero;
            }
            return -pairEnergy(sequence_[left], sequence_[right]) / rt_;
        });
    logPartition_ = partitions.logPartition;
    intervalProbabilities_ = std::move(partitions.probabilities);
    if (!std::isfinite(logPartition_)) {
        throw std::invalid_argument("native accessibility constraints admit no secondary structure");
    }
}


auto NativeAccessibility::cacheOffset(const Interval interval) const -> Index {
    return interval.begin * stride_ + interval.size();
}

auto NativeAccessibility::openingEnergy(const Interval interval) const -> Energy {
    const auto probability = unpairedProbability(interval);
    return probability <= 0.0 ? infinity : -rt_ * std::log(probability);
}

auto NativeAccessibility::unpairedProbability(const Interval interval) const -> double {
    if (!intervalValid(interval, length_)) throw std::out_of_range("accessibility interval is out of range");
    if (std::any_of(blocked_.begin() + static_cast<std::ptrdiff_t>(interval.begin),
                    blocked_.begin() + static_cast<std::ptrdiff_t>(interval.end + 1U),
                    [](const bool value) { return value; })) {
        return 0.0;
    }
    for (Index position = interval.begin; position <= interval.end; ++position) {
        if (constraints_[position] == 'p') return 0.0;
    }
    return intervalProbabilities_[cacheOffset(interval)];
}

auto NativeAccessibility::positionUnpairedProbability(const Index position) const -> double {
    if (position >= length_) throw std::out_of_range("accessibility position is out of range");
    return unpairedProbability({position, position});
}

auto NativeAccessibility::blocked(const Index position) const -> bool {
    if (position >= length_) throw std::out_of_range("accessibility position is out of range");
    return blocked_[position];
}

TableAccessibility::TableAccessibility(
    const Sequence& sequence,
    const SideConfig& config,
    const double temperatureCelsius,
    const bool tableContainsProbabilities)
    : length_(sequence.size()),
      stride_(sequence.size() + 1U),
      rt_(gasConstantKcal * (temperatureCelsius + 273.15)),
      probabilities_(length_ * stride_, std::numeric_limits<double>::quiet_NaN()),
      blocked_(length_, false) {
    const auto constraints = parseConstraint(sequence, config.accessibilityConstraint);
    for (Index index = 0; index < length_; ++index) blocked_[index] = constraints[index] == 'b';

    if (config.accessibilityFile == "STDIN" || config.accessibilityFile == "-") {
        throw std::invalid_argument("table accessibility from STDIN must be supplied through a dedicated stream");
    }
    std::istringstream input(readBoundedAccessibilityFile(config.accessibilityFile));

    std::string line;
    Index implicitRow{};
    while (std::getline(input, line)) {
        const auto first = line.find_first_not_of(" \t\r");
        if (first == std::string::npos || line[first] == '#' || line[first] == '>') continue;
        auto values = parseNumericRow(line);
        if (values.empty()) continue;

        Index row = implicitRow;
        Index valueStart{};
        const auto possibleIndex = static_cast<long long>(std::llround(values.front()));
        if (values.size() > 1U && std::abs(values.front() - static_cast<double>(possibleIndex)) < 1e-9 &&
            possibleIndex >= sequence.firstPosition() &&
            possibleIndex < sequence.firstPosition() + static_cast<long long>(sequence.size())) {
            row = static_cast<Index>(possibleIndex - sequence.firstPosition());
            valueStart = 1U;
        }
        if (row >= length_) continue;
        for (Index column = valueStart; column < values.size(); ++column) {
            const Index intervalLength = column - valueStart + 1U;
            if (row + 1U < intervalLength) continue;
            const Interval interval{row + 1U - intervalLength, row};
            const auto raw = values[column];
            const auto probability = tableContainsProbabilities
                ? raw
                : (std::isfinite(raw) ? std::exp(-raw / rt_) : 0.0);
            probabilities_[offset(interval)] = std::clamp(probability, 0.0, 1.0);
        }
        ++implicitRow;
    }
    if (!input.eof()) {
        throw std::invalid_argument("failed while parsing accessibility file '" +
                                    config.accessibilityFile + "'");
    }
}

TableAccessibility::TableAccessibility(
    const Sequence& sequence,
    const SideConfig& config,
    const Config& globalConfig,
    const bool tableContainsProbabilities)
    : TableAccessibility(
          sequence, config, globalConfig.temperatureCelsius, tableContainsProbabilities) {
    if (!config.shapeFile.empty() || !config.shapeMethod.empty() || !config.shapeConversion.empty()) {
        throw std::invalid_argument("SHAPE data requires computed accessibility");
    }
    if (globalConfig.predictionRequirements.computeMonomerPartition) {
        ensemble_ = makeMonomerSummary(sequence, globalConfig);
    }
}

auto TableAccessibility::offset(const Interval interval) const -> Index {
    if (!intervalValid(interval, length_)) throw std::out_of_range("accessibility interval is out of range");
    return interval.begin * stride_ + interval.size();
}

auto TableAccessibility::openingEnergy(const Interval interval) const -> Energy {
    const auto probability = unpairedProbability(interval);
    return probability <= 0.0 ? infinity : -rt_ * std::log(probability);
}

auto TableAccessibility::unpairedProbability(const Interval interval) const -> double {
    if (!intervalValid(interval, length_)) throw std::out_of_range("accessibility interval is out of range");
    if (std::any_of(blocked_.begin() + static_cast<std::ptrdiff_t>(interval.begin),
                    blocked_.begin() + static_cast<std::ptrdiff_t>(interval.end + 1U),
                    [](const bool value) { return value; })) return 0.0;
    const auto stored = probabilities_[offset(interval)];
    if (std::isfinite(stored)) return stored;
    throw std::runtime_error("accessibility table does not contain the requested joint interval probability");
}

auto TableAccessibility::positionUnpairedProbability(const Index position) const -> double {
    return unpairedProbability({position, position});
}

auto TableAccessibility::blocked(const Index position) const -> bool {
    if (position >= length_) throw std::out_of_range("accessibility position is out of range");
    return blocked_[position];
}

auto TableAccessibility::ensembleLogPartition() const noexcept -> std::optional<double> {
    return ensemble_ ? std::optional<double>{ensemble_->logPartition()} : std::nullopt;
}

auto TableAccessibility::ensembleFreeEnergy() const noexcept -> std::optional<Energy> {
    return ensemble_ ? std::optional<Energy>{ensemble_->ensembleFreeEnergy()} : std::nullopt;
}

auto makeAccessibility(
    const Sequence& sequence,
    const SideConfig& config,
    const double temperatureCelsius) -> std::expected<std::unique_ptr<AccessibilityProvider>, std::string> {
    Config globalConfig;
    globalConfig.temperatureCelsius = temperatureCelsius;
    return makeAccessibility(sequence, config, globalConfig);
}

auto makeAccessibility(
    const Sequence& sequence,
    const SideConfig& config,
    const Config& globalConfig) -> std::expected<std::unique_ptr<AccessibilityProvider>, std::string> {
    try {
        switch (config.accessibility) {
            case AccessibilityKind::disabled:
                return std::make_unique<DisabledAccessibility>(sequence, config, globalConfig);
            case AccessibilityKind::compute:
                return std::make_unique<TurnerAccessibility>(sequence, config, globalConfig);
            case AccessibilityKind::probabilitiesFile:
                return std::make_unique<TableAccessibility>(
                    sequence, config, globalConfig, true);
            case AccessibilityKind::energiesFile:
                return std::make_unique<TableAccessibility>(
                    sequence, config, globalConfig, false);
        }
    } catch (const std::exception& exception) {
        return std::unexpected(exception.what());
    }
    return std::unexpected("unsupported accessibility mode");
}

} // namespace intarnanew
