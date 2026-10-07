#include "intarnanew/folding.hpp"

#include "folding_parameters.hpp"
#include "noncrossing_partition.hpp"

#include <algorithm>
#include <array>
#include <cerrno>
#include <charconv>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <mutex>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace intarnanew {
namespace {

using folding_detail::FoldingParameters;
inline constexpr double logZero = -std::numeric_limits<double>::infinity();
inline constexpr Index minimumHairpinUnpaired = 3U;

enum class PairType : std::size_t { cg, gc, gu, ug, au, ua, ambiguous };

struct ShapeTerms {
    std::vector<Energy> unpaired;
    std::vector<Energy> paired;
    std::vector<Energy> stack;
};

struct ShapeEncoding {
    char method{'D'};
    double first{};
    double second{};
    char conversion{'O'};
    double conversionFirst{};
    double conversionSecond{};
};

[[nodiscard]] auto logAdd(const double left, const double right) noexcept -> double {
    if (left == logZero) return right;
    if (right == logZero) return left;
    const auto high = std::max(left, right);
    return high + std::log1p(std::exp(std::min(left, right) - high));
}

[[nodiscard]] auto pairType(const char left, const char right) noexcept -> PairType {
    const auto first = static_cast<char>(std::toupper(static_cast<unsigned char>(left)));
    const auto second = static_cast<char>(std::toupper(static_cast<unsigned char>(right)));
    if (first == 'C' && second == 'G') return PairType::cg;
    if (first == 'G' && second == 'C') return PairType::gc;
    if (first == 'G' && second == 'U') return PairType::gu;
    if (first == 'U' && second == 'G') return PairType::ug;
    if (first == 'A' && second == 'U') return PairType::au;
    if (first == 'U' && second == 'A') return PairType::ua;
    return PairType::ambiguous;
}

[[nodiscard]] constexpr auto pairIndex(const PairType type) noexcept -> std::size_t {
    return static_cast<std::size_t>(type);
}

[[nodiscard]] constexpr auto reversePair(const PairType type) noexcept -> PairType {
    constexpr std::array reverse{
        PairType::gc, PairType::cg, PairType::ug, PairType::gu,
        PairType::ua, PairType::au, PairType::ambiguous,
    };
    return reverse[pairIndex(type)];
}

[[nodiscard]] constexpr auto weakPair(const PairType type) noexcept -> bool {
    return type == PairType::gu || type == PairType::ug ||
           type == PairType::au || type == PairType::ua;
}

[[nodiscard]] constexpr auto guPair(const PairType type) noexcept -> bool {
    return type == PairType::gu || type == PairType::ug;
}

[[nodiscard]] auto baseIndex(const char base) noexcept -> std::size_t {
    switch (static_cast<char>(std::toupper(static_cast<unsigned char>(base)))) {
        case 'A': return 1U;
        case 'C': return 2U;
        case 'G': return 3U;
        case 'U': return 4U;
        default: return 0U;
    }
}

[[nodiscard]] auto canonicalBaseIndex(const char base) noexcept -> std::optional<std::size_t> {
    switch (static_cast<char>(std::toupper(static_cast<unsigned char>(base)))) {
        case 'A': return 0U;
        case 'C': return 1U;
        case 'G': return 2U;
        case 'U': return 3U;
        default: return std::nullopt;
    }
}

[[nodiscard]] auto tableStack(const FoldingParameters& parameters,
                              const PairType outer, const PairType inner) noexcept -> Energy {
    return parameters.stack[pairIndex(outer) * 7U + pairIndex(inner)];
}

[[nodiscard]] auto tableMismatch(const std::vector<Energy>& table, const PairType type,
                                 const char first, const char second) noexcept -> Energy {
    return table[(pairIndex(type) * 5U + baseIndex(first)) * 5U + baseIndex(second)];
}

[[nodiscard]] auto tableDangle(const std::vector<Energy>& table, const PairType type,
                               const char nucleotide) noexcept -> Energy {
    const auto value = table[pairIndex(type) * 5U + baseIndex(nucleotide)];
    return std::isfinite(value) ? std::min(0.0, value) : 0.0;
}

[[nodiscard]] auto tableInt11(const FoldingParameters& parameters, const PairType outer,
                              const PairType inner, const char first, const char second) noexcept -> Energy {
    const auto index = (((pairIndex(outer) * 7U + pairIndex(inner)) * 5U + baseIndex(first)) * 5U +
                        baseIndex(second));
    return parameters.int11[index];
}

[[nodiscard]] auto tableInt21(const FoldingParameters& parameters, const PairType outer,
                              const PairType inner, const char first, const char second,
                              const char third) noexcept -> Energy {
    const auto index = ((((pairIndex(outer) * 7U + pairIndex(inner)) * 5U + baseIndex(first)) * 5U +
                         baseIndex(second)) * 5U + baseIndex(third));
    return parameters.int21[index];
}

[[nodiscard]] auto tableInt22(const FoldingParameters& parameters, const PairType outer,
                              const PairType inner, const char first, const char second,
                              const char third, const char fourth) noexcept -> std::optional<Energy> {
    if (pairIndex(outer) >= 6U || pairIndex(inner) >= 6U) return std::nullopt;
    const auto a = canonicalBaseIndex(first);
    const auto b = canonicalBaseIndex(second);
    const auto c = canonicalBaseIndex(third);
    const auto d = canonicalBaseIndex(fourth);
    if (!a || !b || !c || !d) return std::nullopt;
    const auto index = (((((pairIndex(outer) * 6U + pairIndex(inner)) * 4U + *a) * 4U + *b) * 4U +
                          *c) * 4U + *d);
    return parameters.int22[index];
}

[[nodiscard]] auto loopInitiation(const std::vector<Energy>& table, const Index size,
                                  const Energy logarithmicSlope) noexcept -> Energy {
    if (size < table.size()) return table[size];
    return table.back() + logarithmicSlope *
        std::log(static_cast<double>(size) / static_cast<double>(table.size() - 1U));
}

[[nodiscard]] auto parseConstraint(const Sequence& sequence, const std::string_view encoding)
    -> std::vector<char> {
    std::vector<char> result(sequence.size(), '.');
    if (encoding.empty()) return result;
    if (encoding.size() == sequence.size() && encoding.find(':') == std::string_view::npos) {
        for (Index index{}; index < sequence.size(); ++index) {
            const auto symbol = static_cast<char>(std::tolower(static_cast<unsigned char>(encoding[index])));
            if (symbol != '.' && symbol != 'x' && symbol != 'p' && symbol != 'b') {
                throw std::invalid_argument("invalid accessibility constraint symbol");
            }
            result[index] = symbol;
        }
        return result;
    }

    char kind{};
    std::size_t cursor{};
    while (cursor < encoding.size()) {
        while (cursor < encoding.size() &&
               (encoding[cursor] == ',' || std::isspace(static_cast<unsigned char>(encoding[cursor])) != 0)) {
            ++cursor;
        }
        if (cursor >= encoding.size()) break;
        if (cursor + 1U >= encoding.size() || encoding[cursor + 1U] != ':') {
            throw std::invalid_argument("each accessibility range requires an x:, p:, or b: prefix");
        }
        kind = static_cast<char>(std::tolower(static_cast<unsigned char>(encoding[cursor])));
        if (kind != 'x' && kind != 'p' && kind != 'b') {
            throw std::invalid_argument("invalid range accessibility constraint type");
        }
        cursor += 2U;
        while (true) {
            const auto start = cursor;
            while (cursor < encoding.size() && encoding[cursor] != ',') ++cursor;
            auto token = encoding.substr(start, cursor - start);
            while (!token.empty() && std::isspace(static_cast<unsigned char>(token.back())) != 0) {
                token.remove_suffix(1U);
            }
            const auto dash = token.find('-', 1U);
            if (dash == std::string_view::npos) {
                throw std::invalid_argument("constraint range must use FROM-TO encoding");
            }
            long long first{};
            long long last{};
            const auto firstText = token.substr(0U, dash);
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
            if (cursor >= encoding.size()) break;
            ++cursor;
            std::size_t lookahead = cursor;
            while (lookahead < encoding.size() &&
                   std::isspace(static_cast<unsigned char>(encoding[lookahead])) != 0) ++lookahead;
            if (lookahead + 1U < encoding.size() && encoding[lookahead + 1U] == ':') {
                cursor = lookahead;
                break;
            }
        }
    }
    return result;
}

[[nodiscard]] auto parseFloating(const std::string_view text, std::size_t& cursor,
                                 const std::string_view context) -> double {
    const std::string owned(text.substr(cursor));
    char* end{};
    errno = 0;
    const auto value = std::strtod(owned.c_str(), &end);
    if (end == owned.c_str() || errno == ERANGE || !std::isfinite(value)) {
        throw std::invalid_argument("invalid numeric value in " + std::string(context));
    }
    cursor += static_cast<std::size_t>(end - owned.c_str());
    return value;
}

[[nodiscard]] auto parseTagged(const std::string_view text, const char kind,
                               const std::string_view tags, std::array<double, 2U> defaults,
                               const std::string_view context) -> std::array<double, 2U> {
    if (text.empty()) return defaults;
    if (static_cast<char>(std::toupper(static_cast<unsigned char>(text.front()))) != kind) {
        throw std::invalid_argument("invalid " + std::string(context) + " encoding");
    }
    std::array<bool, 2U> seen{};
    std::size_t cursor{1U};
    while (cursor < text.size()) {
        const auto tag = static_cast<char>(std::tolower(static_cast<unsigned char>(text[cursor++])));
        const auto position = tags.find(tag);
        if (position == std::string_view::npos || position >= defaults.size() || seen[position]) {
            throw std::invalid_argument("invalid or duplicate tag in " + std::string(context));
        }
        seen[position] = true;
        defaults[position] = parseFloating(text, cursor, context);
    }
    return defaults;
}

[[nodiscard]] auto parseShapeEncoding(const std::string_view methodText,
                                      const std::string_view conversionText) -> ShapeEncoding {
    ShapeEncoding result;
    const auto method = methodText.empty()
        ? 'D' : static_cast<char>(std::toupper(static_cast<unsigned char>(methodText.front())));
    result.method = method;
    if (method == 'D') {
        const auto values = parseTagged(methodText.empty() ? std::string_view{"D"} : methodText,
                                        'D', "mb", {1.8, -0.6}, "SHAPE method");
        result.first = values[0U];
        result.second = values[1U];
    } else if (method == 'Z') {
        const auto values = parseTagged(methodText, 'Z', "b_", {0.89, 0.0}, "SHAPE method");
        result.first = values[0U];
    } else if (method == 'W') {
        if (methodText.size() != 1U) throw std::invalid_argument("SHAPE method W takes no parameters");
    } else {
        throw std::invalid_argument("SHAPE method must be D, Z, or W");
    }

    const auto conversion = conversionText.empty()
        ? 'O' : static_cast<char>(std::toupper(static_cast<unsigned char>(conversionText.front())));
    result.conversion = conversion;
    if (conversion == 'M' || conversion == 'S') {
        if (!conversionText.empty() && conversionText.size() != 1U) {
            throw std::invalid_argument("SHAPE conversion M/S takes no parameters");
        }
    } else if (conversion == 'C') {
        std::size_t cursor{1U};
        result.conversionFirst = conversionText.size() == 1U
            ? 0.25 : parseFloating(conversionText, cursor, "SHAPE conversion");
        if (cursor != conversionText.size() || result.conversionFirst < 0.0) {
            throw std::invalid_argument("invalid SHAPE cutoff conversion");
        }
    } else if (conversion == 'L') {
        const auto values = parseTagged(conversionText, 'L', "si", {0.68, 0.2}, "SHAPE conversion");
        result.conversionFirst = values[0U];
        result.conversionSecond = values[1U];
    } else if (conversion == 'O') {
        const auto values = parseTagged(conversionText.empty() ? std::string_view{"O"} : conversionText,
                                        'O', "si", {1.6, -2.29}, "SHAPE conversion");
        result.conversionFirst = values[0U];
        result.conversionSecond = values[1U];
    } else {
        throw std::invalid_argument("SHAPE conversion must be M, C, S, L, or O");
    }
    if (method != 'Z' && !conversionText.empty()) {
        throw std::invalid_argument("SHAPE conversion is only applicable to method Z");
    }
    return result;
}

[[nodiscard]] auto readShapeValues(const Sequence& sequence, const std::string& path)
    -> std::vector<std::optional<double>> {
    if (path == "STDIN" || path == "-") {
        throw std::invalid_argument("SHAPE data from STDIN requires an explicit input stream");
    }
    std::ifstream input(path);
    if (!input) throw std::invalid_argument("cannot open SHAPE data file '" + path + "'");
    std::vector<std::optional<double>> result(sequence.size());
    std::string line;
    std::size_t lineNumber{};
    while (std::getline(input, line)) {
        ++lineNumber;
        const auto first = line.find_first_not_of(" \t\r");
        if (first == std::string::npos || line[first] == '#') continue;
        std::istringstream fields(line);
        long long position{};
        char nucleotide{};
        double value{};
        std::string extra;
        if (!(fields >> position >> nucleotide >> value) || (fields >> extra) || !std::isfinite(value)) {
            throw std::invalid_argument("malformed SHAPE row " + std::to_string(lineNumber));
        }
        const auto internal = sequence.internalIndex(position);
        // Local folding windows consume a full-sequence SHAPE file. Rows that
        // belong to another window are intentionally ignored.
        if (!internal) continue;
        const auto observed = static_cast<char>(std::toupper(static_cast<unsigned char>(nucleotide)));
        if (observed != sequence[*internal]) {
            throw std::invalid_argument("SHAPE nucleotide does not match the RNA sequence at position " +
                                        std::to_string(position));
        }
        if (result[*internal]) throw std::invalid_argument("duplicate SHAPE position " + std::to_string(position));
        if (value > -999.0) {
            result[*internal] = value;
        }
    }
    if (!input.eof()) throw std::invalid_argument("failed while reading SHAPE data file");
    return result;
}

[[nodiscard]] auto shapeTerms(const Sequence& sequence, const FoldingOptions& options) -> ShapeTerms {
    ShapeTerms result{
        std::vector<Energy>(sequence.size()),
        std::vector<Energy>(sequence.size()),
        std::vector<Energy>(sequence.size()),
    };
    if (options.shapeFile.empty()) {
        if (!options.shapeMethod.empty() || !options.shapeConversion.empty()) {
            throw std::invalid_argument("SHAPE method/conversion requires a SHAPE data file");
        }
        return result;
    }
    const auto encoding = parseShapeEncoding(options.shapeMethod, options.shapeConversion);
    const auto values = readShapeValues(sequence, options.shapeFile);
    for (Index index{}; index < sequence.size(); ++index) {
        if (!values[index]) continue;
        const auto value = *values[index];
        if (encoding.method == 'D') {
            if (value < 0.0) throw std::invalid_argument("Deigan SHAPE reactivities must be nonnegative");
            result.stack[index] = encoding.first * std::log(value + 1.0) + encoding.second;
        } else if (encoding.method == 'W') {
            result.unpaired[index] = value;
        } else {
            if (value < 0.0) throw std::invalid_argument("Zarringhalam SHAPE reactivities must be nonnegative");
            double probability{};
            switch (encoding.conversion) {
                case 'M':
                    // Piecewise-linear Zarringhalam mapping (public ViennaRNA
                    // SHAPE strategy specification).
                    if (value < 0.25) probability = 1.4 * value;
                    else if (value < 0.30) probability = 0.35 + 4.0 * (value - 0.25);
                    else if (value < 0.70) probability = 0.55 + 0.75 * (value - 0.30);
                    else probability = 0.85 + 0.5 * (value - 0.70);
                    break;
                case 'C': probability = value >= encoding.conversionFirst ? 1.0 : 0.0; break;
                case 'S': probability = value; break;
                case 'L': probability = encoding.conversionFirst * value + encoding.conversionSecond; break;
                case 'O': probability = (std::log(std::max(value, std::numeric_limits<double>::min())) -
                                         encoding.conversionSecond) / encoding.conversionFirst; break;
                default: throw std::logic_error("validated SHAPE conversion is unreachable");
            }
            probability = std::clamp(probability, 0.0, 1.0);
            result.unpaired[index] = encoding.first * probability;
            result.paired[index] = encoding.first * (1.0 - probability);
        }
    }
    return result;
}

class BasePairFoldingEnsemble final : public FoldingEnsemble {
public:
    BasePairFoldingEnsemble(const Sequence& sequence, FoldingOptions options)
        : sequence_(sequence.str()),
          options_(std::move(options)),
          constraints_(parseConstraint(sequence, options_.constraint)),
          span_(options_.maximumPairSpan == 0U
              ? sequence.size() : std::min(options_.maximumPairSpan, sequence.size())),
          cache_(sequence.size() * (sequence.size() + 1U),
                 std::numeric_limits<double>::quiet_NaN()) {
        if (!options_.shapeFile.empty() || !options_.shapeMethod.empty() ||
            !options_.shapeConversion.empty()) {
            throw std::invalid_argument("SHAPE pseudo-energies require the nearest-neighbour folding model");
        }
        if (!std::isfinite(options_.partitionScale) || options_.partitionScale < 1.0) {
            throw std::invalid_argument("partition scale must be finite and at least 1");
        }
        if (!options_.noLonelyPairs && !options_.noGuHelixEnds) {
            auto partitions = detail::computeNoncrossingIntervalPartitions(
                sequence_.size(),
                [this](const Index position) { return constraints_[position] != 'p'; },
                [this](const Index left, const Index right) {
                    return pairAllowed(left, right, std::nullopt) ? 1.0 : logZero;
                });
            logPartition_ = partitions.logPartition;
            cache_ = std::move(partitions.probabilities);
        } else {
            // noLP/noGUend make pair admissibility depend on an outer/inner
            // stack state. Preserve the stateful recurrence for these modes.
            logPartition_ = compute(std::nullopt);
        }
        if (!std::isfinite(logPartition_)) {
            throw std::invalid_argument("folding constraints admit no secondary structure");
        }
    }

    [[nodiscard]] auto logPartition() const noexcept -> double override { return logPartition_; }
    [[nodiscard]] auto ensembleFreeEnergy() const noexcept -> Energy override { return -logPartition_; }

    [[nodiscard]] auto jointUnpairedProbability(const Interval interval) const -> double override {
        if (interval.begin > interval.end || interval.end >= sequence_.size()) {
            throw std::out_of_range("folding interval is out of range");
        }
        for (Index position = interval.begin; position <= interval.end; ++position) {
            if (constraints_[position] == 'p') return 0.0;
        }
        const auto offset = interval.begin * (sequence_.size() + 1U) + interval.size();
        {
            const std::scoped_lock lock(cacheMutex_);
            if (!std::isnan(cache_[offset])) return cache_[offset];
        }
        const auto constrained = compute(interval);
        const auto probability = constrained == logZero
            ? 0.0 : std::clamp(std::exp(std::min(0.0, constrained - logPartition_)), 0.0, 1.0);
        {
            const std::scoped_lock lock(cacheMutex_);
            cache_[offset] = probability;
        }
        return probability;
    }

private:
    [[nodiscard]] auto matrixOffset(const Index left, const Index right) const noexcept -> Index {
        return left * sequence_.size() + right;
    }

    [[nodiscard]] auto intervalOffset(const Index begin, const Index end) const noexcept -> Index {
        return begin * (sequence_.size() + 1U) + end;
    }

    [[nodiscard]] auto pairAllowed(const Index left, const Index right,
                                   const std::optional<Interval> forcedUnpaired) const noexcept -> bool {
        return right >= left + minimumHairpinUnpaired + 1U &&
               right - left + 1U <= span_ &&
               !(forcedUnpaired && (forcedUnpaired->contains(left) || forcedUnpaired->contains(right))) &&
               constraints_[left] != 'x' && constraints_[left] != 'b' &&
               constraints_[right] != 'x' && constraints_[right] != 'b' &&
               canPair(sequence_[left], sequence_[right]);
    }

    [[nodiscard]] auto compute(const std::optional<Interval> forcedUnpaired) const -> double {
        if (forcedUnpaired) {
            for (Index position = forcedUnpaired->begin; position <= forcedUnpaired->end; ++position) {
                if (constraints_[position] == 'p') return logZero;
            }
        }
        const auto length = sequence_.size();
        std::vector<double> partition((length + 1U) * (length + 1U), logZero);
        std::vector<double> pairAny(length * length, logZero);
        std::vector<double> pairInnerStack(length * length, logZero);
        for (Index index{}; index <= length; ++index) partition[intervalOffset(index, index)] = 0.0;

        const auto selectedPair = [&](const Index left, const Index right,
                                      const bool outerStacked) noexcept -> double {
            const auto any = pairAny[matrixOffset(left, right)];
            const auto stacked = pairInnerStack[matrixOffset(left, right)];
            if (options_.noGuHelixEnds && isGuPair(sequence_[left], sequence_[right])) {
                return outerStacked ? stacked : logZero;
            }
            if (options_.noLonelyPairs && !outerStacked) return stacked;
            return any;
        };

        for (Index width = 1U; width <= length; ++width) {
            if (width >= minimumHairpinUnpaired + 2U) {
                for (Index left{}; left + width <= length; ++left) {
                    const auto right = left + width - 1U;
                    if (!pairAllowed(left, right, forcedUnpaired)) continue;
                    const auto inside = partition[intervalOffset(left + 1U, right)];
                    if (inside != logZero) pairAny[matrixOffset(left, right)] = 1.0 + inside;
                    if (left + 1U < right && pairAllowed(left + 1U, right - 1U, forcedUnpaired)) {
                        const auto child = selectedPair(left + 1U, right - 1U, true);
                        if (child != logZero) pairInnerStack[matrixOffset(left, right)] = 1.0 + child;
                    }
                }
            }

            for (Index begin{}; begin + width <= length; ++begin) {
                const auto end = begin + width;
                const auto last = end - 1U;
                auto value = constraints_[last] == 'p'
                    ? logZero : partition[intervalOffset(begin, last)];
                if (forcedUnpaired && forcedUnpaired->contains(last)) {
                    partition[intervalOffset(begin, end)] = value;
                    continue;
                }
                for (Index partner = begin; partner + minimumHairpinUnpaired + 1U <= last; ++partner) {
                    const auto pair = selectedPair(partner, last, false);
                    if (pair == logZero) continue;
                    value = logAdd(value, partition[intervalOffset(begin, partner)] + pair);
                }
                partition[intervalOffset(begin, end)] = value;
            }
        }
        return partition[intervalOffset(0U, length)];
    }

    std::string sequence_;
    FoldingOptions options_;
    std::vector<char> constraints_;
    Index span_{};
    double logPartition_{};
    mutable std::vector<double> cache_;
    mutable std::mutex cacheMutex_;
};

class TurnerFoldingEnsemble final : public FoldingEnsemble {
public:
    TurnerFoldingEnsemble(const Sequence& sequence, FoldingOptions options)
        : sequence_(sequence.str()),
          options_(std::move(options)),
          constraints_(parseConstraint(sequence, options_.constraint)),
          shape_(shapeTerms(sequence, options_)),
          parameters_(folding_detail::loadFoldingParameters(
              options_.temperatureCelsius, options_.parameterSet)),
          rt_(gasConstantKcal * (options_.temperatureCelsius + 273.15)),
          span_(options_.maximumPairSpan == 0U
              ? sequence.size() : std::min(options_.maximumPairSpan, sequence.size())),
          cache_(sequence.size() * (sequence.size() + 1U),
                 std::numeric_limits<double>::quiet_NaN()) {
        if (!std::isfinite(options_.partitionScale) || options_.partitionScale < 1.0) {
            throw std::invalid_argument("partition scale must be finite and at least 1");
        }
        if (options_.maximumInternalLoop > 30U) {
            throw std::invalid_argument("the Turner parameter model supports internal loops up to 30 bases");
        }
        logPartition_ = compute(std::nullopt);
        if (!std::isfinite(logPartition_)) {
            throw std::invalid_argument("folding constraints admit no secondary structure");
        }
    }

    [[nodiscard]] auto logPartition() const noexcept -> double override { return logPartition_; }
    [[nodiscard]] auto ensembleFreeEnergy() const noexcept -> Energy override {
        return -rt_ * logPartition_;
    }

    [[nodiscard]] auto jointUnpairedProbability(const Interval interval) const -> double override {
        if (interval.begin > interval.end || interval.end >= sequence_.size()) {
            throw std::out_of_range("folding interval is out of range");
        }
        for (Index position = interval.begin; position <= interval.end; ++position) {
            if (constraints_[position] == 'p') return 0.0;
        }
        const auto offset = interval.begin * (sequence_.size() + 1U) + interval.size();
        {
            const std::scoped_lock lock(cacheMutex_);
            if (!std::isnan(cache_[offset])) return cache_[offset];
        }
        const auto constrained = compute(interval);
        const auto probability = constrained == logZero
            ? 0.0 : std::clamp(std::exp(std::min(0.0, constrained - logPartition_)), 0.0, 1.0);
        {
            const std::scoped_lock lock(cacheMutex_);
            cache_[offset] = probability;
        }
        return probability;
    }

private:
    [[nodiscard]] auto matrixOffset(const Index left, const Index right) const noexcept -> Index {
        return left * sequence_.size() + right;
    }

    [[nodiscard]] auto intervalOffset(const Index begin, const Index end) const noexcept -> Index {
        return begin * (sequence_.size() + 1U) + end;
    }

    [[nodiscard]] auto forced(const Index position, const std::optional<Interval> interval) const noexcept -> bool {
        return interval && interval->contains(position);
    }

    [[nodiscard]] auto unpairedAllowed(const Index position,
                                       const std::optional<Interval>) const noexcept -> bool {
        return constraints_[position] != 'p';
    }

    [[nodiscard]] auto pairAllowed(const Index left, const Index right,
                                   const std::optional<Interval> interval) const noexcept -> bool {
        return right > left && right - left + 1U <= span_ &&
               !forced(left, interval) && !forced(right, interval) &&
               constraints_[left] != 'x' && constraints_[left] != 'b' &&
               constraints_[right] != 'x' && constraints_[right] != 'b' &&
               canPair(sequence_[left], sequence_[right]);
    }

    [[nodiscard]] auto allUnpairedAllowed(const Index begin, const Index end,
                                          const std::optional<Interval> interval) const noexcept -> bool {
        for (Index position = begin; position < end; ++position) {
            if (!unpairedAllowed(position, interval)) return false;
        }
        return true;
    }

    [[nodiscard]] auto unpairedEnergy(const Index begin, const Index end,
                                      const Energy thermodynamicPerBase = 0.0) const noexcept -> Energy {
        Energy result = thermodynamicPerBase * static_cast<Energy>(end - begin);
        for (Index position = begin; position < end; ++position) result += shape_.unpaired[position];
        return result;
    }

    [[nodiscard]] auto pairPseudoEnergy(const Index left, const Index right) const noexcept -> Energy {
        return shape_.paired[left] + shape_.paired[right];
    }

    [[nodiscard]] auto stackPseudoEnergy(const Index outerLeft, const Index outerRight,
                                         const Index innerLeft, const Index innerRight) const noexcept -> Energy {
        return shape_.stack[outerLeft] + shape_.stack[outerRight] +
               shape_.stack[innerLeft] + shape_.stack[innerRight];
    }

    [[nodiscard]] auto terminalPenalty(const PairType type) const noexcept -> Energy {
        return weakPair(type) ? parameters_->terminalAu : 0.0;
    }

    [[nodiscard]] auto hairpinEnergy(const Index left, const Index right) const noexcept -> Energy {
        const auto length = right - left - 1U;
        const auto motif = sequence_.substr(left, right - left + 1U);
        if (const auto special = parameters_->specialHairpins.find(motif);
            special != parameters_->specialHairpins.end()) {
            return special->second;
        }
        const auto type = pairType(sequence_[left], sequence_[right]);
        auto result = loopInitiation(parameters_->hairpin, length, parameters_->logarithmicLoopSlope);
        if (length == 3U) {
            result += terminalPenalty(type);
        } else {
            result += tableMismatch(parameters_->mismatchHairpin, type,
                                    sequence_[left + 1U], sequence_[right - 1U]);
        }
        return result;
    }

    [[nodiscard]] auto interiorEnergy(const Index outerLeft, const Index outerRight,
                                      const Index innerLeft, const Index innerRight) const noexcept -> Energy {
        const auto leftGap = innerLeft - outerLeft - 1U;
        const auto rightGap = outerRight - innerRight - 1U;
        const auto outer = pairType(sequence_[outerLeft], sequence_[outerRight]);
        const auto inner = reversePair(pairType(sequence_[innerLeft], sequence_[innerRight]));
        if (leftGap == 0U && rightGap == 0U) return tableStack(*parameters_, outer, inner);
        const auto total = leftGap + rightGap;
        if (leftGap == 0U || rightGap == 0U) {
            auto result = loopInitiation(parameters_->bulge, total, parameters_->logarithmicLoopSlope);
            if (total == 1U) {
                result += tableStack(*parameters_, outer, inner);
            } else {
                result += terminalPenalty(outer) + terminalPenalty(inner);
            }
            return result;
        }
        const auto outerLeftBase = sequence_[outerLeft + 1U];
        const auto outerRightBase = sequence_[outerRight - 1U];
        const auto innerLeftBase = sequence_[innerLeft - 1U];
        const auto innerRightBase = sequence_[innerRight + 1U];
        if (leftGap == 1U && rightGap == 1U) {
            return tableInt11(*parameters_, outer, inner, outerLeftBase, outerRightBase);
        }
        if (leftGap == 1U && rightGap == 2U) {
            return tableInt21(*parameters_, outer, inner,
                              outerLeftBase, innerRightBase, outerRightBase);
        }
        if (leftGap == 2U && rightGap == 1U) {
            return tableInt21(*parameters_, inner, outer,
                              innerRightBase, outerLeftBase, innerLeftBase);
        }
        if (leftGap == 2U && rightGap == 2U) {
            if (const auto exact = tableInt22(*parameters_, outer, inner,
                                              outerLeftBase, innerLeftBase,
                                              innerRightBase, outerRightBase)) {
                return *exact;
            }
        }
        auto result = loopInitiation(parameters_->internal, total, parameters_->logarithmicLoopSlope);
        const auto asymmetry = leftGap > rightGap ? leftGap - rightGap : rightGap - leftGap;
        result += std::min(parameters_->ninioMaximum,
                           parameters_->ninioSlope * static_cast<Energy>(asymmetry));
        const std::vector<Energy>* mismatch = &parameters_->mismatchInternal;
        if (leftGap == 1U || rightGap == 1U) mismatch = &parameters_->mismatchInternal1n;
        else if ((leftGap == 2U && rightGap == 3U) ||
                 (leftGap == 3U && rightGap == 2U)) mismatch = &parameters_->mismatchInternal23;
        result += tableMismatch(*mismatch, outer, outerLeftBase, outerRightBase);
        result += tableMismatch(*mismatch, inner, innerRightBase, innerLeftBase);
        return result;
    }

    [[nodiscard]] auto compute(const std::optional<Interval> interval) const -> double {
        const auto length = sequence_.size();
        if (interval) {
            for (Index position = interval->begin; position <= interval->end; ++position) {
                if (constraints_[position] == 'p') return logZero;
            }
        }
        std::vector<double> pairAny(length * length, logZero);
        std::vector<double> pairInnerStack(length * length, logZero);
        const auto intervalDimension = length + 1U;
        std::vector<double> multiZero(intervalDimension * intervalDimension, logZero);
        std::vector<double> multiOne(intervalDimension * intervalDimension, logZero);
        std::vector<double> multiTwo(intervalDimension * intervalDimension, logZero);
        for (Index index{}; index <= length; ++index) {
            multiZero[intervalOffset(index, index)] = 0.0;
        }

        const auto selectedPair = [&](const Index left, const Index right,
                                      const bool outerStacked) -> double {
            const auto type = pairType(sequence_[left], sequence_[right]);
            if (options_.noGuHelixEnds && guPair(type)) {
                // A GU pair is internal to a helix only when it has both the
                // already-known outer neighbour and an immediate inner stack.
                return outerStacked ? pairInnerStack[matrixOffset(left, right)] : logZero;
            }
            if (options_.noLonelyPairs && !outerStacked) {
                return pairInnerStack[matrixOffset(left, right)];
            }
            return pairAny[matrixOffset(left, right)];
        };
        const auto boltzmannLog = [this](const Energy energy) noexcept -> double {
            return std::isfinite(energy) ? -energy / rt_ : logZero;
        };

        for (Index width = 1U; width <= length; ++width) {
            if (width >= minimumHairpinUnpaired + 2U) {
                for (Index left{}; left + width <= length; ++left) {
                    const auto right = left + width - 1U;
                    if (!pairAllowed(left, right, interval)) continue;
                    const auto outerType = pairType(sequence_[left], sequence_[right]);
                    const auto pairSoft = pairPseudoEnergy(left, right);
                    auto any = logZero;
                    auto stackOnly = logZero;

                    if ((!options_.noGuHelixEnds || !guPair(outerType)) &&
                        allUnpairedAllowed(left + 1U, right, interval)) {
                        const auto energy = hairpinEnergy(left, right) + pairSoft +
                                            unpairedEnergy(left + 1U, right);
                        any = logAdd(any, boltzmannLog(energy));
                    }

                    for (Index leftGap{}; leftGap <= options_.maximumInternalLoop; ++leftGap) {
                        const auto innerLeft = left + 1U + leftGap;
                        if (innerLeft >= right) break;
                        for (Index rightGap{}; leftGap + rightGap <= options_.maximumInternalLoop; ++rightGap) {
                            const auto innerRight = right - 1U - rightGap;
                            if (innerLeft >= innerRight) break;
                            if (!allUnpairedAllowed(left + 1U, innerLeft, interval) ||
                                !allUnpairedAllowed(innerRight + 1U, right, interval)) continue;
                            const bool isStack = leftGap == 0U && rightGap == 0U;
                            if (options_.noGuHelixEnds && guPair(outerType) && !isStack) continue;
                            const auto child = selectedPair(innerLeft, innerRight, isStack);
                            if (child == logZero) continue;
                            auto energy = interiorEnergy(left, right, innerLeft, innerRight) + pairSoft +
                                unpairedEnergy(left + 1U, innerLeft) +
                                unpairedEnergy(innerRight + 1U, right);
                            if (isStack) {
                                energy += stackPseudoEnergy(left, right, innerLeft, innerRight);
                            }
                            const auto term = child + boltzmannLog(energy);
                            any = logAdd(any, term);
                            if (isStack) stackOnly = logAdd(stackOnly, term);
                        }
                    }

                    if (!options_.noGuHelixEnds || !guPair(outerType)) {
                        const auto content = multiTwo[intervalOffset(left + 1U, right)];
                        if (content != logZero) {
                            const auto closingType = reversePair(outerType);
                            auto energy = parameters_->multiloopClosing + parameters_->multiloopStem +
                                          terminalPenalty(closingType) + pairSoft;
                            if (options_.includeDangles) {
                                energy += tableMismatch(
                                    parameters_->mismatchMulti, closingType,
                                    sequence_[right - 1U], sequence_[left + 1U]);
                            }
                            any = logAdd(any, content + boltzmannLog(energy));
                        }
                    }
                    pairAny[matrixOffset(left, right)] = any;
                    pairInnerStack[matrixOffset(left, right)] = stackOnly;
                }
            }

            for (Index begin{}; begin + width <= length; ++begin) {
                const auto end = begin + width;
                auto zero = logZero;
                auto one = logZero;
                auto two = logZero;
                if (unpairedAllowed(begin, interval)) {
                    const auto unpaired = boltzmannLog(parameters_->multiloopUnpaired + shape_.unpaired[begin]);
                    zero = multiZero[intervalOffset(begin + 1U, end)] + unpaired;
                    if (multiOne[intervalOffset(begin + 1U, end)] != logZero) {
                        one = multiOne[intervalOffset(begin + 1U, end)] + unpaired;
                    }
                    if (multiTwo[intervalOffset(begin + 1U, end)] != logZero) {
                        two = multiTwo[intervalOffset(begin + 1U, end)] + unpaired;
                    }
                }
                for (Index right = begin + minimumHairpinUnpaired + 1U; right < end; ++right) {
                    const auto pair = selectedPair(begin, right, false);
                    if (pair == logZero) continue;
                    const auto type = pairType(sequence_[begin], sequence_[right]);
                    auto branchEnergy = parameters_->multiloopStem + terminalPenalty(type);
                    if (options_.includeDangles) {
                        // Vienna's dangle-2 model decorates a multiloop stem
                        // from its direct sequence neighbours even when a
                        // neighbour is the endpoint of an adjacent stem.
                        if (begin > 0U && right + 1U < length) {
                            branchEnergy += tableMismatch(
                                parameters_->mismatchMulti, type,
                                sequence_[begin - 1U], sequence_[right + 1U]);
                        } else if (begin > 0U) {
                            branchEnergy += tableDangle(
                                parameters_->dangle5, type, sequence_[begin - 1U]);
                        } else if (right + 1U < length) {
                            branchEnergy += tableDangle(
                                parameters_->dangle3, type, sequence_[right + 1U]);
                        }
                    }
                    const auto branch = pair + boltzmannLog(branchEnergy);
                    const auto suffixZero = multiZero[intervalOffset(right + 1U, end)];
                    if (suffixZero != logZero) one = logAdd(one, branch + suffixZero);
                    const auto suffixAtLeastOne = logAdd(multiOne[intervalOffset(right + 1U, end)],
                                                        multiTwo[intervalOffset(right + 1U, end)]);
                    if (suffixAtLeastOne != logZero) two = logAdd(two, branch + suffixAtLeastOne);
                }
                multiZero[intervalOffset(begin, end)] = zero;
                multiOne[intervalOffset(begin, end)] = one;
                multiTwo[intervalOffset(begin, end)] = two;
            }
        }

        // Exterior structures have a unique left-to-right tokenization. In
        // the dangle-2 model each exterior stem is decorated from its direct
        // sequence neighbours, irrespective of whether a neighbour is
        // unpaired or is the endpoint of an adjacent stem.
        std::vector<double> exterior(length + 1U, logZero);
        exterior[0U] = 0.0;
        for (Index cursor{}; cursor < length; ++cursor) {
            const auto source = exterior[cursor];
            if (source == logZero) continue;
            if (unpairedAllowed(cursor, interval)) {
                exterior[cursor + 1U] = logAdd(
                    exterior[cursor + 1U], source + boltzmannLog(shape_.unpaired[cursor]));
            }
            for (Index right = cursor + minimumHairpinUnpaired + 1U; right < length; ++right) {
                const auto pair = selectedPair(cursor, right, false);
                if (pair == logZero) continue;
                const auto type = pairType(sequence_[cursor], sequence_[right]);
                auto branchEnergy = terminalPenalty(type);
                if (options_.includeDangles) {
                    const bool hasLeft = cursor > 0U;
                    const bool hasRight = right + 1U < length;
                    if (hasLeft && hasRight) {
                        branchEnergy += tableMismatch(
                            parameters_->mismatchExterior, type,
                            sequence_[cursor - 1U], sequence_[right + 1U]);
                    } else if (hasLeft) {
                        branchEnergy += tableDangle(
                            parameters_->dangle5, type, sequence_[cursor - 1U]);
                    } else if (hasRight) {
                        branchEnergy += tableDangle(
                            parameters_->dangle3, type, sequence_[right + 1U]);
                    }
                }
                exterior[right + 1U] = logAdd(
                    exterior[right + 1U], source + pair + boltzmannLog(branchEnergy));
            }
        }
        return exterior[length];
    }

    std::string sequence_;
    FoldingOptions options_;
    std::vector<char> constraints_;
    ShapeTerms shape_;
    std::shared_ptr<const FoldingParameters> parameters_;
    double rt_{};
    Index span_{};
    double logPartition_{};
    mutable std::vector<double> cache_;
    mutable std::mutex cacheMutex_;
};

} // namespace

void validateShapeEncoding(const std::string_view method, const std::string_view conversion) {
    static_cast<void>(parseShapeEncoding(method, conversion));
}

auto makeTurnerFoldingEnsemble(const Sequence& sequence, const FoldingOptions& options)
    -> std::unique_ptr<FoldingEnsemble> {
    return std::make_unique<TurnerFoldingEnsemble>(sequence, options);
}

auto makeBasePairFoldingEnsemble(const Sequence& sequence, const FoldingOptions& options)
    -> std::unique_ptr<FoldingEnsemble> {
    return std::make_unique<BasePairFoldingEnsemble>(sequence, options);
}

} // namespace intarnanew
