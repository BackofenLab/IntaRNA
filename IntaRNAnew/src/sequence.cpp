#include "intarnanew/sequence.hpp"

#include "intarnanew/compression.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace intarnanew {
namespace {

[[nodiscard]] auto normalized(std::string sequence) -> std::expected<std::string, std::string> {
    std::string result;
    result.reserve(sequence.size());
    for (const unsigned char raw : sequence) {
        if (std::isspace(raw) != 0) {
            continue;
        }
        auto symbol = static_cast<char>(std::toupper(raw));
        if (symbol == 'T') {
            symbol = 'U';
        }
        if (nucleotideMask(symbol) == 0U) {
            return std::unexpected("invalid IUPAC nucleotide '" + std::string(1, symbol) + "'");
        }
        if (symbol != 'A' && symbol != 'C' && symbol != 'G' && symbol != 'U') {
            symbol = 'N';
        }
        result.push_back(symbol);
    }
    if (result.empty()) {
        return std::unexpected("RNA sequence is empty");
    }
    return result;
}

[[nodiscard]] auto trimmed(std::string value) -> std::string {
    const auto first = std::find_if_not(value.begin(), value.end(), [](const unsigned char c) {
        return std::isspace(c) != 0;
    });
    const auto last = std::find_if_not(value.rbegin(), value.rend(), [](const unsigned char c) {
        return std::isspace(c) != 0;
    }).base();
    return first < last ? std::string(first, last) : std::string{};
}

[[nodiscard]] auto complementMask(const std::uint8_t mask, const bool allowGu) noexcept -> std::uint8_t {
    // A=1, C=2, G=4, U=8
    std::uint8_t result{};
    if ((mask & 0x1U) != 0U) {
        result |= 0x8U;
    }
    if ((mask & 0x2U) != 0U) {
        result |= 0x4U;
    }
    if ((mask & 0x4U) != 0U) {
        result |= static_cast<std::uint8_t>(0x2U | (allowGu ? 0x8U : 0U));
    }
    if ((mask & 0x8U) != 0U) {
        result |= static_cast<std::uint8_t>(0x1U | (allowGu ? 0x4U : 0U));
    }
    return result;
}

[[nodiscard]] auto coordinateSpanFits(
    const long long firstPosition,
    const std::size_t length) noexcept -> bool {
    if (length == 0U) return false;
    if (!std::in_range<unsigned long long>(length - 1U)) return false;
    const auto lastOffset = static_cast<unsigned long long>(length - 1U);
    constexpr auto signedMaximum = static_cast<unsigned long long>(
        std::numeric_limits<long long>::max());
    if (firstPosition >= 0) {
        return lastOffset <= signedMaximum - static_cast<unsigned long long>(firstPosition);
    }

    // Unsigned negation obtains |LLONG_MIN| without overflowing signed arithmetic.
    const auto negativePositions = 0ULL - static_cast<unsigned long long>(firstPosition);
    if (lastOffset < negativePositions) return true;
    return lastOffset - negativePositions <= signedMaximum - 1U;
}

[[nodiscard]] auto readBounded(
    std::istream& input,
    const std::size_t maximum,
    const std::string_view description) -> std::expected<std::string, std::string> {
    std::string bytes;
    std::array<char, 65536U> buffer{};
    while (input) {
        input.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
        const auto count = input.gcount();
        if (count < 0) return std::unexpected("failed to read " + std::string(description));
        const auto amount = static_cast<std::size_t>(count);
        if (amount > maximum - bytes.size()) {
            return std::unexpected(std::string(description) +
                                   " exceeds the configured input-byte limit");
        }
        bytes.append(buffer.data(), amount);
    }
    if (!input.eof()) return std::unexpected("failed to read " + std::string(description));
    return bytes;
}

[[nodiscard]] auto hasGzipSuffix(const std::string_view value) noexcept -> bool {
    if (value.size() < 3U) return false;
    const auto offset = value.size() - 3U;
    return value[offset] == '.' &&
        static_cast<char>(std::tolower(static_cast<unsigned char>(value[offset + 1U]))) == 'g' &&
        static_cast<char>(std::tolower(static_cast<unsigned char>(value[offset + 2U]))) == 'z';
}

[[nodiscard]] auto parsePossiblyCompressed(
    std::string bytes,
    const bool gzipExpected,
    const std::string_view description,
    const std::string_view fallbackId,
    const long long firstPosition) -> std::expected<std::vector<Sequence>, std::string> {
    if (hasGzipMagic(bytes)) {
        auto decoded = gzipDecompress(bytes);
        if (!decoded) {
            return std::unexpected("cannot decode " + std::string(description) + ": " +
                                   decoded.error());
        }
        bytes = std::move(*decoded);
    } else if (gzipExpected) {
        return std::unexpected(std::string(description) +
                               " has a .gz suffix but no gzip signature");
    }
    std::istringstream input(std::move(bytes));
    return SequenceReader::parseFasta(input, fallbackId, firstPosition);
}

} // namespace

Sequence::Sequence(std::string identifier, std::string nucleotides, const long long firstPosition)
    : identifier_(std::move(identifier)), firstPosition_(firstPosition) {
    auto parsed = normalized(std::move(nucleotides));
    if (!parsed) {
        throw std::invalid_argument(parsed.error());
    }
    nucleotides_ = std::move(*parsed);
    if (!coordinateSpanFits(firstPosition_, nucleotides_.size())) {
        throw std::out_of_range("sequence coordinates exceed the signed 64-bit range");
    }
    if (identifier_.empty()) {
        identifier_ = "sequence";
    }
}

auto Sequence::externalIndex(const Index index) const -> long long {
    if (index >= size()) throw std::out_of_range("sequence index is outside the sequence");
    if (firstPosition_ >= 0) {
        return firstPosition_ + static_cast<long long>(index);
    }
    const auto offset = static_cast<unsigned long long>(index);
    const auto negativePositions = 0ULL - static_cast<unsigned long long>(firstPosition_);
    if (offset < negativePositions) {
        return firstPosition_ + static_cast<long long>(offset);
    }
    return static_cast<long long>(offset - negativePositions + 1U);
}

auto Sequence::internalIndex(const long long external) const -> std::expected<Index, std::string> {
    if (firstPosition_ < 0 && external == 0) {
        return std::unexpected("sequence index zero is not used for a negative coordinate origin");
    }
    const auto adjusted = firstPosition_ < 0 && external > 0 ? external - 1LL : external;
    if (adjusted < firstPosition_) {
        return std::unexpected("sequence index is before the first configured position");
    }
    const auto distance = static_cast<unsigned long long>(adjusted) -
                          static_cast<unsigned long long>(firstPosition_);
    if (distance >= size()) {
        return std::unexpected("sequence index is beyond the sequence end");
    }
    return static_cast<Index>(distance);
}

auto SequenceReader::read(
    const std::string_view specification,
    const std::string_view fallbackId,
    const long long firstPosition,
    std::istream& standardInput) -> std::expected<std::vector<Sequence>, std::string> {
    if (specification.empty()) {
        return std::unexpected("missing RNA sequence input");
    }
    if (specification == "STDIN" || specification == "-") {
        GzipLimits limits;
        auto bytes = readBounded(standardInput, limits.maxCompressedBytes, "standard input");
        if (!bytes) return std::unexpected(bytes.error());
        return parsePossiblyCompressed(
            std::move(*bytes), false, "standard input", fallbackId, firstPosition);
    }

    std::error_code error;
    const std::filesystem::path path{specification};
    if (std::filesystem::is_regular_file(path, error)) {
        std::ifstream input(path, std::ios::binary);
        if (!input) {
            return std::unexpected("cannot open sequence file '" + path.string() + "'");
        }
        GzipLimits limits;
        auto bytes = readBounded(
            input, limits.maxCompressedBytes, "sequence file '" + path.string() + "'");
        if (!bytes) return std::unexpected(bytes.error());
        return parsePossiblyCompressed(
            std::move(*bytes), hasGzipSuffix(path.string()),
            "sequence file '" + path.string() + "'", fallbackId, firstPosition);
    }

    auto parsed = normalized(std::string(specification));
    if (!parsed) {
        return std::unexpected(
            "input is neither a readable file nor an RNA sequence: " + parsed.error());
    }
    try {
        std::vector<Sequence> result;
        result.emplace_back(std::string(fallbackId), std::move(*parsed), firstPosition);
        return result;
    } catch (const std::exception& exception) {
        return std::unexpected(exception.what());
    }
}

auto SequenceReader::parseFasta(
    std::istream& input,
    const std::string_view fallbackId,
    const long long firstPosition) -> std::expected<std::vector<Sequence>, std::string> {
    std::vector<Sequence> result;
    std::string identifier{fallbackId};
    std::string bases;
    std::string line;
    bool sawHeader = false;
    bool openHeader = false;

    const auto commit = [&](const bool requireSequence) -> std::expected<void, std::string> {
        if (bases.empty()) {
            if (requireSequence) {
                return std::unexpected("FASTA record '" + identifier + "' has no RNA sequence");
            }
            return {};
        }
        auto parsed = normalized(std::move(bases));
        bases.clear();
        if (!parsed) {
            return std::unexpected(parsed.error());
        }
        try {
            result.emplace_back(identifier, std::move(*parsed), firstPosition);
        } catch (const std::exception& exception) {
            return std::unexpected(exception.what());
        }
        return {};
    };

    while (std::getline(input, line)) {
        line = trimmed(std::move(line));
        if (line.empty() || line.front() == ';') {
            continue;
        }
        if (line.front() == '>') {
            if (auto status = commit(openHeader); !status) {
                return std::unexpected(status.error());
            }
            identifier = trimmed(line.substr(1));
            if (const auto blank = identifier.find_first_of(" \t"); blank != std::string::npos) {
                identifier.erase(blank);
            }
            if (identifier.empty()) {
                identifier = std::string(fallbackId) + std::to_string(result.size() + 1U);
            }
            sawHeader = true;
            openHeader = true;
        } else {
            bases += line;
        }
    }
    if (auto status = commit(openHeader); !status) {
        return std::unexpected(status.error());
    }
    if (result.empty()) {
        return std::unexpected("no RNA sequence found in input");
    }
    if (!sawHeader && result.size() == 1U && result.front().id().empty()) {
        return std::unexpected("unable to assign sequence identifier");
    }
    return result;
}

auto SequenceReader::select(
    std::vector<Sequence> sequences,
    const std::string_view specification) -> std::expected<std::vector<Sequence>, std::string> {
    if (specification.empty()) {
        return sequences;
    }

    std::vector<bool> selected(sequences.size(), false);
    const auto parseRecordIndex = [&](const std::string_view text)
        -> std::expected<std::size_t, std::string> {
        std::size_t value{};
        const auto [position, error] = std::from_chars(text.data(), text.data() + text.size(), value);
        if (text.empty() || error != std::errc{} || position != text.data() + text.size() ||
            value == 0U || value > sequences.size()) {
            return std::unexpected(
                "sequence-set index '" + std::string(text) + "' is outside 1.." +
                std::to_string(sequences.size()));
        }
        return value;
    };

    std::size_t cursor{};
    while (cursor <= specification.size()) {
        const auto comma = specification.find(',', cursor);
        const auto rawToken = specification.substr(
            cursor, comma == std::string_view::npos ? std::string_view::npos : comma - cursor);
        const auto tokenStorage = trimmed(std::string(rawToken));
        const std::string_view token{tokenStorage};
        if (token.empty()) {
            return std::unexpected("sequence-set specification contains an empty item");
        }

        const auto dash = token.find('-');
        if (dash == std::string_view::npos) {
            auto index = parseRecordIndex(token);
            if (!index) return std::unexpected(index.error());
            selected[*index - 1U] = true;
        } else {
            if (token.find('-', dash + 1U) != std::string_view::npos) {
                return std::unexpected("sequence-set range must use FROM-TO encoding");
            }
            auto first = parseRecordIndex(token.substr(0U, dash));
            auto last = parseRecordIndex(token.substr(dash + 1U));
            if (!first) return std::unexpected(first.error());
            if (!last) return std::unexpected(last.error());
            if (*first > *last) {
                return std::unexpected("sequence-set range start exceeds its end");
            }
            std::fill(selected.begin() + static_cast<std::ptrdiff_t>(*first - 1U),
                      selected.begin() + static_cast<std::ptrdiff_t>(*last), true);
        }

        if (comma == std::string_view::npos) break;
        cursor = comma + 1U;
    }

    std::vector<Sequence> result;
    result.reserve(sequences.size());
    for (std::size_t index = 0U; index < sequences.size(); ++index) {
        if (selected[index]) result.push_back(std::move(sequences[index]));
    }
    return result;
}

auto nucleotideMask(const char symbol) noexcept -> std::uint8_t {
    switch (static_cast<char>(std::toupper(static_cast<unsigned char>(symbol)))) {
        case 'A': return 0x1U;
        case 'C': return 0x2U;
        case 'G': return 0x4U;
        case 'U': case 'T': return 0x8U;
        case 'R': return 0x5U;
        case 'Y': return 0xAU;
        case 'S': return 0x6U;
        case 'W': return 0x9U;
        case 'K': return 0xCU;
        case 'M': return 0x3U;
        case 'B': return 0xEU;
        case 'D': return 0xDU;
        case 'H': return 0xBU;
        case 'V': return 0x7U;
        case 'N': return 0xFU;
        default: return 0U;
    }
}

auto canPair(const char target, const char query, const bool allowGu) noexcept -> bool {
    const auto canonicalMask = [](const char symbol) noexcept -> std::uint8_t {
        switch (static_cast<char>(std::toupper(static_cast<unsigned char>(symbol)))) {
            case 'A': return 0x1U;
            case 'C': return 0x2U;
            case 'G': return 0x4U;
            case 'U': case 'T': return 0x8U;
            default: return 0U;
        }
    };
    const auto targetOptions = canonicalMask(target);
    const auto queryOptions = canonicalMask(query);
    if (targetOptions == 0U || queryOptions == 0U) return false;
    return (complementMask(targetOptions, allowGu) & queryOptions) != 0U;
}

auto isGuPair(const char target, const char query) noexcept -> bool {
    const auto left = static_cast<char>(std::toupper(static_cast<unsigned char>(target)));
    const auto right = static_cast<char>(std::toupper(static_cast<unsigned char>(query)));
    return (left == 'G' && (right == 'U' || right == 'T')) ||
           ((left == 'U' || left == 'T') && right == 'G');
}

auto reverseComplement(const std::string_view sequence) -> std::string {
    std::string result;
    result.reserve(sequence.size());
    for (auto iterator = sequence.rbegin(); iterator != sequence.rend(); ++iterator) {
        switch (*iterator) {
            case 'A': result.push_back('U'); break;
            case 'C': result.push_back('G'); break;
            case 'G': result.push_back('C'); break;
            case 'U': case 'T': result.push_back('A'); break;
            default: result.push_back('N'); break;
        }
    }
    return result;
}

} // namespace intarnanew
