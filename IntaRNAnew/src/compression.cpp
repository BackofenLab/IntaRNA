#include "intarnanew/compression.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace intarnanew {
namespace {

constexpr std::uint8_t gzipId1{0x1fU};
constexpr std::uint8_t gzipId2{0x8bU};
constexpr std::uint8_t deflateMethod{8U};

[[nodiscard]] constexpr auto makeCrcTable() noexcept -> std::array<std::uint32_t, 256U> {
    std::array<std::uint32_t, 256U> table{};
    for (std::uint32_t value = 0U; value < table.size(); ++value) {
        auto remainder = value;
        for (unsigned bit = 0U; bit < 8U; ++bit) {
            remainder = (remainder >> 1U) ^
                ((remainder & 1U) != 0U ? 0xedb88320U : 0U);
        }
        table[value] = remainder;
    }
    return table;
}

constexpr auto crcTable = makeCrcTable();

[[nodiscard]] auto crc32Range(
    const std::string_view bytes,
    std::uint32_t state = 0xffffffffU) noexcept -> std::uint32_t {
    for (const unsigned char byte : bytes) {
        state = crcTable[(state ^ byte) & 0xffU] ^ (state >> 8U);
    }
    return state;
}

[[nodiscard]] auto byteAt(const std::string_view bytes, const std::size_t position) noexcept
    -> std::uint8_t {
    return static_cast<std::uint8_t>(static_cast<unsigned char>(bytes[position]));
}

[[nodiscard]] auto little16(const std::string_view bytes, const std::size_t position) noexcept
    -> std::uint16_t {
    return static_cast<std::uint16_t>(byteAt(bytes, position)) |
           static_cast<std::uint16_t>(static_cast<std::uint16_t>(byteAt(bytes, position + 1U)) << 8U);
}

[[nodiscard]] auto little32(const std::string_view bytes, const std::size_t position) noexcept
    -> std::uint32_t {
    return static_cast<std::uint32_t>(byteAt(bytes, position)) |
           (static_cast<std::uint32_t>(byteAt(bytes, position + 1U)) << 8U) |
           (static_cast<std::uint32_t>(byteAt(bytes, position + 2U)) << 16U) |
           (static_cast<std::uint32_t>(byteAt(bytes, position + 3U)) << 24U);
}

void append16(std::string& output, const std::uint16_t value) {
    output.push_back(static_cast<char>(value & 0xffU));
    output.push_back(static_cast<char>((value >> 8U) & 0xffU));
}

void append32(std::string& output, const std::uint32_t value) {
    output.push_back(static_cast<char>(value & 0xffU));
    output.push_back(static_cast<char>((value >> 8U) & 0xffU));
    output.push_back(static_cast<char>((value >> 16U) & 0xffU));
    output.push_back(static_cast<char>((value >> 24U) & 0xffU));
}

class BitReader {
public:
    BitReader(const std::string_view bytes, const std::size_t start)
        : bytes_(bytes), bytePosition_(start) {}

    [[nodiscard]] auto read(const unsigned count) -> std::optional<std::uint32_t> {
        if (count > 24U) return std::nullopt;
        while (bufferedBits_ < count) {
            if (bytePosition_ >= bytes_.size()) return std::nullopt;
            bitBuffer_ |= static_cast<std::uint64_t>(byteAt(bytes_, bytePosition_)) << bufferedBits_;
            bufferedBits_ += 8U;
            ++bytePosition_;
        }
        const auto mask = count == 0U ? 0U : (std::uint64_t{1U} << count) - 1U;
        const auto result = static_cast<std::uint32_t>(bitBuffer_ & mask);
        bitBuffer_ >>= count;
        bufferedBits_ -= count;
        return result;
    }

    void alignToByte() noexcept {
        bitBuffer_ = 0U;
        bufferedBits_ = 0U;
    }

    [[nodiscard]] auto bytePosition() const noexcept -> std::size_t {
        return bytePosition_;
    }

    [[nodiscard]] auto bitPosition() const noexcept -> unsigned {
        return bufferedBits_ == 0U ? 0U : 8U - bufferedBits_;
    }

private:
    std::string_view bytes_;
    std::size_t bytePosition_{};
    std::uint64_t bitBuffer_{};
    unsigned bufferedBits_{};
};

class HuffmanTree {
public:
    [[nodiscard]] static auto build(
        const std::span<const std::uint8_t> lengths,
        const unsigned permittedMaximum,
        const bool permitEmpty,
        const bool permitSingleSymbol) -> std::expected<HuffmanTree, std::string> {
        HuffmanTree result;
        result.counts_.assign(permittedMaximum + 1U, 0U);
        for (const auto length : lengths) {
            if (length > permittedMaximum) {
                return std::unexpected("DEFLATE Huffman code length exceeds its limit");
            }
            if (length != 0U) ++result.counts_[length];
        }

        unsigned symbolCount{};
        for (unsigned length = 1U; length <= permittedMaximum; ++length) {
            symbolCount += result.counts_[length];
        }
        if (symbolCount == 0U) {
            if (permitEmpty) return result;
            return std::unexpected("DEFLATE Huffman alphabet has no symbols");
        }

        int available{1};
        for (unsigned length = 1U; length <= permittedMaximum; ++length) {
            available = available * 2 - static_cast<int>(result.counts_[length]);
            if (available < 0) {
                return std::unexpected("DEFLATE Huffman alphabet is oversubscribed");
            }
        }
        if (available != 0 &&
            !(permitSingleSymbol && symbolCount == 1U && result.counts_[1U] == 1U)) {
            return std::unexpected("DEFLATE Huffman alphabet is incomplete");
        }

        result.firstCodes_.assign(permittedMaximum + 1U, 0U);
        result.firstOffsets_.assign(permittedMaximum + 1U, 0U);
        unsigned code{};
        unsigned offset{};
        for (unsigned length = 1U; length <= permittedMaximum; ++length) {
            code = (code + result.counts_[length - 1U]) << 1U;
            result.firstCodes_[length] = code;
            result.firstOffsets_[length] = offset;
            offset += result.counts_[length];
        }
        result.symbols_.resize(symbolCount);
        auto nextOffset = result.firstOffsets_;
        for (std::size_t symbol = 0U; symbol < lengths.size(); ++symbol) {
            const auto length = lengths[symbol];
            if (length == 0U) continue;
            result.symbols_[nextOffset[length]++] = static_cast<unsigned>(symbol);
            result.maximumLength_ = std::max(result.maximumLength_, static_cast<unsigned>(length));
        }
        return result;
    }

    [[nodiscard]] auto decode(BitReader& input) const -> std::expected<unsigned, std::string> {
        if (maximumLength_ == 0U) {
            return std::unexpected("DEFLATE distance alphabet is empty but a match was requested");
        }
        unsigned code{};
        for (unsigned length = 1U; length <= maximumLength_; ++length) {
            const auto bit = input.read(1U);
            if (!bit) return std::unexpected("truncated DEFLATE Huffman code");
            code = (code << 1U) | *bit;
            if (code >= firstCodes_[length]) {
                const auto relative = code - firstCodes_[length];
                if (relative < counts_[length]) {
                    return symbols_[firstOffsets_[length] + relative];
                }
            }
        }
        return std::unexpected("invalid DEFLATE Huffman code");
    }

private:
    std::vector<unsigned> counts_;
    std::vector<unsigned> firstCodes_;
    std::vector<unsigned> firstOffsets_;
    std::vector<unsigned> symbols_;
    unsigned maximumLength_{};
};

struct DeflateTrees {
    HuffmanTree literals;
    HuffmanTree distances;
};

[[nodiscard]] auto fixedTrees() -> std::expected<DeflateTrees, std::string> {
    std::array<std::uint8_t, 288U> literalLengths{};
    std::fill(literalLengths.begin(), literalLengths.begin() + 144, 8U);
    std::fill(literalLengths.begin() + 144, literalLengths.begin() + 256, 9U);
    std::fill(literalLengths.begin() + 256, literalLengths.begin() + 280, 7U);
    std::fill(literalLengths.begin() + 280, literalLengths.end(), 8U);
    std::array<std::uint8_t, 32U> distanceLengths{};
    distanceLengths.fill(5U);
    auto literals = HuffmanTree::build(literalLengths, 15U, false, false);
    if (!literals) return std::unexpected(literals.error());
    auto distances = HuffmanTree::build(distanceLengths, 15U, false, false);
    if (!distances) return std::unexpected(distances.error());
    return DeflateTrees{std::move(*literals), std::move(*distances)};
}

[[nodiscard]] auto dynamicTrees(BitReader& input) -> std::expected<DeflateTrees, std::string> {
    const auto rawLiteralCount = input.read(5U);
    const auto rawDistanceCount = input.read(5U);
    const auto rawCodeLengthCount = input.read(4U);
    if (!rawLiteralCount || !rawDistanceCount || !rawCodeLengthCount) {
        return std::unexpected("truncated dynamic DEFLATE header");
    }
    const auto literalCount = static_cast<std::size_t>(*rawLiteralCount) + 257U;
    const auto distanceCount = static_cast<std::size_t>(*rawDistanceCount) + 1U;
    const auto codeLengthCount = static_cast<std::size_t>(*rawCodeLengthCount) + 4U;
    if (*rawLiteralCount > 29U) {
        return std::unexpected("dynamic DEFLATE header declares reserved literal/length symbols");
    }

    constexpr std::array<unsigned, 19U> order{
        16U, 17U, 18U, 0U, 8U, 7U, 9U, 6U, 10U, 5U,
        11U, 4U, 12U, 3U, 13U, 2U, 14U, 1U, 15U,
    };
    std::array<std::uint8_t, 19U> codeLengths{};
    for (std::size_t index = 0U; index < codeLengthCount; ++index) {
        const auto length = input.read(3U);
        if (!length) return std::unexpected("truncated dynamic DEFLATE code-length alphabet");
        codeLengths[order[index]] = static_cast<std::uint8_t>(*length);
    }
    auto codeLengthTree = HuffmanTree::build(codeLengths, 7U, false, false);
    if (!codeLengthTree) {
        return std::unexpected("invalid dynamic DEFLATE code-length alphabet: " +
                               codeLengthTree.error());
    }

    std::vector<std::uint8_t> lengths;
    lengths.reserve(literalCount + distanceCount);
    while (lengths.size() < literalCount + distanceCount) {
        auto symbol = codeLengthTree->decode(input);
        if (!symbol) return std::unexpected(symbol.error());
        if (*symbol <= 15U) {
            lengths.push_back(static_cast<std::uint8_t>(*symbol));
            continue;
        }

        std::size_t repeat{};
        std::uint8_t value{};
        if (*symbol == 16U) {
            if (lengths.empty()) {
                return std::unexpected("dynamic DEFLATE repeat has no preceding code length");
            }
            const auto extra = input.read(2U);
            if (!extra) return std::unexpected("truncated dynamic DEFLATE repeat");
            repeat = static_cast<std::size_t>(*extra) + 3U;
            value = lengths.back();
        } else if (*symbol == 17U) {
            const auto extra = input.read(3U);
            if (!extra) return std::unexpected("truncated dynamic DEFLATE zero repeat");
            repeat = static_cast<std::size_t>(*extra) + 3U;
        } else if (*symbol == 18U) {
            const auto extra = input.read(7U);
            if (!extra) return std::unexpected("truncated dynamic DEFLATE long zero repeat");
            repeat = static_cast<std::size_t>(*extra) + 11U;
        } else {
            return std::unexpected("invalid dynamic DEFLATE code-length symbol");
        }
        if (repeat > literalCount + distanceCount - lengths.size()) {
            return std::unexpected("dynamic DEFLATE code-length repeat exceeds its alphabets");
        }
        lengths.insert(lengths.end(), repeat, value);
    }

    if (lengths[256U] == 0U) {
        return std::unexpected("dynamic DEFLATE literal alphabet omits end-of-block");
    }
    auto literals = HuffmanTree::build(
        std::span<const std::uint8_t>{lengths.data(), literalCount}, 15U, false, true);
    if (!literals) {
        return std::unexpected("invalid dynamic DEFLATE literal alphabet: " + literals.error());
    }
    auto distances = HuffmanTree::build(
        std::span<const std::uint8_t>{lengths.data() + literalCount, distanceCount},
        15U, true, true);
    if (!distances) {
        return std::unexpected("invalid dynamic DEFLATE distance alphabet: " + distances.error());
    }
    return DeflateTrees{std::move(*literals), std::move(*distances)};
}

constexpr std::array<unsigned, 29U> lengthBase{
    3U, 4U, 5U, 6U, 7U, 8U, 9U, 10U,
    11U, 13U, 15U, 17U, 19U, 23U, 27U, 31U,
    35U, 43U, 51U, 59U, 67U, 83U, 99U, 115U,
    131U, 163U, 195U, 227U, 258U,
};
constexpr std::array<unsigned, 29U> lengthExtra{
    0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
    1U, 1U, 1U, 1U, 2U, 2U, 2U, 2U,
    3U, 3U, 3U, 3U, 4U, 4U, 4U, 4U,
    5U, 5U, 5U, 5U, 0U,
};
constexpr std::array<unsigned, 30U> distanceBase{
    1U, 2U, 3U, 4U, 5U, 7U, 9U, 13U,
    17U, 25U, 33U, 49U, 65U, 97U, 129U, 193U,
    257U, 385U, 513U, 769U, 1025U, 1537U, 2049U, 3073U,
    4097U, 6145U, 8193U, 12289U, 16385U, 24577U,
};
constexpr std::array<unsigned, 30U> distanceExtra{
    0U, 0U, 0U, 0U, 1U, 1U, 2U, 2U,
    3U, 3U, 4U, 4U, 5U, 5U, 6U, 6U,
    7U, 7U, 8U, 8U, 9U, 9U, 10U, 10U,
    11U, 11U, 12U, 12U, 13U, 13U,
};

[[nodiscard]] auto appendDecodedByte(
    std::string& output,
    const char value,
    const GzipLimits& limits) -> std::expected<void, std::string> {
    if (output.size() >= limits.maxDecompressedBytes) {
        return std::unexpected("gzip output exceeds the configured decompressed-byte limit");
    }
    output.push_back(value);
    return {};
}

[[nodiscard]] auto decodeCompressedBlock(
    BitReader& input,
    const DeflateTrees& trees,
    std::string& output,
    const std::size_t memberStart,
    const GzipLimits& limits) -> std::expected<void, std::string> {
    while (true) {
        auto literal = trees.literals.decode(input);
        if (!literal) return std::unexpected(literal.error());
        if (*literal < 256U) {
            auto status = appendDecodedByte(output, static_cast<char>(*literal), limits);
            if (!status) return status;
            continue;
        }
        if (*literal == 256U) return {};
        if (*literal < 257U || *literal > 285U) {
            return std::unexpected("invalid DEFLATE length symbol");
        }
        const auto lengthIndex = *literal - 257U;
        auto length = lengthBase[lengthIndex];
        if (const auto extraCount = lengthExtra[lengthIndex]; extraCount != 0U) {
            const auto extra = input.read(extraCount);
            if (!extra) return std::unexpected("truncated DEFLATE match length");
            length += *extra;
        }

        auto distanceSymbol = trees.distances.decode(input);
        if (!distanceSymbol) return std::unexpected(distanceSymbol.error());
        if (*distanceSymbol >= distanceBase.size()) {
            return std::unexpected("invalid DEFLATE distance symbol");
        }
        auto distance = distanceBase[*distanceSymbol];
        if (const auto extraCount = distanceExtra[*distanceSymbol]; extraCount != 0U) {
            const auto extra = input.read(extraCount);
            if (!extra) return std::unexpected("truncated DEFLATE match distance");
            distance += *extra;
        }
        const auto producedInMember = output.size() - memberStart;
        if (distance == 0U || distance > producedInMember) {
            return std::unexpected("DEFLATE match distance precedes the current member");
        }
        if (length > limits.maxDecompressedBytes - output.size()) {
            return std::unexpected("gzip output exceeds the configured decompressed-byte limit");
        }
        for (unsigned index = 0U; index < length; ++index) {
            output.push_back(output[output.size() - distance]);
        }
    }
}

[[nodiscard]] auto decodeDeflate(
    const std::string_view bytes,
    const std::size_t start,
    std::string& output,
    const std::size_t memberStart,
    const GzipLimits& limits) -> std::expected<std::size_t, std::string> {
    BitReader input(bytes, start);
    bool finalBlock{};
    std::size_t blockCount{};
    while (!finalBlock) {
        if (++blockCount > limits.maxDeflateBlocks) {
            return std::unexpected("gzip member exceeds the configured DEFLATE block limit");
        }
        const auto final = input.read(1U);
        const auto type = input.read(2U);
        if (!final || !type) return std::unexpected("truncated DEFLATE block header");
        finalBlock = *final != 0U;
        if (*type == 0U) {
            input.alignToByte();
            const auto position = input.bytePosition();
            if (position > bytes.size() || bytes.size() - position < 4U) {
                return std::unexpected("truncated stored DEFLATE block header");
            }
            const auto length = little16(bytes, position);
            const auto inverse = little16(bytes, position + 2U);
            if (static_cast<std::uint16_t>(length ^ 0xffffU) != inverse) {
                return std::unexpected("stored DEFLATE block has inconsistent LEN/NLEN");
            }
            const auto payload = position + 4U;
            if (payload > bytes.size() || bytes.size() - payload < length) {
                return std::unexpected("truncated stored DEFLATE block payload");
            }
            if (static_cast<std::size_t>(length) > limits.maxDecompressedBytes - output.size()) {
                return std::unexpected("gzip output exceeds the configured decompressed-byte limit");
            }
            output.append(bytes.substr(payload, length));
            input = BitReader(bytes, payload + length);
        } else if (*type == 1U) {
            static const auto fixed = fixedTrees();
            if (!fixed) return std::unexpected(fixed.error());
            auto status = decodeCompressedBlock(
                input, *fixed, output, memberStart, limits);
            if (!status) return std::unexpected(status.error());
        } else if (*type == 2U) {
            auto trees = dynamicTrees(input);
            if (!trees) return std::unexpected(trees.error());
            auto status = decodeCompressedBlock(
                input, *trees, output, memberStart, limits);
            if (!status) return std::unexpected(status.error());
        } else {
            return std::unexpected("DEFLATE block uses the reserved block type");
        }
    }
    if (input.bitPosition() != 0U) input.alignToByte();
    return input.bytePosition();
}

[[nodiscard]] auto advanceZeroTerminated(
    const std::string_view bytes,
    std::size_t& position,
    const std::size_t headerStart,
    const GzipLimits& limits,
    const std::string_view fieldName) -> std::expected<void, std::string> {
    while (position < bytes.size() && byteAt(bytes, position) != 0U) {
        ++position;
        if (position - headerStart > limits.maxHeaderBytes) {
            return std::unexpected("gzip header exceeds the configured header-byte limit");
        }
    }
    if (position >= bytes.size()) {
        return std::unexpected("truncated gzip " + std::string(fieldName));
    }
    ++position;
    return {};
}

[[nodiscard]] auto parseHeader(
    const std::string_view bytes,
    const std::size_t memberStart,
    const GzipLimits& limits) -> std::expected<std::size_t, std::string> {
    if (memberStart > bytes.size() || bytes.size() - memberStart < 10U) {
        return std::unexpected("truncated gzip member header");
    }
    if (byteAt(bytes, memberStart) != gzipId1 || byteAt(bytes, memberStart + 1U) != gzipId2) {
        return std::unexpected("gzip member is missing the 1f 8b signature");
    }
    if (byteAt(bytes, memberStart + 2U) != deflateMethod) {
        return std::unexpected("gzip member uses an unsupported compression method");
    }
    const auto flags = byteAt(bytes, memberStart + 3U);
    if ((flags & 0xe0U) != 0U) {
        return std::unexpected("gzip member sets reserved header flags");
    }
    std::size_t position = memberStart + 10U;
    if ((flags & 0x04U) != 0U) {
        if (bytes.size() - position < 2U) {
            return std::unexpected("truncated gzip extra-field length");
        }
        const auto extraLength = static_cast<std::size_t>(little16(bytes, position));
        position += 2U;
        if (position > bytes.size() || bytes.size() - position < extraLength) {
            return std::unexpected("truncated gzip extra field");
        }
        position += extraLength;
    }
    if ((flags & 0x08U) != 0U) {
        auto status = advanceZeroTerminated(
            bytes, position, memberStart, limits, "original filename");
        if (!status) return std::unexpected(status.error());
    }
    if ((flags & 0x10U) != 0U) {
        auto status = advanceZeroTerminated(bytes, position, memberStart, limits, "comment");
        if (!status) return std::unexpected(status.error());
    }
    if (position - memberStart > limits.maxHeaderBytes) {
        return std::unexpected("gzip header exceeds the configured header-byte limit");
    }
    if ((flags & 0x02U) != 0U) {
        if (position > bytes.size() || bytes.size() - position < 2U) {
            return std::unexpected("truncated gzip header CRC");
        }
        const auto expected = little16(bytes, position);
        const auto computed = static_cast<std::uint16_t>(
            (crc32Range(bytes.substr(memberStart, position - memberStart)) ^ 0xffffffffU) & 0xffffU);
        if (computed != expected) return std::unexpected("gzip header CRC mismatch");
        position += 2U;
    }
    if (position - memberStart > limits.maxHeaderBytes) {
        return std::unexpected("gzip header exceeds the configured header-byte limit");
    }
    return position;
}

} // namespace

auto hasGzipMagic(const std::string_view bytes) noexcept -> bool {
    return bytes.size() >= 2U && byteAt(bytes, 0U) == gzipId1 && byteAt(bytes, 1U) == gzipId2;
}

auto crc32(const std::string_view bytes) noexcept -> std::uint32_t {
    return crc32Range(bytes) ^ 0xffffffffU;
}

auto gzipDecompress(const std::string_view bytes, const GzipLimits& limits)
    -> std::expected<std::string, std::string> {
    if (bytes.empty()) return std::unexpected("gzip input is empty");
    if (bytes.size() > limits.maxCompressedBytes) {
        return std::unexpected("gzip input exceeds the configured compressed-byte limit");
    }
    if (limits.maxMembers == 0U) {
        return std::unexpected("gzip member limit is zero");
    }

    std::string output;
    output.reserve(std::min(bytes.size() * (bytes.size() <= limits.maxDecompressedBytes / 2U ? 2U : 1U),
                            limits.maxDecompressedBytes));
    std::size_t position{};
    std::size_t memberCount{};
    while (position < bytes.size()) {
        if (++memberCount > limits.maxMembers) {
            return std::unexpected("gzip input exceeds the configured member limit");
        }
        const auto memberStart = position;
        auto deflateStart = parseHeader(bytes, memberStart, limits);
        if (!deflateStart) {
            return std::unexpected("gzip member " + std::to_string(memberCount) + ": " +
                                   deflateStart.error());
        }
        const auto outputStart = output.size();
        auto trailerStart = decodeDeflate(bytes, *deflateStart, output, outputStart, limits);
        if (!trailerStart) {
            return std::unexpected("gzip member " + std::to_string(memberCount) + ": " +
                                   trailerStart.error());
        }
        if (*trailerStart > bytes.size() || bytes.size() - *trailerStart < 8U) {
            return std::unexpected("gzip member " + std::to_string(memberCount) +
                                   ": truncated trailer");
        }
        const auto expectedCrc = little32(bytes, *trailerStart);
        const auto expectedSize = little32(bytes, *trailerStart + 4U);
        const std::string_view memberOutput(output.data() + outputStart, output.size() - outputStart);
        if (crc32(memberOutput) != expectedCrc) {
            return std::unexpected("gzip member " + std::to_string(memberCount) +
                                   ": payload CRC32 mismatch");
        }
        if (static_cast<std::uint32_t>(memberOutput.size() & 0xffffffffU) != expectedSize) {
            return std::unexpected("gzip member " + std::to_string(memberCount) +
                                   ": uncompressed size mismatch");
        }
        position = *trailerStart + 8U;
        if (position < bytes.size() && !hasGzipMagic(bytes.substr(position))) {
            return std::unexpected("trailing bytes after gzip member are not another member");
        }
    }
    return output;
}

auto gzipCompress(const std::string_view bytes) -> std::expected<std::string, std::string> {
    constexpr std::size_t blockMaximum{65535U};
    const auto blocks = std::max<std::size_t>(
        1U, bytes.size() / blockMaximum + (bytes.size() % blockMaximum != 0U ? 1U : 0U));
    constexpr std::size_t wrapperBytes{18U};
    if (bytes.size() > std::numeric_limits<std::size_t>::max() - wrapperBytes) {
        return std::unexpected("gzip output size is not representable");
    }
    if (blocks > (std::numeric_limits<std::size_t>::max() - wrapperBytes - bytes.size()) / 5U) {
        return std::unexpected("gzip output size is not representable");
    }
    std::string output;
    output.reserve(wrapperBytes + bytes.size() + blocks * 5U);
    output.push_back(static_cast<char>(gzipId1));
    output.push_back(static_cast<char>(gzipId2));
    output.push_back(static_cast<char>(deflateMethod));
    output.push_back('\0');
    append32(output, 0U);
    output.push_back('\0');
    output.push_back(static_cast<char>(0xffU));

    std::size_t position{};
    do {
        const auto length = std::min(blockMaximum, bytes.size() - position);
        const auto final = position + length == bytes.size();
        output.push_back(final ? '\x01' : '\0');
        append16(output, static_cast<std::uint16_t>(length));
        append16(output, static_cast<std::uint16_t>(static_cast<std::uint16_t>(length) ^ 0xffffU));
        output.append(bytes.substr(position, length));
        position += length;
    } while (position < bytes.size());

    append32(output, crc32(bytes));
    append32(output, static_cast<std::uint32_t>(bytes.size() & 0xffffffffU));
    return output;
}

} // namespace intarnanew
