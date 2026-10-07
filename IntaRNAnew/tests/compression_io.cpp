#include "intarnanew/compression.hpp"
#include "intarnanew/output.hpp"
#include "intarnanew/sequence.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

int failureCount{};

void check(const bool condition, const std::string_view description) {
    if (!condition) {
        std::cerr << "FAILED: " << description << '\n';
        ++failureCount;
    }
}

[[nodiscard]] auto independentCrc32(const std::string_view bytes) noexcept -> std::uint32_t {
    std::uint32_t result{0xffffffffU};
    for (const unsigned char byte : bytes) {
        result ^= byte;
        for (unsigned bit = 0U; bit < 8U; ++bit) {
            result = (result >> 1U) ^ ((result & 1U) != 0U ? 0xedb88320U : 0U);
        }
    }
    return result ^ 0xffffffffU;
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

class IndependentBitWriter {
public:
    void bit(const unsigned value) {
        if (bitPosition_ == 0U) bytes_.push_back('\0');
        bytes_.back() = static_cast<char>(
            static_cast<unsigned char>(bytes_.back()) |
            static_cast<unsigned char>((value & 1U) << bitPosition_));
        bitPosition_ = (bitPosition_ + 1U) % 8U;
    }

    void leastSignificantFirst(const unsigned value, const unsigned count) {
        for (unsigned index = 0U; index < count; ++index) bit(value >> index);
    }

    void huffman(const unsigned symbol, const std::vector<std::uint8_t>& lengths) {
        constexpr unsigned maximumLength{15U};
        std::array<unsigned, maximumLength + 1U> counts{};
        for (const auto length : lengths) {
            if (length > maximumLength) throw std::logic_error("test code length is too large");
            if (length != 0U) ++counts[length];
        }
        std::array<unsigned, maximumLength + 1U> next{};
        unsigned code{};
        for (unsigned length = 1U; length <= maximumLength; ++length) {
            code = (code + counts[length - 1U]) << 1U;
            next[length] = code;
        }
        for (unsigned candidate = 0U; candidate < symbol; ++candidate) {
            if (lengths[candidate] != 0U) ++next[lengths[candidate]];
        }
        const auto length = static_cast<unsigned>(lengths.at(symbol));
        if (length == 0U) throw std::logic_error("test symbol has no Huffman code");
        const auto assigned = next[length];
        for (unsigned remaining = length; remaining != 0U; --remaining) {
            bit(assigned >> (remaining - 1U));
        }
    }

    [[nodiscard]] auto take() && -> std::string { return std::move(bytes_); }

private:
    std::string bytes_;
    unsigned bitPosition_{};
};

[[nodiscard]] auto fixedLiteralLengths() -> std::vector<std::uint8_t> {
    std::vector<std::uint8_t> lengths(288U, 0U);
    std::fill(lengths.begin(), lengths.begin() + 144, 8U);
    std::fill(lengths.begin() + 144, lengths.begin() + 256, 9U);
    std::fill(lengths.begin() + 256, lengths.begin() + 280, 7U);
    std::fill(lengths.begin() + 280, lengths.end(), 8U);
    return lengths;
}

[[nodiscard]] auto fixedDeflate() -> std::string {
    IndependentBitWriter output;
    output.leastSignificantFirst(1U, 1U); // final block
    output.leastSignificantFirst(1U, 2U); // fixed Huffman block
    const auto literals = fixedLiteralLengths();
    const std::vector<std::uint8_t> distances(32U, 5U);
    output.huffman(static_cast<unsigned>('A'), literals);
    output.huffman(static_cast<unsigned>('B'), literals);
    output.huffman(static_cast<unsigned>('C'), literals);
    output.huffman(260U, literals); // six bytes
    output.huffman(2U, distances);  // distance three
    output.huffman(256U, literals);
    return std::move(output).take();
}

[[nodiscard]] auto dynamicDeflate() -> std::string {
    IndependentBitWriter output;
    output.leastSignificantFirst(1U, 1U);  // final block
    output.leastSignificantFirst(2U, 2U);  // dynamic Huffman block
    output.leastSignificantFirst(3U, 5U);  // 260 literal/length symbols
    output.leastSignificantFirst(0U, 5U);  // one distance symbol
    output.leastSignificantFirst(14U, 4U); // 18 code-length symbols

    constexpr std::array<unsigned, 19U> order{
        16U, 17U, 18U, 0U, 8U, 7U, 9U, 6U, 10U, 5U,
        11U, 4U, 12U, 3U, 13U, 2U, 14U, 1U, 15U,
    };
    std::vector<std::uint8_t> codeLengths(19U, 0U);
    codeLengths[0U] = 2U;
    codeLengths[1U] = 2U;
    codeLengths[2U] = 2U;
    codeLengths[18U] = 2U;
    for (std::size_t index = 0U; index < 18U; ++index) {
        output.leastSignificantFirst(codeLengths[order[index]], 3U);
    }

    output.huffman(18U, codeLengths); // 65 zero lengths
    output.leastSignificantFirst(54U, 7U);
    output.huffman(1U, codeLengths);  // literal A: length one
    output.huffman(18U, codeLengths); // 138 zero lengths
    output.leastSignificantFirst(127U, 7U);
    output.huffman(18U, codeLengths); // 52 zero lengths
    output.leastSignificantFirst(41U, 7U);
    output.huffman(2U, codeLengths);  // end-of-block: length two
    output.huffman(0U, codeLengths);  // symbols 257 and 258 absent
    output.huffman(0U, codeLengths);
    output.huffman(2U, codeLengths);  // length-five symbol: length two
    output.huffman(1U, codeLengths);  // distance-one symbol: length one

    std::vector<std::uint8_t> literals(260U, 0U);
    literals[65U] = 1U;
    literals[256U] = 2U;
    literals[259U] = 2U;
    const std::vector<std::uint8_t> distances{1U};
    output.huffman(65U, literals);
    output.huffman(259U, literals); // five-byte match
    output.huffman(0U, distances);  // distance one
    output.huffman(256U, literals);
    return std::move(output).take();
}

struct IndependentMember {
    std::string bytes;
    std::size_t headerCrcOffset{std::string::npos};
};

[[nodiscard]] auto member(
    const std::string_view deflate,
    const std::string_view uncompressed,
    const bool optionalHeader = false) -> IndependentMember {
    IndependentMember result;
    result.bytes.push_back(static_cast<char>(0x1fU));
    result.bytes.push_back(static_cast<char>(0x8bU));
    result.bytes.push_back(static_cast<char>(8U));
    result.bytes.push_back(static_cast<char>(optionalHeader ? 0x1eU : 0U));
    append32(result.bytes, 0x12345678U);
    result.bytes.push_back('\0');
    result.bytes.push_back(static_cast<char>(3U));
    if (optionalHeader) {
        append16(result.bytes, 3U);
        result.bytes.append("xyz", 3U);
        result.bytes.append("fixture", 7U);
        result.bytes.push_back('\0');
        result.bytes.append("generated", 9U);
        result.bytes.push_back('\0');
        result.headerCrcOffset = result.bytes.size();
        append16(result.bytes, static_cast<std::uint16_t>(independentCrc32(result.bytes)));
    }
    result.bytes.append(deflate);
    append32(result.bytes, independentCrc32(uncompressed));
    append32(result.bytes, static_cast<std::uint32_t>(uncompressed.size() & 0xffffffffU));
    return result;
}

[[nodiscard]] auto readAll(const std::filesystem::path& path) -> std::string {
    std::ifstream input(path, std::ios::binary);
    return std::string(
        std::istreambuf_iterator<char>{input}, std::istreambuf_iterator<char>{});
}

void writeAll(const std::filesystem::path& path, const std::string_view bytes) {
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    output.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
}

[[nodiscard]] auto sameSequences(
    const std::vector<intarnanew::Sequence>& left,
    const std::vector<intarnanew::Sequence>& right) -> bool {
    if (left.size() != right.size()) return false;
    for (std::size_t index = 0U; index < left.size(); ++index) {
        if (left[index].id() != right[index].id() ||
            left[index].str() != right[index].str() ||
            left[index].firstPosition() != right[index].firstPosition()) {
            return false;
        }
    }
    return true;
}

} // namespace

auto main() -> int {
    using namespace intarnanew;

    check(crc32("123456789") == 0xcbf43926U,
          "CRC32 matches the standard check vector");

    const std::string knownHello{
        "\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\x03"
        "\xcb\x48\xcd\xc9\xc9\xe7\x02\x00"
        "\x20\x30\x3a\x36\x06\x00\x00\x00",
        26U};
    auto knownDecoded = gzipDecompress(knownHello);
    check(knownDecoded && *knownDecoded == "hello\n",
          "decoder accepts an externally specified canonical gzip check member");

    const auto fixed = member(fixedDeflate(), "ABCABCABC").bytes;
    auto fixedDecoded = gzipDecompress(fixed);
    check(fixedDecoded && *fixedDecoded == "ABCABCABC",
          "fixed-Huffman DEFLATE literals and back-reference decode");

    const auto dynamicMember = member(dynamicDeflate(), "AAAAAA", true);
    auto dynamicDecoded = gzipDecompress(dynamicMember.bytes);
    check(dynamicDecoded && *dynamicDecoded == "AAAAAA",
          "dynamic-Huffman DEFLATE and all optional gzip header fields decode");
    bool everyOptionalHeaderTruncationRejected = true;
    for (std::size_t size = 0U; size < dynamicMember.bytes.size(); ++size) {
        everyOptionalHeaderTruncationRejected = everyOptionalHeaderTruncationRejected &&
            !gzipDecompress(std::string_view(dynamicMember.bytes).substr(0U, size));
    }
    check(everyOptionalHeaderTruncationRejected,
          "every truncation of a member with optional header fields is rejected");

    auto concatenated = gzipDecompress(fixed + dynamicMember.bytes);
    check(concatenated && *concatenated == "ABCABCABCAAAAAA",
          "concatenated gzip members decode in order");

    std::string binary;
    binary.reserve(150000U);
    for (std::size_t index = 0U; index < 150000U; ++index) {
        binary.push_back(static_cast<char>((index * 73U + 19U) & 0xffU));
    }
    auto encoded = gzipCompress(binary);
    check(encoded.has_value(), "stored-block gzip encoder accepts a multi-block binary payload");
    auto decoded = encoded ? gzipDecompress(*encoded)
                           : std::expected<std::string, std::string>{std::unexpected("encode failed")};
    check(decoded && *decoded == binary, "stored-block gzip round-trip preserves every byte");
    auto encodedAgain = gzipCompress(binary);
    check(encoded && encodedAgain && *encoded == *encodedAgain,
          "gzip encoder is deterministic");

    auto emptyEncoded = gzipCompress("");
    auto emptyDecoded = emptyEncoded ? gzipDecompress(*emptyEncoded)
                                     : std::expected<std::string, std::string>{
                                           std::unexpected("encode failed")};
    check(emptyDecoded && emptyDecoded->empty(), "empty gzip member round-trips");

    auto truncationMember = gzipCompress("truncation fixture");
    bool everyTruncationRejected = truncationMember.has_value();
    if (truncationMember) {
        for (std::size_t size = 0U; size < truncationMember->size(); ++size) {
            everyTruncationRejected = everyTruncationRejected &&
                !gzipDecompress(std::string_view(*truncationMember).substr(0U, size));
        }
    }
    check(everyTruncationRejected, "every truncation of a valid gzip member is rejected");

    if (encoded) {
        auto badCrc = *encoded;
        badCrc[badCrc.size() - 8U] = static_cast<char>(
            static_cast<unsigned char>(badCrc[badCrc.size() - 8U]) ^ 0x80U);
        auto badCrcResult = gzipDecompress(badCrc);
        check(!badCrcResult && badCrcResult.error().find("CRC32") != std::string::npos,
              "corrupt payload CRC32 is diagnosed");

        auto badSize = *encoded;
        badSize.back() = static_cast<char>(static_cast<unsigned char>(badSize.back()) ^ 1U);
        auto badSizeResult = gzipDecompress(badSize);
        check(!badSizeResult && badSizeResult.error().find("size") != std::string::npos,
              "corrupt gzip ISIZE is diagnosed");

        auto badStoredLength = *encoded;
        badStoredLength[13U] = static_cast<char>(
            static_cast<unsigned char>(badStoredLength[13U]) ^ 1U);
        auto badLengthResult = gzipDecompress(badStoredLength);
        check(!badLengthResult && badLengthResult.error().find("LEN/NLEN") != std::string::npos,
              "corrupt stored-block length complement is diagnosed");

        GzipLimits compressedLimit;
        compressedLimit.maxCompressedBytes = encoded->size() - 1U;
        check(!gzipDecompress(*encoded, compressedLimit),
              "compressed-byte limit is enforced");

        GzipLimits decompressedLimit;
        decompressedLimit.maxDecompressedBytes = binary.size() - 1U;
        check(!gzipDecompress(*encoded, decompressedLimit),
              "decompressed-byte limit is enforced");

        GzipLimits blockLimit;
        blockLimit.maxDeflateBlocks = 1U;
        check(!gzipDecompress(*encoded, blockLimit), "DEFLATE block limit is enforced");
    }

    GzipLimits memberLimit;
    memberLimit.maxMembers = 1U;
    check(!gzipDecompress(fixed + dynamicMember.bytes, memberLimit),
          "concatenated-member limit is enforced");

    GzipLimits headerLimit;
    headerLimit.maxHeaderBytes = 16U;
    check(!gzipDecompress(dynamicMember.bytes, headerLimit),
          "optional gzip header-byte limit is enforced");
    headerLimit.maxHeaderBytes = dynamicMember.headerCrcOffset;
    check(!gzipDecompress(dynamicMember.bytes, headerLimit),
          "gzip header-byte limit includes the optional header CRC itself");

    auto badHeaderCrc = dynamicMember.bytes;
    badHeaderCrc[dynamicMember.headerCrcOffset] = static_cast<char>(
        static_cast<unsigned char>(badHeaderCrc[dynamicMember.headerCrcOffset]) ^ 1U);
    auto headerCrcResult = gzipDecompress(badHeaderCrc);
    check(!headerCrcResult && headerCrcResult.error().find("header CRC") != std::string::npos,
          "corrupt optional gzip header CRC is diagnosed");

    auto reservedFlags = fixed;
    reservedFlags[3U] = static_cast<char>(0x20U);
    check(!gzipDecompress(reservedFlags), "reserved gzip header flags are rejected");

    auto trailing = fixed + "junk";
    auto trailingResult = gzipDecompress(trailing);
    check(!trailingResult && trailingResult.error().find("trailing bytes") != std::string::npos,
          "non-member trailing bytes are rejected");

    auto reservedBlock = member(std::string(1U, '\x07'), "").bytes;
    auto reservedBlockResult = gzipDecompress(reservedBlock);
    check(!reservedBlockResult && reservedBlockResult.error().find("reserved block") != std::string::npos,
          "reserved DEFLATE block type is rejected");

    {
        std::istringstream consecutive(">empty\n>valid\nACGU\n");
        auto parsed = SequenceReader::parseFasta(consecutive, "fallback", 1);
        check(!parsed && parsed.error().find("empty") != std::string::npos,
              "empty FASTA record before another header is rejected");

        std::istringstream trailingHeader(">valid\nACGU\n>empty\n");
        check(!SequenceReader::parseFasta(trailingHeader, "fallback", 1),
              "empty trailing FASTA record is rejected");
    }

    {
        bool spanRejected{};
        try {
            const Sequence invalid(
                "overflow", "AA", std::numeric_limits<long long>::max());
            (void)invalid;
        } catch (const std::out_of_range&) {
            spanRejected = true;
        }
        check(spanRejected, "coordinate span beyond LLONG_MAX is rejected at construction");

        const Sequence maximum("maximum", "A", std::numeric_limits<long long>::max());
        check(maximum.externalIndex(0U) == std::numeric_limits<long long>::max(),
              "maximum signed coordinate remains representable");
        bool indexRejected{};
        try {
            (void)maximum.externalIndex(1U);
        } catch (const std::out_of_range&) {
            indexRejected = true;
        }
        check(indexRejected, "external coordinate conversion rejects an out-of-range index");

        const Sequence minimum("minimum", "AA", std::numeric_limits<long long>::min());
        check(minimum.externalIndex(0U) == std::numeric_limits<long long>::min() &&
              minimum.externalIndex(1U) == std::numeric_limits<long long>::min() + 1LL,
              "minimum signed origin advances without signed overflow");
    }

    {
        const auto unique = std::to_string(
            std::chrono::steady_clock::now().time_since_epoch().count());
        const auto directory = std::filesystem::temp_directory_path() /
                               ("intarnanew-gzip-test-" + unique);
        std::error_code error;
        std::filesystem::create_directory(directory, error);
        check(!error, "gzip I/O test directory was created");

        const std::string fasta{">one description\nACGT\n>two\nNNNN\n"};
        const auto plainPath = directory / "input.fa";
        const auto gzipPath = directory / "input.fa.gz";
        const auto magicPath = directory / "input.bin";
        writeAll(plainPath, fasta);
        writeAll(gzipPath, "old destination contents");
        auto writeStatus = writeOutput(gzipPath.string(), fasta);
        check(writeStatus.has_value(), "transactional .gz output replaces a regular file");
        const auto compressedBytes = readAll(gzipPath);
        check(hasGzipMagic(compressedBytes), ".gz output has the gzip signature");
        auto writtenPayload = gzipDecompress(compressedBytes);
        check(writtenPayload && *writtenPayload == fasta,
              "transactional .gz output decompresses to the requested content");

        std::istringstream unused;
        auto plainSequences = SequenceReader::read(plainPath.string(), "fallback", -2, unused);
        auto gzipSequences = SequenceReader::read(gzipPath.string(), "fallback", -2, unused);
        check(plainSequences && gzipSequences && sameSequences(*plainSequences, *gzipSequences),
              "plain and gzip FASTA files parse to byte-equivalent sequences");
        if (plainSequences && plainSequences->size() == 2U) {
            check((*plainSequences)[0].str() == "ACGU" && (*plainSequences)[1].str() == "NNNN",
                  "gzip integration preserves normal FASTA normalization");
        }

        writeAll(magicPath, compressedBytes);
        auto magicSequences = SequenceReader::read(magicPath.string(), "fallback", -2, unused);
        check(gzipSequences && magicSequences && sameSequences(*gzipSequences, *magicSequences),
              "gzip input is detected by signature without relying on its suffix");

        std::istringstream gzipInput(compressedBytes);
        auto stdinSequences = SequenceReader::read("STDIN", "fallback", -2, gzipInput);
        check(gzipSequences && stdinSequences && sameSequences(*gzipSequences, *stdinSequences),
              "gzip-compressed standard input parses identically to a gzip file");

        const auto falseGzipPath = directory / "plain.fa.gz";
        writeAll(falseGzipPath, fasta);
        auto falseGzip = SequenceReader::read(falseGzipPath.string(), "fallback", 1, unused);
        check(!falseGzip && falseGzip.error().find("no gzip signature") != std::string::npos,
              ".gz input without a gzip signature has a clear diagnostic");

        const auto corruptPath = directory / "corrupt.fa.gz";
        auto corruptBytes = compressedBytes;
        corruptBytes[corruptBytes.size() - 8U] = static_cast<char>(
            static_cast<unsigned char>(corruptBytes[corruptBytes.size() - 8U]) ^ 1U);
        writeAll(corruptPath, corruptBytes);
        auto corruptInput = SequenceReader::read(corruptPath.string(), "fallback", 1, unused);
        check(!corruptInput && corruptInput.error().find("CRC32") != std::string::npos,
              "corrupt gzip FASTA reports the member validation failure");

        const auto secondPath = directory / "second.fa.gz";
        auto secondWrite = writeOutput(secondPath.string(), fasta);
        check(secondWrite && readAll(secondPath) == compressedBytes,
              "file-level gzip output is deterministic across destinations");

        const auto directoryDestination = directory / "not-a-file.gz";
        std::filesystem::create_directory(directoryDestination, error);
        auto rejectedOutput = writeOutput(directoryDestination.string(), fasta);
        check(!rejectedOutput && std::filesystem::is_directory(directoryDestination),
              "failed transactional .gz commit leaves a directory destination intact");

        std::filesystem::remove_all(directory, error);
    }

    if (failureCount == 0) {
        std::cout << "All IntaRNAnew compression and gzip I/O tests passed.\n";
    }
    return failureCount == 0 ? 0 : 1;
}
