#pragma once

#include <cstddef>
#include <cstdint>
#include <expected>
#include <string>
#include <string_view>

namespace intarnanew {

struct GzipLimits {
    // Byte limits apply to the complete input and concatenated output. Header
    // and block limits apply independently to each member.
    std::size_t maxCompressedBytes{1U << 30U};
    std::size_t maxDecompressedBytes{1U << 30U};
    std::size_t maxHeaderBytes{1U << 20U};
    std::size_t maxMembers{1024U};
    std::size_t maxDeflateBlocks{1U << 24U};
};

[[nodiscard]] auto hasGzipMagic(std::string_view bytes) noexcept -> bool;

// Decodes all concatenated RFC 1952 members. The returned string contains the
// members' uncompressed payloads in input order.
[[nodiscard]] auto gzipDecompress(
    std::string_view bytes,
    const GzipLimits& limits = {}) -> std::expected<std::string, std::string>;

// Emits one deterministic RFC 1952 member using stored RFC 1951 blocks.
[[nodiscard]] auto gzipCompress(std::string_view bytes)
    -> std::expected<std::string, std::string>;

[[nodiscard]] auto crc32(std::string_view bytes) noexcept -> std::uint32_t;

} // namespace intarnanew
