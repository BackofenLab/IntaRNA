#pragma once

#include "intarnanew/types.hpp"

#include <cstdint>
#include <expected>
#include <istream>
#include <string>
#include <string_view>
#include <vector>

namespace intarnanew {

class Sequence {
public:
    Sequence(std::string identifier, std::string nucleotides, long long firstPosition = 1);

    [[nodiscard]] auto id() const noexcept -> const std::string& { return identifier_; }
    [[nodiscard]] auto str() const noexcept -> const std::string& { return nucleotides_; }
    [[nodiscard]] auto size() const noexcept -> Index { return nucleotides_.size(); }
    [[nodiscard]] auto empty() const noexcept -> bool { return nucleotides_.empty(); }
    [[nodiscard]] auto operator[](Index index) const -> char { return nucleotides_.at(index); }
    [[nodiscard]] auto externalIndex(Index index) const -> long long;
    [[nodiscard]] auto internalIndex(long long external) const -> std::expected<Index, std::string>;
    [[nodiscard]] auto firstPosition() const noexcept -> long long { return firstPosition_; }

private:
    std::string identifier_;
    std::string nucleotides_;
    long long firstPosition_;
};

class SequenceReader {
public:
    // File and standard-input payloads are decoded as gzip when the RFC 1952
    // signature is present; a .gz file suffix requires that signature.
    [[nodiscard]] static auto read(
        std::string_view specification,
        std::string_view fallbackId,
        long long firstPosition,
        std::istream& standardInput) -> std::expected<std::vector<Sequence>, std::string>;

    [[nodiscard]] static auto parseFasta(
        std::istream& input,
        std::string_view fallbackId,
        long long firstPosition) -> std::expected<std::vector<Sequence>, std::string>;

    // Selects one-based FASTA record indices while preserving input order.
    // The grammar is a comma-separated list of indices or closed ranges (e.g. 1,3-5).
    [[nodiscard]] static auto select(
        std::vector<Sequence> sequences,
        std::string_view specification) -> std::expected<std::vector<Sequence>, std::string>;
};

[[nodiscard]] auto nucleotideMask(char symbol) noexcept -> std::uint8_t;
[[nodiscard]] auto canPair(char target, char query, bool allowGu = true) noexcept -> bool;
[[nodiscard]] auto isGuPair(char target, char query) noexcept -> bool;
[[nodiscard]] auto reverseComplement(std::string_view sequence) -> std::string;

} // namespace intarnanew
