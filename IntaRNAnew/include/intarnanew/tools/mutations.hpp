#pragma once

#include "intarnanew/sequence.hpp"

#include <expected>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace intarnanew::tools {

enum class MutationGenerator { flip, any };
enum class CandidateFilter { removeGu, removeAu, removeCg };

struct MutationCandidate {
    Index queryIndex{};
    Index targetIndex{};
    char wildQuery{};
    char wildTarget{};
    char mutatedQuery{};
    char mutatedTarget{};

    [[nodiscard]] auto encoding(const Sequence& query, const Sequence& target) const -> std::string;
    friend auto operator==(const MutationCandidate&, const MutationCandidate&) -> bool = default;
};

// Enumerates compensatory point mutations for supplied intermolecular pairs.
// flip swaps the wild-type bases (GC -> CG); any enumerates every canonical or
// GU pair for which both bases differ from wild type, in lexical A,C,G,U order.
// Filters remove candidate wild-type pair classes, matching the public CopomuS
// terminology.
[[nodiscard]] auto enumerateMutations(
    const Sequence& query,
    const Sequence& target,
    std::span<const BasePair> interactionPairs,
    MutationGenerator generator,
    std::span<const CandidateFilter> filters = {})
    -> std::expected<std::vector<MutationCandidate>, std::string>;

// Parses Q-base/index/mutant '&' target-base/index/mutant, e.g. G1C&U7G.
// The encoded wild-type bases must match the sequences at their external
// coordinates, and the two mutant bases must form a canonical or GU pair.
[[nodiscard]] auto parseMutationEncoding(
    std::string_view encoding,
    const Sequence& query,
    const Sequence& target) -> std::expected<MutationCandidate, std::string>;

struct MutationSequences {
    std::string wildQuery;
    std::string wildTarget;
    std::string mutatedQuery;
    std::string mutatedTarget;
};

[[nodiscard]] auto applyMutation(
    const Sequence& query,
    const Sequence& target,
    const MutationCandidate& mutation) -> std::expected<MutationSequences, std::string>;

// Generates a mono-nucleotide preserving shuffle using a stable, specified
// SplitMix64/Fisher-Yates implementation; results are cross-platform for a seed.
[[nodiscard]] auto shuffleMononucleotides(std::string_view sequence, std::uint64_t seed)
    -> std::expected<std::string, std::string>;

// Generates a dinucleotide-preserving random Euler trail over the sequence's
// directed adjacency multigraph. It preserves the first and last base plus
// exact directed dinucleotide counts. Outgoing edge choices are deterministically
// shuffled before a complete Euler trail is constructed.
[[nodiscard]] auto shuffleDinucleotides(std::string_view sequence, std::uint64_t seed)
    -> std::expected<std::string, std::string>;

} // namespace intarnanew::tools
