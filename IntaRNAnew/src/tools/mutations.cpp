#include "intarnanew/tools/mutations.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cctype>
#include <cstdint>
#include <limits>
#include <set>

namespace intarnanew::tools {
namespace {

[[nodiscard]] auto normalizedBase(const char symbol) noexcept -> char {
    const auto upper = static_cast<char>(std::toupper(static_cast<unsigned char>(symbol)));
    return upper == 'T' ? 'U' : upper;
}

[[nodiscard]] auto pairClassFiltered(
    const char query,
    const char target,
    const std::span<const CandidateFilter> filters) noexcept -> bool {
    const auto q = normalizedBase(query);
    const auto t = normalizedBase(target);
    for (const auto filter : filters) {
        if (filter == CandidateFilter::removeGu && isGuPair(t, q)) {
            return true;
        }
        if (filter == CandidateFilter::removeAu &&
            ((q == 'A' && t == 'U') || (q == 'U' && t == 'A'))) {
            return true;
        }
        if (filter == CandidateFilter::removeCg &&
            ((q == 'C' && t == 'G') || (q == 'G' && t == 'C'))) {
            return true;
        }
    }
    return false;
}

[[nodiscard]] auto baseName(const char value, const long long coordinate, const char mutation)
    -> std::string {
    return std::string(1U, value) + std::to_string(coordinate) + std::string(1U, mutation);
}

struct ParsedHalf {
    char wild{};
    long long index{};
    char mutated{};
};

[[nodiscard]] auto parseHalf(const std::string_view text) -> std::expected<ParsedHalf, std::string> {
    if (text.size() < 3U) {
        return std::unexpected("mutation half '" + std::string{text} + "' is too short");
    }
    ParsedHalf result;
    result.wild = normalizedBase(text.front());
    result.mutated = normalizedBase(text.back());
    constexpr std::string_view bases{"ACGU"};
    if (bases.find(result.wild) == std::string_view::npos ||
        bases.find(result.mutated) == std::string_view::npos) {
        return std::unexpected("mutation bases must be A, C, G, or U");
    }
    const auto coordinate = text.substr(1U, text.size() - 2U);
    const auto [end, error] = std::from_chars(
        coordinate.data(), coordinate.data() + coordinate.size(), result.index);
    if (error != std::errc{} || end != coordinate.data() + coordinate.size()) {
        return std::unexpected("invalid mutation coordinate '" + std::string{coordinate} + "'");
    }
    return result;
}

class SplitMix64 {
public:
    explicit SplitMix64(const std::uint64_t seed) : state_(seed) {}

    [[nodiscard]] auto next() noexcept -> std::uint64_t {
        std::uint64_t result = (state_ += 0x9E3779B97F4A7C15ULL);
        result = (result ^ (result >> 30U)) * 0xBF58476D1CE4E5B9ULL;
        result = (result ^ (result >> 27U)) * 0x94D049BB133111EBULL;
        return result ^ (result >> 31U);
    }

    [[nodiscard]] auto bounded(const std::uint64_t bound) noexcept -> std::uint64_t {
        if (bound <= 1U) {
            return 0U;
        }
        const auto threshold = static_cast<std::uint64_t>(-bound) % bound;
        for (;;) {
            const auto value = next();
            if (value >= threshold) {
                return value % bound;
            }
        }
    }

private:
    std::uint64_t state_{};
};

template <typename Value>
void shuffle(std::vector<Value>& values, SplitMix64& random) {
    for (std::size_t remaining = values.size(); remaining > 1U; --remaining) {
        const auto chosen = static_cast<std::size_t>(random.bounded(remaining));
        std::swap(values[remaining - 1U], values[chosen]);
    }
}

[[nodiscard]] auto canonicalSequence(const std::string_view sequence)
    -> std::expected<std::string, std::string> {
    if (sequence.empty()) {
        return std::unexpected("sequence must not be empty");
    }
    std::string result;
    result.reserve(sequence.size());
    for (const char symbol : sequence) {
        const char base = normalizedBase(symbol);
        if (std::string_view{"ACGU"}.find(base) == std::string_view::npos) {
            return std::unexpected("shuffle input must contain only A, C, G, U, or T");
        }
        result.push_back(base);
    }
    return result;
}

[[nodiscard]] auto baseIndex(const char base) noexcept -> std::size_t {
    switch (base) {
        case 'A': return 0U;
        case 'C': return 1U;
        case 'G': return 2U;
        default: return 3U;
    }
}

} // namespace

auto MutationCandidate::encoding(const Sequence& query, const Sequence& target) const -> std::string {
    if (queryIndex >= query.size() || targetIndex >= target.size()) {
        return {};
    }
    return baseName(wildQuery, query.externalIndex(queryIndex), mutatedQuery) + '&' +
           baseName(wildTarget, target.externalIndex(targetIndex), mutatedTarget);
}

auto enumerateMutations(
    const Sequence& query,
    const Sequence& target,
    const std::span<const BasePair> interactionPairs,
    const MutationGenerator generator,
    const std::span<const CandidateFilter> filters)
    -> std::expected<std::vector<MutationCandidate>, std::string> {
    constexpr std::array<char, 4U> bases{'A', 'C', 'G', 'U'};
    std::vector<MutationCandidate> result;
    std::set<std::array<std::size_t, 4U>> seen;
    for (const auto pair : interactionPairs) {
        if (pair.query >= query.size() || pair.target >= target.size()) {
            return std::unexpected("interaction pair is outside its sequence");
        }
        const char wildQuery = query[pair.query];
        const char wildTarget = target[pair.target];
        if (!canPair(wildTarget, wildQuery)) {
            return std::unexpected("interaction contains a non-pairing base combination");
        }
        if (pairClassFiltered(wildQuery, wildTarget, filters)) {
            continue;
        }
        if (generator == MutationGenerator::flip) {
            const char mutatedQuery = wildTarget;
            const char mutatedTarget = wildQuery;
            if (mutatedQuery != wildQuery && mutatedTarget != wildTarget &&
                canPair(mutatedTarget, mutatedQuery)) {
                const auto key = std::array{
                    pair.query, pair.target,
                    static_cast<std::size_t>(mutatedQuery),
                    static_cast<std::size_t>(mutatedTarget)};
                if (seen.insert(key).second) {
                    result.push_back({
                        pair.query, pair.target, wildQuery, wildTarget,
                        mutatedQuery, mutatedTarget});
                }
            }
            continue;
        }
        for (const char mutatedQuery : bases) {
            for (const char mutatedTarget : bases) {
                if (mutatedQuery == wildQuery || mutatedTarget == wildTarget ||
                    !canPair(mutatedTarget, mutatedQuery)) {
                    continue;
                }
                const auto key = std::array{
                    pair.query, pair.target,
                    static_cast<std::size_t>(mutatedQuery),
                    static_cast<std::size_t>(mutatedTarget)};
                if (seen.insert(key).second) {
                    result.push_back({
                        pair.query, pair.target, wildQuery, wildTarget,
                        mutatedQuery, mutatedTarget});
                }
            }
        }
    }
    return result;
}

auto parseMutationEncoding(
    const std::string_view encoding,
    const Sequence& query,
    const Sequence& target) -> std::expected<MutationCandidate, std::string> {
    const auto separator = encoding.find('&');
    if (separator == std::string_view::npos ||
        encoding.find('&', separator + 1U) != std::string_view::npos) {
        return std::unexpected("mutation encoding must contain exactly one '&'");
    }
    auto queryHalf = parseHalf(encoding.substr(0U, separator));
    auto targetHalf = parseHalf(encoding.substr(separator + 1U));
    if (!queryHalf) return std::unexpected(queryHalf.error());
    if (!targetHalf) return std::unexpected(targetHalf.error());
    auto queryIndex = query.internalIndex(queryHalf->index);
    auto targetIndex = target.internalIndex(targetHalf->index);
    if (!queryIndex) return std::unexpected("query mutation: " + queryIndex.error());
    if (!targetIndex) return std::unexpected("target mutation: " + targetIndex.error());
    if (query[*queryIndex] != queryHalf->wild) {
        return std::unexpected("encoded query wild-type base does not match the sequence");
    }
    if (target[*targetIndex] != targetHalf->wild) {
        return std::unexpected("encoded target wild-type base does not match the sequence");
    }
    if (queryHalf->mutated == queryHalf->wild || targetHalf->mutated == targetHalf->wild) {
        return std::unexpected("a compensatory candidate must mutate both sequence positions");
    }
    if (!canPair(targetHalf->mutated, queryHalf->mutated)) {
        return std::unexpected("encoded mutant bases do not form a canonical or GU pair");
    }
    return MutationCandidate{
        *queryIndex,
        *targetIndex,
        queryHalf->wild,
        targetHalf->wild,
        queryHalf->mutated,
        targetHalf->mutated};
}

auto applyMutation(
    const Sequence& query,
    const Sequence& target,
    const MutationCandidate& mutation) -> std::expected<MutationSequences, std::string> {
    if (mutation.queryIndex >= query.size() || mutation.targetIndex >= target.size()) {
        return std::unexpected("mutation position is outside its sequence");
    }
    if (query[mutation.queryIndex] != mutation.wildQuery ||
        target[mutation.targetIndex] != mutation.wildTarget) {
        return std::unexpected("mutation wild-type bases do not match the sequences");
    }
    if (mutation.mutatedQuery == mutation.wildQuery ||
        mutation.mutatedTarget == mutation.wildTarget) {
        return std::unexpected("a compensatory candidate must mutate both sequence positions");
    }
    if (!canPair(mutation.mutatedTarget, mutation.mutatedQuery)) {
        return std::unexpected("mutated bases do not form a canonical or GU pair");
    }
    MutationSequences result{query.str(), target.str(), query.str(), target.str()};
    result.mutatedQuery[mutation.queryIndex] = mutation.mutatedQuery;
    result.mutatedTarget[mutation.targetIndex] = mutation.mutatedTarget;
    return result;
}

auto shuffleMononucleotides(const std::string_view sequence, const std::uint64_t seed)
    -> std::expected<std::string, std::string> {
    auto parsed = canonicalSequence(sequence);
    if (!parsed) {
        return std::unexpected(parsed.error());
    }
    SplitMix64 random(seed);
    std::vector<char> symbols(parsed->begin(), parsed->end());
    shuffle(symbols, random);
    return std::string(symbols.begin(), symbols.end());
}

auto shuffleDinucleotides(const std::string_view sequence, const std::uint64_t seed)
    -> std::expected<std::string, std::string> {
    auto parsed = canonicalSequence(sequence);
    if (!parsed) {
        return std::unexpected(parsed.error());
    }
    if (parsed->size() < 2U) {
        return *parsed;
    }
    std::array<std::vector<char>, 4U> originalEdges;
    for (std::size_t index = 0U; index + 1U < parsed->size(); ++index) {
        originalEdges[baseIndex((*parsed)[index])].push_back((*parsed)[index + 1U]);
    }

    // Shuffling adjacency order and following a complete Hierholzer traversal
    // samples from valid Euler trails. The traversal always uses every edge for
    // this graph, whose degree conditions derive from an existing sequence.
    SplitMix64 random(seed);
    auto edges = originalEdges;
    for (auto& outgoing : edges) {
        shuffle(outgoing, random);
    }
    std::array<std::size_t, 4U> used{};
    std::vector<char> stack{parsed->front()};
    std::vector<char> reversedTrail;
    reversedTrail.reserve(parsed->size());
    while (!stack.empty()) {
        const char current = stack.back();
        auto& offset = used[baseIndex(current)];
        const auto& outgoing = edges[baseIndex(current)];
        if (offset < outgoing.size()) {
            stack.push_back(outgoing[offset++]);
        } else {
            reversedTrail.push_back(current);
            stack.pop_back();
        }
    }
    if (reversedTrail.size() != parsed->size()) {
        return std::unexpected("could not construct a complete dinucleotide-preserving trail");
    }
    std::reverse(reversedTrail.begin(), reversedTrail.end());
    if (reversedTrail.front() != parsed->front() || reversedTrail.back() != parsed->back()) {
        return std::unexpected("dinucleotide shuffle did not preserve sequence endpoints");
    }
    return std::string(reversedTrail.begin(), reversedTrail.end());
}

} // namespace intarnanew::tools
