#pragma once

#include "intarnanew/config.hpp"

#include <cstddef>
#include <expected>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace intarnanew {

struct OutputGroupKey {
    std::size_t targetIndex{};
    std::size_t queryIndex{};
    std::size_t targetRegionIndex{};
    std::size_t queryRegionIndex{};

    [[nodiscard]] friend constexpr auto operator==(
        const OutputGroupKey&, const OutputGroupKey&) noexcept -> bool = default;
};

enum class OutputPartKind {
    primary,
    pairAuxiliary,
    queryAccessibility,
    targetAccessibility,
};

struct OutputPart {
    OutputPartKind kind{OutputPartKind::primary};
    std::size_t descriptorIndex{};
    std::size_t groupIndex{};
    std::size_t sequenceIndex{};
};

struct OutputPublication {
    std::string destination;
    std::vector<OutputPart> parts;
};

struct OutputPlan {
    std::vector<OutputPublication> publications;
};

// Resolves every output destination before prediction starts. Parts belonging
// to one descriptor and one stream are deliberately collected into one writer;
// independent descriptors may never resolve to the same stream or real path.
[[nodiscard]] auto planOutputs(
    const Config& config,
    std::size_t targetCount,
    std::size_t queryCount,
    std::span<const OutputGroupKey> groups) -> std::expected<OutputPlan, std::string>;

[[nodiscard]] auto isAuxiliaryOutput(std::string_view descriptor) noexcept -> bool;
[[nodiscard]] auto auxiliaryOutputDestination(std::string_view descriptor) -> std::string;

struct OutputArtifact {
    std::string destination;
    std::string content;
};

// Stages every regular-file artifact first, then commits all files as one
// rollback-capable batch. Standard streams are written only after file commit;
// streams themselves cannot be rolled back by an operating system.
[[nodiscard]] auto publishOutputs(std::span<const OutputArtifact> artifacts)
    -> std::expected<void, std::string>;

} // namespace intarnanew
