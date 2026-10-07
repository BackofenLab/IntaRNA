#include "intarnanew/output_plan.hpp"

#include "intarnanew/compression.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <ranges>
#include <set>
#include <system_error>
#include <tuple>
#include <utility>

namespace intarnanew {
namespace {

enum class AuxiliaryScope { pair, query, target };
enum class StreamKind { none, standardOutput, standardError };

struct DescriptorInfo {
    std::size_t index{};
    std::string kind;
    std::string destination;
    AuxiliaryScope scope{AuxiliaryScope::pair};
};

struct DestinationIdentity {
    StreamKind stream{StreamKind::none};
    std::filesystem::path path;
};

struct StagedFile {
    std::filesystem::path destination;
    std::filesystem::path temporary;
    std::filesystem::path backup;
    bool hadOriginal{};
    bool committed{};
};

[[nodiscard]] auto lowerAscii(const std::string_view value) -> std::string {
    std::string result(value);
    std::ranges::transform(result, result.begin(), [](const unsigned char character) {
        return static_cast<char>(std::tolower(character));
    });
    return result;
}

[[nodiscard]] auto descriptorKind(const std::string_view descriptor) -> std::string {
    return lowerAscii(descriptor.substr(0U, descriptor.find(':')));
}

[[nodiscard]] auto streamKind(const std::string_view destination) noexcept -> StreamKind {
    if (destination.empty() || destination == "-") return StreamKind::standardOutput;
    const auto normalized = lowerAscii(destination);
    if (normalized == "stdout") return StreamKind::standardOutput;
    if (normalized == "stderr") return StreamKind::standardError;
    return StreamKind::none;
}

[[nodiscard]] auto normalizedStreamName(const StreamKind stream) -> std::string {
    return stream == StreamKind::standardError ? "STDERR" : "STDOUT";
}

[[nodiscard]] auto auxiliaryScope(const std::string_view kind) -> AuxiliaryScope {
    if (kind == "qacc" || kind == "qpu") return AuxiliaryScope::query;
    if (kind == "tacc" || kind == "tpu") return AuxiliaryScope::target;
    return AuxiliaryScope::pair;
}

[[nodiscard]] auto withSuffix(
    const std::string_view destination,
    const std::string_view suffix) -> std::string {
    if (suffix.empty() || streamKind(destination) != StreamKind::none) {
        return std::string(destination);
    }
    const std::filesystem::path path(destination);
    const auto extension = path.extension();
    const auto filename = path.stem().string() + std::string(suffix) + extension.string();
    return (path.parent_path() / filename).string();
}

[[nodiscard]] auto destinationIdentity(const std::string_view destination)
    -> std::expected<DestinationIdentity, std::string> {
    const auto stream = streamKind(destination);
    if (stream != StreamKind::none) return DestinationIdentity{stream, {}};

    std::error_code error;
    auto absolute = std::filesystem::absolute(std::filesystem::path(destination), error);
    if (error) {
        return std::unexpected(
            "cannot resolve output destination '" + std::string(destination) + "': " +
            error.message());
    }
    auto resolved = std::filesystem::weakly_canonical(absolute, error);
    if (error) resolved = absolute.lexically_normal();

    auto status = std::filesystem::symlink_status(resolved, error);
    if (error == std::make_error_code(std::errc::no_such_file_or_directory)) {
        error.clear();
        status = std::filesystem::file_status{std::filesystem::file_type::not_found};
    } else if (error) {
        return std::unexpected(
            "cannot inspect output destination '" + std::string(destination) + "': " +
            error.message());
    }
    if (status.type() != std::filesystem::file_type::not_found &&
        !std::filesystem::is_regular_file(resolved, error)) {
        return std::unexpected(
            "output destination '" + std::string(destination) + "' is not a regular file");
    }
    if (error) {
        return std::unexpected(
            "cannot inspect output destination '" + std::string(destination) + "': " +
            error.message());
    }

    auto parent = resolved.parent_path();
    if (parent.empty()) parent = std::filesystem::current_path(error);
    if (error || !std::filesystem::is_directory(parent, error)) {
        return std::unexpected(
            "parent directory for output destination '" + std::string(destination) +
            "' does not exist or is not a directory");
    }
    return DestinationIdentity{StreamKind::none, std::move(resolved)};
}

[[nodiscard]] auto sameDestination(
    const DestinationIdentity& first,
    const DestinationIdentity& second) -> bool {
    if (first.stream != StreamKind::none || second.stream != StreamKind::none) {
        return first.stream != StreamKind::none && first.stream == second.stream;
    }
    if (first.path == second.path) return true;
    std::error_code error;
    const auto equivalent = std::filesystem::equivalent(first.path, second.path, error);
    return !error && equivalent;
}

[[nodiscard]] auto validatePublications(const std::vector<OutputPublication>& publications)
    -> std::expected<void, std::string> {
    std::vector<DestinationIdentity> identities;
    identities.reserve(publications.size());
    for (std::size_t index = 0U; index < publications.size(); ++index) {
        auto identity = destinationIdentity(publications[index].destination);
        if (!identity) return std::unexpected(identity.error());
        for (std::size_t previous = 0U; previous < identities.size(); ++previous) {
            if (!sameDestination(*identity, identities[previous])) continue;
            const auto stream = identity->stream;
            if (stream != StreamKind::none) {
                return std::unexpected(
                    "multiple output descriptors write " + normalizedStreamName(stream));
            }
            return std::unexpected(
                "output destinations '" + publications[previous].destination + "' and '" +
                publications[index].destination + "' resolve to the same file");
        }
        identities.push_back(std::move(*identity));
    }
    return {};
}

[[nodiscard]] auto uniqueSibling(
    const std::filesystem::path& destination,
    const std::string_view marker) -> std::filesystem::path {
    static std::atomic_size_t sequence{};
    const auto timestamp = std::chrono::steady_clock::now().time_since_epoch().count();
    for (;;) {
        auto candidate = destination;
        candidate += "." + std::string(marker) + "." + std::to_string(timestamp) + "." +
                     std::to_string(sequence.fetch_add(1U, std::memory_order_relaxed));
        std::error_code error;
        if (!std::filesystem::exists(candidate, error)) return candidate;
    }
}

void removeIfPresent(const std::filesystem::path& path) noexcept {
    if (path.empty()) return;
    std::error_code ignored;
    std::filesystem::remove(path, ignored);
}

void cleanStaging(std::vector<StagedFile>& files) noexcept {
    for (auto& file : files) {
        removeIfPresent(file.temporary);
        if (!file.committed) removeIfPresent(file.backup);
    }
}

[[nodiscard]] auto rollback(std::vector<StagedFile>& files, const std::string& reason)
    -> std::expected<void, std::string> {
    std::string rollbackError;
    for (auto iterator = files.rbegin(); iterator != files.rend(); ++iterator) {
        auto& file = *iterator;
        if (file.committed) {
            removeIfPresent(file.destination);
            if (file.hadOriginal) {
                std::error_code error;
                std::filesystem::rename(file.backup, file.destination, error);
                if (error && rollbackError.empty()) {
                    rollbackError = " (also failed to restore '" + file.destination.string() +
                                    "': " + error.message() + ")";
                }
            }
        } else if (file.hadOriginal && !file.backup.empty()) {
            std::error_code error;
            std::filesystem::rename(file.backup, file.destination, error);
            if (error && rollbackError.empty()) {
                rollbackError = " (also failed to restore '" + file.destination.string() +
                                "': " + error.message() + ")";
            }
        }
        removeIfPresent(file.temporary);
    }
    return std::unexpected(reason + rollbackError);
}

} // namespace

auto isAuxiliaryOutput(const std::string_view descriptor) noexcept -> bool {
    const auto kind = descriptorKind(descriptor);
    return kind == "qmine" || kind == "tmine" || kind == "qspotprob" ||
           kind == "tspotprob" || kind == "qacc" || kind == "tacc" ||
           kind == "qpu" || kind == "tpu" || kind == "pmine" || kind == "spotprob";
}

auto auxiliaryOutputDestination(const std::string_view descriptor) -> std::string {
    const auto first = descriptor.find(':');
    if (first == std::string_view::npos) return "STDOUT";
    if (descriptorKind(descriptor) == "spotprob") {
        const auto second = descriptor.find(':', first + 1U);
        return second == std::string_view::npos
            ? std::string(descriptor.substr(first + 1U))
            : std::string(descriptor.substr(second + 1U));
    }
    return std::string(descriptor.substr(first + 1U));
}

auto planOutputs(
    const Config& config,
    const std::size_t targetCount,
    const std::size_t queryCount,
    const std::span<const OutputGroupKey> groups) -> std::expected<OutputPlan, std::string> {
    if (targetCount == 0U || queryCount == 0U) {
        return std::unexpected("output planning requires at least one target and one query");
    }
    for (const auto& group : groups) {
        if (group.targetIndex >= targetCount || group.queryIndex >= queryCount) {
            return std::unexpected("output group references a sequence outside the selected set");
        }
    }

    std::vector<std::size_t> primaryIndices;
    std::vector<DescriptorInfo> auxiliary;
    for (std::size_t index = 0U; index < config.output.destinations.size(); ++index) {
        const auto& descriptor = config.output.destinations[index];
        if (!isAuxiliaryOutput(descriptor)) {
            primaryIndices.push_back(index);
            continue;
        }
        const auto kind = descriptorKind(descriptor);
        auxiliary.push_back({
            index, kind, auxiliaryOutputDestination(descriptor), auxiliaryScope(kind),
        });
    }
    if (primaryIndices.size() > 1U) {
        return std::unexpected("primary output destination was specified more than once");
    }

    OutputPlan plan;
    OutputPublication primary;
    primary.destination = primaryIndices.empty()
        ? "STDOUT" : config.output.destinations[primaryIndices.front()];
    primary.parts.reserve(groups.size());
    for (std::size_t groupIndex = 0U; groupIndex < groups.size(); ++groupIndex) {
        primary.parts.push_back({OutputPartKind::primary, 0U, groupIndex, 0U});
    }
    plan.publications.push_back(std::move(primary));

    std::set<std::pair<std::size_t, std::size_t>> distinctPairs;
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> groupsPerPair;
    for (const auto& group : groups) {
        const auto pair = std::pair{group.targetIndex, group.queryIndex};
        distinctPairs.insert(pair);
        ++groupsPerPair[pair];
    }
    const auto multiplePairs = distinctPairs.size() > 1U;

    for (const auto& descriptor : auxiliary) {
        const auto stream = streamKind(descriptor.destination);
        if (descriptor.scope == AuxiliaryScope::pair) {
            if (stream != StreamKind::none) {
                OutputPublication publication{normalizedStreamName(stream), {}};
                publication.parts.reserve(groups.size());
                for (std::size_t groupIndex = 0U; groupIndex < groups.size(); ++groupIndex) {
                    publication.parts.push_back({
                        OutputPartKind::pairAuxiliary, descriptor.index, groupIndex, 0U,
                    });
                }
                plan.publications.push_back(std::move(publication));
                continue;
            }
            for (std::size_t groupIndex = 0U; groupIndex < groups.size(); ++groupIndex) {
                const auto& group = groups[groupIndex];
                std::string suffix;
                if (multiplePairs) {
                    suffix += "-t" + std::to_string(group.targetIndex + 1U) +
                              "q" + std::to_string(group.queryIndex + 1U);
                }
                const auto pair = std::pair{group.targetIndex, group.queryIndex};
                if (groupsPerPair[pair] > 1U) {
                    suffix += "-rT" + std::to_string(group.targetRegionIndex + 1U) +
                              "Q" + std::to_string(group.queryRegionIndex + 1U);
                }
                plan.publications.push_back({
                    withSuffix(descriptor.destination, suffix),
                    {{OutputPartKind::pairAuxiliary, descriptor.index, groupIndex, 0U}},
                });
            }
            continue;
        }

        const auto sequenceCount = descriptor.scope == AuxiliaryScope::query
            ? queryCount : targetCount;
        const auto partKind = descriptor.scope == AuxiliaryScope::query
            ? OutputPartKind::queryAccessibility : OutputPartKind::targetAccessibility;
        if (stream != StreamKind::none) {
            OutputPublication publication{normalizedStreamName(stream), {}};
            publication.parts.reserve(sequenceCount);
            for (std::size_t sequenceIndex = 0U; sequenceIndex < sequenceCount; ++sequenceIndex) {
                publication.parts.push_back({
                    partKind, descriptor.index, 0U, sequenceIndex,
                });
            }
            plan.publications.push_back(std::move(publication));
            continue;
        }
        for (std::size_t sequenceIndex = 0U; sequenceIndex < sequenceCount; ++sequenceIndex) {
            const auto suffix = sequenceCount > 1U
                ? "-s" + std::to_string(sequenceIndex + 1U) : std::string{};
            plan.publications.push_back({
                withSuffix(descriptor.destination, suffix),
                {{partKind, descriptor.index, 0U, sequenceIndex}},
            });
        }
    }

    if (auto status = validatePublications(plan.publications); !status) {
        return std::unexpected(status.error());
    }
    return plan;
}

auto publishOutputs(const std::span<const OutputArtifact> artifacts)
    -> std::expected<void, std::string> {
    std::vector<OutputPublication> publications;
    publications.reserve(artifacts.size());
    for (const auto& artifact : artifacts) {
        publications.push_back({artifact.destination, {}});
    }
    if (auto status = validatePublications(publications); !status) {
        return std::unexpected(status.error());
    }

    std::vector<StagedFile> staged;
    staged.reserve(artifacts.size());
    for (const auto& artifact : artifacts) {
        if (streamKind(artifact.destination) != StreamKind::none) continue;
        const std::filesystem::path destination(artifact.destination);
        std::string compressed;
        std::string_view payload = artifact.content;
        if (lowerAscii(artifact.destination).ends_with(".gz")) {
            auto encoded = gzipCompress(artifact.content);
            if (!encoded) {
                cleanStaging(staged);
                return std::unexpected(
                    "cannot gzip output file '" + artifact.destination + "': " + encoded.error());
            }
            compressed = std::move(*encoded);
            payload = compressed;
        }
        if (payload.size() > static_cast<std::size_t>(
                std::numeric_limits<std::streamsize>::max())) {
            cleanStaging(staged);
            return std::unexpected("output file '" + artifact.destination + "' is too large to write");
        }

        StagedFile file;
        file.destination = destination;
        file.temporary = uniqueSibling(destination, "intarnanew-stage");
        {
            std::ofstream output(file.temporary, std::ios::binary | std::ios::trunc);
            if (!output) {
                cleanStaging(staged);
                return std::unexpected(
                    "cannot create temporary output for '" + artifact.destination + "'");
            }
            output.write(payload.data(), static_cast<std::streamsize>(payload.size()));
            output.close();
            if (!output) {
                removeIfPresent(file.temporary);
                cleanStaging(staged);
                return std::unexpected("failed to stage output file '" + artifact.destination + "'");
            }
        }
        staged.push_back(std::move(file));
    }

    for (auto& file : staged) {
        std::error_code error;
        auto status = std::filesystem::symlink_status(file.destination, error);
        if (error == std::make_error_code(std::errc::no_such_file_or_directory)) {
            error.clear();
            status = std::filesystem::file_status{std::filesystem::file_type::not_found};
        } else if (error) {
            return rollback(staged,
                "cannot inspect output file '" + file.destination.string() + "': " + error.message());
        }
        if (status.type() != std::filesystem::file_type::not_found) {
            if (!std::filesystem::is_regular_file(file.destination, error) || error) {
                return rollback(staged,
                    "cannot replace non-regular output '" + file.destination.string() + "'");
            }
            file.hadOriginal = true;
            file.backup = uniqueSibling(file.destination, "intarnanew-backup");
            std::filesystem::rename(file.destination, file.backup, error);
            if (error) {
                return rollback(staged,
                    "cannot preserve output file '" + file.destination.string() + "': " +
                    error.message());
            }
        }
        error.clear();
        std::filesystem::rename(file.temporary, file.destination, error);
        if (error) {
            return rollback(staged,
                "cannot commit output file '" + file.destination.string() + "': " + error.message());
        }
        file.temporary.clear();
        file.committed = true;
    }

    for (auto& file : staged) removeIfPresent(file.backup);

    for (const auto& artifact : artifacts) {
        const auto stream = streamKind(artifact.destination);
        if (stream == StreamKind::standardOutput) {
            std::cout << artifact.content;
            if (!std::cout) return std::unexpected("failed to write standard output");
        } else if (stream == StreamKind::standardError) {
            std::cerr << artifact.content;
            if (!std::cerr) return std::unexpected("failed to write standard error");
        }
    }
    return {};
}

} // namespace intarnanew
