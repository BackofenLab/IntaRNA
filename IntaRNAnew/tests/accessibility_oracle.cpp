#include "intarnanew/accessibility.hpp"
#include "intarnanew/compression.hpp"
#include "intarnanew/sequence.hpp"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace {

using intarnanew::Index;
using intarnanew::Interval;
using intarnanew::NativeAccessibility;
using intarnanew::Sequence;
using intarnanew::SideConfig;
using intarnanew::gasConstantKcal;

struct Edge {
    Index left{};
    Index right{};
    double energy{};
};

auto failures = 0;

void check(const bool condition, const std::string_view message) {
    if (!condition) {
        std::cerr << "FAILED: " << message << '\n';
        ++failures;
    }
}

[[nodiscard]] auto close(const double left, const double right, const double tolerance = 1e-11) -> bool {
    const auto scale = std::max({1.0, std::abs(left), std::abs(right)});
    return std::abs(left - right) <= tolerance * scale;
}

[[nodiscard]] auto pairEnergy(const char left, const char right) -> double {
    const auto a = static_cast<char>(std::toupper(static_cast<unsigned char>(left)));
    const auto b = static_cast<char>(std::toupper(static_cast<unsigned char>(right)));
    if ((a == 'G' && b == 'C') || (a == 'C' && b == 'G')) return -2.4;
    if ((a == 'A' && b == 'U') || (a == 'U' && b == 'A')) return -1.1;
    if ((a == 'G' && b == 'U') || (a == 'U' && b == 'G')) return -0.5;
    return -0.8;
}

[[nodiscard]] auto crosses(const Edge& left, const Edge& right) -> bool {
    return (left.left < right.left && right.left < left.right && left.right < right.right) ||
           (right.left < left.left && left.left < right.right && right.right < left.right);
}

[[nodiscard]] auto brutePartition(
    const Sequence& sequence,
    const std::string_view constraints,
    const Index maxPairSpan,
    const std::optional<Interval> forcedUnpaired = std::nullopt) -> double {
    std::vector<Edge> edges;
    for (Index left = 0; left < sequence.size(); ++left) {
        for (Index right = left + 4U; right < sequence.size(); ++right) {
            if (right - left > maxPairSpan ||
                constraints[left] == 'x' || constraints[left] == 'b' ||
                constraints[right] == 'x' || constraints[right] == 'b' ||
                (forcedUnpaired && (forcedUnpaired->contains(left) || forcedUnpaired->contains(right))) ||
                !intarnanew::canPair(sequence[left], sequence[right])) {
                continue;
            }
            edges.push_back({left, right, pairEnergy(sequence[left], sequence[right])});
        }
    }

    const auto rt = gasConstantKcal * (37.0 + 273.15);
    std::vector<bool> paired(sequence.size(), false);
    std::vector<Edge> selected;
    double partition{};

    const std::function<void(Index, double)> enumerate = [&](const Index edgeIndex, const double energy) {
        if (edgeIndex == edges.size()) {
            for (Index position = 0; position < sequence.size(); ++position) {
                if (constraints[position] == 'p' && !paired[position]) return;
            }
            partition += std::exp(-energy / rt);
            return;
        }

        enumerate(edgeIndex + 1U, energy);

        const auto edge = edges[edgeIndex];
        if (paired[edge.left] || paired[edge.right] ||
            std::ranges::any_of(selected, [&](const Edge chosen) { return crosses(edge, chosen); })) {
            return;
        }
        paired[edge.left] = true;
        paired[edge.right] = true;
        selected.push_back(edge);
        enumerate(edgeIndex + 1U, energy + edge.energy);
        selected.pop_back();
        paired[edge.left] = false;
        paired[edge.right] = false;
    };
    enumerate(0U, 0.0);
    return partition;
}

void compareEveryInterval(
    const std::string_view bases,
    const std::string_view constraints,
    const Index maxPairSpan) {
    const Sequence sequence("oracle", std::string(bases));
    SideConfig config;
    config.accessibilitySpan = maxPairSpan;
    config.accessibilityConstraint = std::string(constraints);
    const NativeAccessibility provider(sequence, config, 37.0);

    const auto denominator = brutePartition(sequence, constraints, maxPairSpan);
    check(denominator > 0.0, "enumerated constrained ensemble is nonempty");
    for (Index begin = 0; begin < sequence.size(); ++begin) {
        for (Index end = begin; end < sequence.size(); ++end) {
            const Interval interval{begin, end};
            const bool containsBlocked = std::string_view(constraints).substr(
                begin, end - begin + 1U).find('b') != std::string_view::npos;
            if (containsBlocked) {
                check(provider.unpairedProbability(interval) == 0.0,
                      "blocked interval is excluded from interaction");
                continue;
            }
            const auto expected = brutePartition(
                sequence, constraints, maxPairSpan, interval) / denominator;
            const auto observed = provider.unpairedProbability(interval);
            check(close(observed, expected), "joint interval probability matches exhaustive enumeration");
            check(close(provider.positionUnpairedProbability(begin),
                        brutePartition(sequence, constraints, maxPairSpan, Interval{begin, begin}) /
                            denominator),
                  "single-position marginal matches exhaustive enumeration");
            const auto opening = provider.openingEnergy(interval);
            if (expected == 0.0) {
                check(std::isinf(opening), "zero joint probability has infinite opening energy");
            } else {
                const auto rt = gasConstantKcal * (37.0 + 273.15);
                check(close(opening, -rt * std::log(expected)),
                      "opening energy is -RT log of the joint probability");
            }
        }
    }
}

} // namespace

auto main() -> int {
    using namespace intarnanew;

    check(NativeAccessibility::modelName == "native-noncrossing-pair-v1",
          "native model is explicitly identified");

    compareEveryInterval("GCAAAUGC", "........", 8U);
    compareEveryInterval("GCAAAUGC", "p.......", 8U);
    compareEveryInterval("GCAAAUGC", "x.......", 8U);
    compareEveryInterval("GCAAAUGC", ".......b", 8U);
    compareEveryInterval("GCAAAUGC", "........", 4U);

    {
        const Sequence sequence("correlation", "GAAAC");
        SideConfig config;
        config.accessibilitySpan = sequence.size();
        NativeAccessibility provider(sequence, config, 37.0);
        const auto singleLeft = provider.unpairedProbability({0U, 0U});
        const auto singleRight = provider.unpairedProbability({4U, 4U});
        const auto joint = provider.unpairedProbability({0U, 4U});
        const auto expected = brutePartition(sequence, ".....", sequence.size(), Interval{0U, 4U}) /
                              brutePartition(sequence, ".....", sequence.size());
        check(close(joint, expected), "correlated endpoint joint probability is exact");
        check(!close(joint, singleLeft * singleRight, 1e-6),
              "joint probability is not fabricated from marginal products");
    }

    {
        const Sequence sequence("forced-pair", "GAAAC");
        SideConfig config;
        config.accessibilitySpan = sequence.size();
        config.accessibilityConstraint = "p...p";
        NativeAccessibility provider(sequence, config, 37.0);
        check(provider.positionUnpairedProbability(0U) == 0.0,
              "p forces its position to be paired");
        check(provider.positionUnpairedProbability(4U) == 0.0,
              "both forced-pair endpoints remain paired");
        check(provider.positionUnpairedProbability(2U) == 1.0,
              "unpairable hairpin interior remains unpaired");
    }

    {
        const Sequence sequence("ranges", "GAAAC", 10);
        SideConfig config;
        config.accessibilitySpan = sequence.size();
        config.accessibilityConstraint = "x:10-10,b:14-14";
        NativeAccessibility provider(sequence, config, 37.0);
        check(provider.positionUnpairedProbability(0U) == 1.0,
              "x range forces a position unpaired");
        check(provider.blocked(4U), "b range marks a position blocked");
        check(provider.positionUnpairedProbability(4U) == 0.0,
              "blocked position is excluded from interaction");
        check(std::isinf(provider.openingEnergy({4U, 4U})),
              "blocked position has infinite interaction opening energy");
    }

    {
        const Sequence sequence("negative-ranges", "GAAAC", -5);
        SideConfig config;
        config.accessibilitySpan = sequence.size();
        config.accessibilityConstraint = "x:-5--3,b:-1--1";
        NativeAccessibility provider(sequence, config, 37.0);
        check(provider.positionUnpairedProbability(0U) == 1.0,
              "signed coordinate ranges are parsed without confusing the separator");
        check(provider.blocked(4U), "negative blocked range maps to the final position");
    }

    {
        const Sequence sequence("temperature", "GAAAC");
        SideConfig config;
        bool rejected{};
        try {
            const NativeAccessibility provider(sequence, config, -273.15);
            static_cast<void>(provider);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        check(rejected, "nonphysical absolute temperature is rejected");
    }

    {
        const Sequence sequence("impossible", "AAAAA");
        SideConfig config;
        config.accessibilityConstraint = "p....";
        bool rejected{};
        try {
            const NativeAccessibility provider(sequence, config, 37.0);
            static_cast<void>(provider);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        check(rejected, "infeasible forced-pair constraints are rejected");
    }

    {
        const auto tablePath = std::filesystem::temp_directory_path() /
                               "intarnanew-accessibility-oracle-table.txt";
        const auto gzipPath = std::filesystem::temp_directory_path() /
                              "intarnanew-accessibility-oracle-table.txt.gz";
        const auto magicPath = std::filesystem::temp_directory_path() /
                              "intarnanew-accessibility-oracle-table.bin";
        const auto falseGzipPath = std::filesystem::temp_directory_path() /
                                   "intarnanew-accessibility-oracle-false.txt.gz";
        const auto corruptGzipPath = std::filesystem::temp_directory_path() /
                                     "intarnanew-accessibility-oracle-corrupt.txt.gz";
        const std::string tableText{"1 0.5\n2 0.5\n3 0.5\n4 0.5\n5 0.5\n"};
        {
            std::ofstream table(tablePath);
            table << tableText;
        }
        const Sequence sequence("table", "GAAAC");
        SideConfig config;
        config.accessibilityFile = tablePath.string();
        TableAccessibility provider(sequence, config, 37.0, true);
        check(close(provider.unpairedProbability({0U, 0U}), 0.5),
              "table-backed single-position probability is read");
        bool rejected{};
        try {
            static_cast<void>(provider.unpairedProbability({0U, 1U}));
        } catch (const std::runtime_error&) {
            rejected = true;
        }
        check(rejected,
              "missing table joint probability is not replaced by a marginal product");

        const auto compressed = gzipCompress(tableText);
        check(compressed.has_value(), "accessibility table gzip fixture is encoded");
        if (compressed) {
            {
                std::ofstream gzipFile(gzipPath, std::ios::binary);
                gzipFile.write(compressed->data(), static_cast<std::streamsize>(compressed->size()));
            }
            config.accessibilityFile = gzipPath.string();
            TableAccessibility gzipProvider(sequence, config, 37.0, true);
            check(close(gzipProvider.positionUnpairedProbability(2U), 0.5),
                  "gzip accessibility table is decoded through the native compression API");

            {
                std::ofstream magicFile(magicPath, std::ios::binary);
                magicFile.write(compressed->data(), static_cast<std::streamsize>(compressed->size()));
            }
            config.accessibilityFile = magicPath.string();
            TableAccessibility magicProvider(sequence, config, 37.0, true);
            check(close(magicProvider.positionUnpairedProbability(4U), 0.5),
                  "gzip accessibility table is detected by signature without a .gz suffix");

            auto corrupt = *compressed;
            corrupt[corrupt.size() - 8U] = static_cast<char>(
                static_cast<unsigned char>(corrupt[corrupt.size() - 8U]) ^ 1U);
            {
                std::ofstream corruptFile(corruptGzipPath, std::ios::binary);
                corruptFile.write(corrupt.data(), static_cast<std::streamsize>(corrupt.size()));
            }
            config.accessibilityFile = corruptGzipPath.string();
            rejected = false;
            try {
                const TableAccessibility corruptProvider(sequence, config, 37.0, true);
                static_cast<void>(corruptProvider);
            } catch (const std::invalid_argument& error) {
                rejected = std::string_view(error.what()).find("CRC32") != std::string_view::npos;
            }
            check(rejected, "corrupt gzip accessibility table reports its checksum failure");
        }

        {
            std::ofstream falseGzip(falseGzipPath, std::ios::binary);
            falseGzip << tableText;
        }
        config.accessibilityFile = falseGzipPath.string();
        rejected = false;
        try {
            const TableAccessibility falseGzipProvider(sequence, config, 37.0, true);
            static_cast<void>(falseGzipProvider);
        } catch (const std::invalid_argument& error) {
            rejected = std::string_view(error.what()).find("no gzip signature") != std::string_view::npos;
        }
        check(rejected, ".gz accessibility input without gzip magic has a clear diagnostic");

        std::error_code removalError;
        for (const auto& path : {tablePath, gzipPath, magicPath, falseGzipPath, corruptGzipPath}) {
            std::filesystem::remove(path, removalError);
            check(!removalError, "temporary accessibility table fixture is removed");
            removalError.clear();
        }
    }

    {
        std::string bases;
        bases.reserve(120U);
        for (Index index = 0; index < 120U; ++index) {
            bases.push_back(index % 2U == 0U ? 'G' : 'C');
        }
        const Sequence sequence("stability", bases);
        SideConfig config;
        config.accessibilitySpan = 80U;
        NativeAccessibility provider(sequence, config, 37.0);
        const auto probability = provider.unpairedProbability({45U, 74U});
        check(std::isfinite(probability) && probability >= 0.0 && probability <= 1.0,
              "log-space recurrence remains finite on a high-weight ensemble");
        const auto repeat = provider.unpairedProbability({45U, 74U});
        check(repeat == probability, "cached probability is deterministic");
    }

    if (failures == 0) {
        std::cout << "All native accessibility oracle tests passed.\n";
    }
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
