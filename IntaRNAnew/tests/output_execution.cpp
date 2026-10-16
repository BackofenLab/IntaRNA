#include "intarnanew/compression.hpp"
#include "intarnanew/output_plan.hpp"
#include "intarnanew/parallel.hpp"

#include <atomic>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

int failures{};

void check(const bool condition, const std::string_view description) {
    if (!condition) {
        std::cerr << "FAILED: " << description << '\n';
        ++failures;
    }
}

[[nodiscard]] auto readAll(const std::filesystem::path& path) -> std::string {
    std::ifstream input(path, std::ios::binary);
    return {std::istreambuf_iterator<char>{input}, std::istreambuf_iterator<char>{}};
}

[[nodiscard]] auto publication(
    const intarnanew::OutputPlan& plan,
    const std::string_view destination) -> const intarnanew::OutputPublication* {
    for (const auto& item : plan.publications) {
        if (item.destination == destination) return &item;
    }
    return nullptr;
}

} // namespace

auto main() -> int {
    using namespace intarnanew;

    {
        Config config;
        config.output.destinations = {
            "primary.csv",
            "qMinE:pair.csv.gz",
            "qAcc:qacc.csv.gz",
            "tPu:tpu.txt",
        };
        const std::vector<OutputGroupKey> groups{
            {0U, 0U, 0U, 0U}, {0U, 1U, 0U, 0U},
            {1U, 0U, 0U, 0U}, {1U, 1U, 0U, 0U},
        };
        const auto plan = planOutputs(config, 2U, 2U, groups);
        check(plan.has_value(), "multi-sequence output plan is valid");
        if (plan) {
            const auto* primary = publication(*plan, "primary.csv");
            check(primary != nullptr && primary->parts.size() == 4U,
                  "primary output combines every pair in deterministic order");
            check(publication(*plan, "pair.csv-t1q1.gz") != nullptr &&
                      publication(*plan, "pair.csv-t1q2.gz") != nullptr &&
                      publication(*plan, "pair.csv-t2q1.gz") != nullptr &&
                      publication(*plan, "pair.csv-t2q2.gz") != nullptr,
                  "pair suffixes precede the gzip extension and identify target/query");
            const auto* query1 = publication(*plan, "qacc.csv-s1.gz");
            const auto* query2 = publication(*plan, "qacc.csv-s2.gz");
            check(query1 != nullptr && query2 != nullptr &&
                      query1->parts.size() == 1U && query2->parts.size() == 1U &&
                      query1->parts.front().sequenceIndex == 0U &&
                      query2->parts.front().sequenceIndex == 1U,
                  "query accessibility is emitted once per query with legacy -s suffixes");
            check(publication(*plan, "tpu-s1.txt") != nullptr &&
                      publication(*plan, "tpu-s2.txt") != nullptr,
                  "target accessibility is emitted once per target");
        }
    }

    {
        Config config;
        config.output.destinations = {"primary.csv", "qMinE:pair.csv"};
        const std::vector<OutputGroupKey> groups{
            {0U, 0U, 0U, 0U}, {0U, 0U, 1U, 0U},
        };
        const auto plan = planOutputs(config, 1U, 1U, groups);
        check(plan && publication(*plan, "pair-rT1Q1.csv") != nullptr &&
                  publication(*plan, "pair-rT2Q1.csv") != nullptr,
              "per-region pair artifacts receive deterministic collision-free suffixes");
    }

    {
        Config streams;
        streams.output.destinations = {"STDOUT", "qMinE:STDOUT"};
        const std::vector<OutputGroupKey> groups{{0U, 0U, 0U, 0U}};
        const auto plan = planOutputs(streams, 1U, 1U, groups);
        check(!plan && plan.error().find("multiple output descriptors") != std::string::npos,
              "conflicting standard-stream writers are rejected before prediction");

        Config files;
        files.output.destinations = {"same.csv", "qMinE:same.csv"};
        const auto collision = planOutputs(files, 1U, 1U, groups);
        check(!collision && collision.error().find("same file") != std::string::npos,
              "colliding real-file destinations are rejected before prediction");

        Config oneQuery;
        oneQuery.output.destinations = {"main.csv", "qAcc:qacc.csv"};
        const std::vector<OutputGroupKey> twoPairs{
            {0U, 0U, 0U, 0U}, {1U, 0U, 0U, 0U},
        };
        const auto oneQueryPlan = planOutputs(oneQuery, 2U, 1U, twoPairs);
        const auto* query = oneQueryPlan ? publication(*oneQueryPlan, "qacc.csv") : nullptr;
        check(query != nullptr && query->parts.size() == 1U,
              "one query accessibility artifact is not duplicated across target pairs");
    }

    {
        const auto unique = std::to_string(
            std::chrono::steady_clock::now().time_since_epoch().count());
        const auto directory = std::filesystem::temp_directory_path() /
                               ("intarnanew-publication-test-" + unique);
        std::error_code error;
        std::filesystem::create_directory(directory, error);
        check(!error, "batch-publication test directory is created");

        const auto preserved = directory / "preserved.txt";
        {
            std::ofstream output(preserved);
            output << "old\n";
        }
        const std::vector<OutputArtifact> invalidBatch{
            {preserved.string(), "new\n"},
            {(directory / "missing" / "later.txt").string(), "later\n"},
        };
        const auto rejected = publishOutputs(invalidBatch);
        check(!rejected && readAll(preserved) == "old\n",
              "a later invalid destination cannot partially publish an earlier file");

        const auto plain = directory / "plain.txt";
        const auto compressed = directory / "compressed.txt.gz";
        const std::vector<OutputArtifact> validBatch{
            {plain.string(), "plain payload\n"},
            {compressed.string(), "compressed payload\n"},
        };
        const auto published = publishOutputs(validBatch);
        const auto decoded = published ? gzipDecompress(readAll(compressed))
                                       : std::expected<std::string, std::string>{
                                             std::unexpected("publication failed")};
        check(published && readAll(plain) == "plain payload\n" && decoded &&
                  *decoded == "compressed payload\n",
              "batch publication commits plain and gzip artifacts together");

        const std::vector<OutputArtifact> duplicate{
            {plain.string(), "first"}, {plain.string(), "second"},
        };
        check(!publishOutputs(duplicate) && readAll(plain) == "plain payload\n",
              "duplicate batch destinations fail without changing existing content");

        std::filesystem::remove_all(directory, error);
        check(!error, "batch-publication test directory is removed");
    }

    {
        std::vector<std::size_t> values(64U);
        const auto success = runParallelIndexed(
            values.size(), 8U,
            [&](const std::size_t, const std::size_t task, const std::stop_token token) {
                if (token.stop_requested()) throw std::runtime_error("unexpected cancellation");
                values[task] = task + 1U;
            });
        bool complete = success.has_value();
        for (std::size_t index = 0U; index < values.size(); ++index) {
            complete = complete && values[index] == index + 1U;
        }
        check(complete, "parallel indexed execution completes every task exactly once");

        std::atomic_size_t claimed{};
        const auto failed = runParallelIndexed(
            100'000U, 8U,
            [&](const std::size_t, const std::size_t task, const std::stop_token) {
                claimed.fetch_add(1U, std::memory_order_relaxed);
                if (task == 3U) throw std::runtime_error("intentional worker failure");
            });
        check(!failed && failed.error().taskIndex == 3U &&
                  failed.error().message == "intentional worker failure",
              "parallel failures retain a deterministic task index and diagnostic");
        check(claimed.load(std::memory_order_relaxed) < 100'000U,
              "the first worker failure stops new task claims cooperatively");
    }

    if (failures == 0) {
        std::cout << "All output planning and parallel execution tests passed.\n";
    }
    return failures == 0 ? 0 : 1;
}
