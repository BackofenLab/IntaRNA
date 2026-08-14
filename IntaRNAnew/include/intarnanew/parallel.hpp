#pragma once

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <expected>
#include <functional>
#include <limits>
#include <stop_token>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace intarnanew {

struct ParallelFailure {
    std::size_t taskIndex{std::numeric_limits<std::size_t>::max()};
    std::string message;
};

// Executes indices in increasing claim order. Once one callback fails, the
// shared stop source prevents workers from claiming further work. Callbacks
// already in flight may finish; failures are reported by the lowest task index
// so diagnostics do not depend on thread completion order.
template <class Function>
[[nodiscard]] auto runParallelIndexed(
    const std::size_t taskCount,
    const std::size_t requestedWorkers,
    Function&& function) -> std::expected<void, ParallelFailure> {
    if (taskCount == 0U) return {};

    const auto workerCount = std::min(
        taskCount, std::max<std::size_t>(1U, requestedWorkers));
    std::atomic_size_t nextTask{};
    std::stop_source cancellation;
    const auto cancellationToken = cancellation.get_token();
    std::vector<std::string> failures(taskCount);
    std::vector<std::jthread> workers;
    workers.reserve(workerCount);

    try {
        for (std::size_t workerIndex = 0U; workerIndex < workerCount; ++workerIndex) {
            workers.emplace_back([&, workerIndex](const std::stop_token threadToken) {
                while (!cancellationToken.stop_requested() &&
                       !threadToken.stop_requested()) {
                    const auto taskIndex = nextTask.fetch_add(1U, std::memory_order_relaxed);
                    if (taskIndex >= taskCount) return;
                    if (cancellationToken.stop_requested() || threadToken.stop_requested()) return;
                    try {
                        std::invoke(function, workerIndex, taskIndex, cancellationToken);
                    } catch (const std::exception& exception) {
                        failures[taskIndex] = exception.what();
                        if (failures[taskIndex].empty()) {
                            failures[taskIndex] = "exception with no diagnostic";
                        }
                        cancellation.request_stop();
                        return;
                    } catch (...) {
                        failures[taskIndex] = "unknown exception";
                        cancellation.request_stop();
                        return;
                    }
                }
            });
        }
    } catch (const std::exception& exception) {
        cancellation.request_stop();
        return std::unexpected(ParallelFailure{
            std::numeric_limits<std::size_t>::max(),
            "failed to start worker threads: " + std::string(exception.what()),
        });
    }

    // Explicit join avoids jthread destruction requesting an otherwise
    // successful worker to stop before it has exhausted the task queue.
    for (auto& worker : workers) worker.join();

    for (std::size_t taskIndex = 0U; taskIndex < failures.size(); ++taskIndex) {
        if (!failures[taskIndex].empty()) {
            return std::unexpected(ParallelFailure{taskIndex, std::move(failures[taskIndex])});
        }
    }
    return {};
}

} // namespace intarnanew
