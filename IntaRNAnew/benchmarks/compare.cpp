#include <algorithm>
#include <cerrno>
#include <charconv>
#include <chrono>
#include <csignal>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <expected>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <locale>
#include <numeric>
#include <optional>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

#include <fcntl.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

namespace {

using Clock = std::chrono::steady_clock;

struct Arguments {
    std::filesystem::path legacy;
    std::filesystem::path current;
    std::filesystem::path cases;
    std::size_t warmups{1U};
    std::size_t repetitions{7U};
    unsigned int threads{1U};
    std::size_t timeoutSeconds{300U};
    std::optional<std::size_t> expectedCases;
    bool verifyOnly{};
    bool skipWithoutLegacy{};
};

struct RunResult {
    std::chrono::nanoseconds elapsed{};
    long peakRssKiB{};
    int exitCode{};
    std::optional<int> signal;
    std::string output;
};

[[nodiscard]] auto usage() -> std::string_view {
    return
        "Usage: IntaRNAnewBenchmark [--legacy PATH] --new PATH --cases FILE|DIR\n"
        "       [--warmups N] [--repetitions N] [--threads N] [--timeout N]\n"
        "       [--verify-only] [--expect-cases N] [--skip-without-legacy]\n\n"
        "--legacy may be omitted when INTARNA_LEGACY_BIN names the executable.\n"
        "Every .parameter case is first run through both executables and compared\n"
        "byte-for-byte. A mismatching case is never timed. Timed repetitions are\n"
        "alternated to reduce order bias; stdout and stderr are discarded. Use\n"
        "--verify-only for a compatibility gate without timing.\n";
}

template <typename Integer>
[[nodiscard]] auto parseInteger(
    const std::string_view text,
    const std::string_view option,
    const bool allowZero) -> std::expected<Integer, std::string> {
    Integer value{};
    const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
    if (text.empty() || error != std::errc{} || end != text.data() + text.size() ||
        (!allowZero && value == 0)) {
        return std::unexpected(std::string(option) + " requires " +
                               (allowZero ? "a nonnegative" : "a positive") + " integer");
    }
    return value;
}

[[nodiscard]] auto parseArguments(
    const int argc,
    char** argv) -> std::expected<Arguments, std::string> {
    Arguments result;
    for (int index = 1; index < argc; ++index) {
        const std::string_view option{argv[index]};
        const auto value = [&]() -> std::expected<std::string_view, std::string> {
            if (index + 1 >= argc) {
                return std::unexpected("missing value for " + std::string(option));
            }
            return std::string_view{argv[++index]};
        };
        if (option == "--help" || option == "-h") {
            return std::unexpected(std::string{});
        }
        if (option == "--verify-only") {
            result.verifyOnly = true;
            continue;
        }
        if (option == "--skip-without-legacy") {
            result.skipWithoutLegacy = true;
            continue;
        }
        if (option == "--legacy") {
            auto parsed = value();
            if (!parsed) return std::unexpected(parsed.error());
            result.legacy = *parsed;
        } else if (option == "--new") {
            auto parsed = value();
            if (!parsed) return std::unexpected(parsed.error());
            result.current = *parsed;
        } else if (option == "--cases") {
            auto parsed = value();
            if (!parsed) return std::unexpected(parsed.error());
            result.cases = *parsed;
        } else if (option == "--warmups") {
            auto parsed = value();
            if (!parsed) return std::unexpected(parsed.error());
            auto number = parseInteger<std::size_t>(*parsed, option, true);
            if (!number) return std::unexpected(number.error());
            result.warmups = *number;
        } else if (option == "--repetitions") {
            auto parsed = value();
            if (!parsed) return std::unexpected(parsed.error());
            auto number = parseInteger<std::size_t>(*parsed, option, false);
            if (!number) return std::unexpected(number.error());
            result.repetitions = *number;
        } else if (option == "--threads") {
            auto parsed = value();
            if (!parsed) return std::unexpected(parsed.error());
            auto number = parseInteger<unsigned int>(*parsed, option, false);
            if (!number) return std::unexpected(number.error());
            result.threads = *number;
        } else if (option == "--timeout") {
            auto parsed = value();
            if (!parsed) return std::unexpected(parsed.error());
            auto number = parseInteger<std::size_t>(*parsed, option, false);
            if (!number) return std::unexpected(number.error());
            result.timeoutSeconds = *number;
        } else if (option == "--expect-cases") {
            auto parsed = value();
            if (!parsed) return std::unexpected(parsed.error());
            auto number = parseInteger<std::size_t>(*parsed, option, false);
            if (!number) return std::unexpected(number.error());
            result.expectedCases = *number;
        } else {
            return std::unexpected("unknown option '" + std::string(option) + "'");
        }
    }
    if (result.legacy.empty()) {
        if (const char* environment = std::getenv("INTARNA_LEGACY_BIN");
            environment != nullptr && *environment != '\0') {
            result.legacy = environment;
        }
    }
    if (result.current.empty() || result.cases.empty()) {
        return std::unexpected("--new and --cases are required");
    }
    if (result.legacy.empty() && !result.skipWithoutLegacy) {
        return std::unexpected("--legacy or INTARNA_LEGACY_BIN is required");
    }
    const auto executable = [](const std::filesystem::path& path) {
        std::error_code error;
        return std::filesystem::is_regular_file(path, error) &&
               ::access(path.c_str(), X_OK) == 0;
    };
    if (!result.legacy.empty() && !executable(result.legacy)) {
        return std::unexpected("legacy executable is not a runnable regular file");
    }
    if (!executable(result.current)) {
        return std::unexpected("new executable is not a runnable regular file");
    }
    return result;
}

[[nodiscard]] auto benchmarkCases(
    const std::filesystem::path& source) -> std::expected<std::vector<std::filesystem::path>, std::string> {
    std::error_code error;
    const auto absolute = [](const std::filesystem::path& path)
        -> std::expected<std::filesystem::path, std::string> {
        std::error_code pathError;
        auto result = std::filesystem::absolute(path, pathError);
        if (pathError) {
            return std::unexpected("cannot resolve benchmark case path: " + pathError.message());
        }
        return result;
    };
    if (std::filesystem::is_regular_file(source, error)) {
        if (source.extension() != ".parameter") {
            return std::unexpected("benchmark case must have a .parameter suffix");
        }
        auto resolved = absolute(source);
        if (!resolved) return std::unexpected(resolved.error());
        return std::vector<std::filesystem::path>{std::move(*resolved)};
    }
    if (!std::filesystem::is_directory(source, error)) {
        return std::unexpected("benchmark case path is neither a file nor a directory");
    }
    std::vector<std::filesystem::path> result;
    for (std::filesystem::directory_iterator iterator(source, error), end;
         !error && iterator != end; iterator.increment(error)) {
        if (iterator->is_regular_file(error) && iterator->path().extension() == ".parameter") {
            auto resolved = absolute(iterator->path());
            if (!resolved) return std::unexpected(resolved.error());
            result.push_back(std::move(*resolved));
        }
    }
    if (error) return std::unexpected("cannot enumerate benchmark cases: " + error.message());
    std::ranges::sort(result);
    if (result.empty()) return std::unexpected("no .parameter benchmark cases found");
    return result;
}

[[nodiscard]] auto temporaryCapture() -> std::expected<std::pair<int, std::filesystem::path>, std::string> {
    std::error_code error;
    const auto directory = std::filesystem::temp_directory_path(error);
    if (error) {
        return std::unexpected("cannot locate temporary directory: " + error.message());
    }
    auto pattern = (directory / "intarnanew-benchmark-output-XXXXXX").string();
    std::vector<char> writable(pattern.begin(), pattern.end());
    writable.push_back('\0');
    const int descriptor = ::mkstemp(writable.data());
    if (descriptor < 0) {
        return std::unexpected("cannot create benchmark capture file: " +
                               std::error_code(errno, std::generic_category()).message());
    }
    return std::pair{descriptor, std::filesystem::path{writable.data()}};
}

[[nodiscard]] auto readCapture(const std::filesystem::path& path)
    -> std::expected<std::string, std::string> {
    std::ifstream input(path, std::ios::binary);
    if (!input) return std::unexpected("cannot read benchmark capture file");
    std::string bytes(
        std::istreambuf_iterator<char>{input}, std::istreambuf_iterator<char>{});
    if (input.bad()) return std::unexpected("failed while reading benchmark capture file");
    return bytes;
}

[[nodiscard]] auto run(
    const std::filesystem::path& executable,
    const std::filesystem::path& parameterFile,
    const unsigned int threads,
    const std::size_t timeoutSeconds,
    const bool capture) -> std::expected<RunResult, std::string> {
    int outputDescriptor{-1};
    std::filesystem::path outputPath;
    if (capture) {
        auto temporary = temporaryCapture();
        if (!temporary) return std::unexpected(temporary.error());
        outputDescriptor = temporary->first;
        outputPath = std::move(temporary->second);
    } else {
        outputDescriptor = ::open("/dev/null", O_WRONLY);
        if (outputDescriptor < 0) return std::unexpected("cannot open /dev/null for benchmark output");
    }
    const int errorDescriptor = ::open("/dev/null", O_WRONLY);
    if (errorDescriptor < 0) {
        ::close(outputDescriptor);
        if (capture) std::filesystem::remove(outputPath);
        return std::unexpected("cannot open /dev/null for benchmark diagnostics");
    }

    std::vector<std::string> storage{
        executable.string(),
        "--parameterFile=" + parameterFile.string(),
        "--threads=" + std::to_string(threads),
        "--default-log-file=/dev/null",
        "--out=STDOUT",
    };
    std::vector<char*> childArguments;
    childArguments.reserve(storage.size() + 1U);
    for (auto& argument : storage) childArguments.push_back(argument.data());
    childArguments.push_back(nullptr);

    const auto started = Clock::now();
    const auto child = ::fork();
    if (child < 0) {
        const auto message = std::error_code(errno, std::generic_category()).message();
        ::close(outputDescriptor);
        ::close(errorDescriptor);
        if (capture) std::filesystem::remove(outputPath);
        return std::unexpected("cannot fork benchmark process: " + message);
    }
    if (child == 0) {
        if (::dup2(outputDescriptor, STDOUT_FILENO) < 0 ||
            ::dup2(errorDescriptor, STDERR_FILENO) < 0) {
            ::_exit(126);
        }
        ::close(outputDescriptor);
        ::close(errorDescriptor);
        ::execv(storage.front().c_str(), childArguments.data());
        ::_exit(127);
    }
    ::close(outputDescriptor);
    ::close(errorDescriptor);

    int status{};
    rusage resources{};
    pid_t waited{};
    const auto deadline = Clock::now() + std::chrono::seconds(timeoutSeconds);
    do {
        waited = ::wait4(child, &status, WNOHANG, &resources);
        if (waited == 0 && Clock::now() >= deadline) {
            static_cast<void>(::kill(child, SIGKILL));
            do {
                waited = ::wait4(child, &status, 0, &resources);
            } while (waited < 0 && errno == EINTR);
            break;
        }
        if (waited == 0) ::usleep(1'000U);
    } while (waited == 0 || (waited < 0 && errno == EINTR));
    const auto stopped = Clock::now();
    if (waited < 0) {
        if (capture) std::filesystem::remove(outputPath);
        return std::unexpected("cannot wait for benchmark process: " +
                               std::error_code(errno, std::generic_category()).message());
    }

    RunResult result;
    result.elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(stopped - started);
#if defined(__APPLE__)
    result.peakRssKiB = resources.ru_maxrss / 1024L;
#else
    result.peakRssKiB = resources.ru_maxrss;
#endif
    if (WIFEXITED(status)) result.exitCode = WEXITSTATUS(status);
    else if (WIFSIGNALED(status)) result.signal = WTERMSIG(status);

    if (capture) {
        auto bytes = readCapture(outputPath);
        std::error_code ignored;
        std::filesystem::remove(outputPath, ignored);
        if (!bytes) return std::unexpected(bytes.error());
        result.output = std::move(*bytes);
    }
    return result;
}

[[nodiscard]] auto failure(const RunResult& result) -> std::optional<std::string> {
    if (result.signal) return "terminated by signal " + std::to_string(*result.signal);
    if (result.exitCode != 0) return "exited with status " + std::to_string(result.exitCode);
    return std::nullopt;
}

template <typename Value>
[[nodiscard]] auto median(std::vector<Value> values) -> double {
    std::ranges::sort(values);
    const auto middle = values.size() / 2U;
    if (values.size() % 2U != 0U) return static_cast<double>(values[middle]);
    return (static_cast<double>(values[middle - 1U]) + static_cast<double>(values[middle])) / 2.0;
}

[[nodiscard]] auto medianMilliseconds(
    const std::vector<std::chrono::nanoseconds>& values) -> double {
    std::vector<long long> counts;
    counts.reserve(values.size());
    for (const auto value : values) counts.push_back(value.count());
    return median(std::move(counts)) / 1'000'000.0;
}

[[nodiscard]] auto firstDifference(
    const std::string_view left,
    const std::string_view right) noexcept -> std::size_t {
    const auto common = std::min(left.size(), right.size());
    for (std::size_t index{}; index < common; ++index) {
        if (left[index] != right[index]) return index;
    }
    return common;
}

[[nodiscard]] auto csvField(const std::string_view value) -> std::string {
    if (value.find_first_of(",\"\r\n") == std::string_view::npos) return std::string(value);
    std::string result{"\""};
    for (const char character : value) {
        if (character == '"') result.push_back('"');
        result.push_back(character);
    }
    result.push_back('"');
    return result;
}

} // namespace

auto main(const int argc, char** argv) -> int {
    auto arguments = parseArguments(argc, argv);
    if (!arguments) {
        if (!arguments.error().empty()) std::cerr << "IntaRNAnewBenchmark: " << arguments.error() << '\n';
        std::cerr << usage();
        return arguments.error().empty() ? 0 : 2;
    }
    auto cases = benchmarkCases(arguments->cases);
    if (!cases) {
        std::cerr << "IntaRNAnewBenchmark: " << cases.error() << '\n';
        return 2;
    }
    if (arguments->expectedCases && cases->size() != *arguments->expectedCases) {
        std::cerr << "IntaRNAnewBenchmark: expected " << *arguments->expectedCases
                  << " parameter cases but found " << cases->size() << '\n';
        return 2;
    }
    if (arguments->legacy.empty()) {
        std::cout << "SKIP: INTARNA_LEGACY_BIN is not set\n";
        return 0;
    }

    std::cout.imbue(std::locale::classic());
    if (!arguments->verifyOnly) {
        std::cout << "case,legacy_median_ms,new_median_ms,speedup,legacy_peak_rss_kib,new_peak_rss_kib\n";
    }
    std::vector<double> speedups;
    bool failed{};
    for (const auto& parameterFile : *cases) {
        auto legacyOracle = run(arguments->legacy, parameterFile, arguments->threads,
                                arguments->timeoutSeconds, true);
        auto newOracle = run(arguments->current, parameterFile, arguments->threads,
                             arguments->timeoutSeconds, true);
        if (!legacyOracle || !newOracle) {
            std::cerr << parameterFile.filename().string() << ": verification launch failed: "
                      << (!legacyOracle ? legacyOracle.error() : newOracle.error()) << '\n';
            failed = true;
            continue;
        }
        if (const auto reason = failure(*legacyOracle)) {
            std::cerr << parameterFile.filename().string() << ": legacy " << *reason << '\n';
            failed = true;
            continue;
        }
        if (const auto reason = failure(*newOracle)) {
            std::cerr << parameterFile.filename().string() << ": new " << *reason << '\n';
            failed = true;
            continue;
        }
        if (legacyOracle->output != newOracle->output) {
            const auto offset = firstDifference(legacyOracle->output, newOracle->output);
            std::cerr << parameterFile.filename().string()
                      << ": output mismatch at byte " << offset
                      << " (legacy " << legacyOracle->output.size()
                      << " bytes, new " << newOracle->output.size()
                      << " bytes); case was not timed\n";
            failed = true;
            continue;
        }

        if (arguments->verifyOnly) {
            std::cout << "PASS " << parameterFile.filename().string() << '\n';
            speedups.push_back(1.0);
            continue;
        }

        bool warmupFailed{};
        for (std::size_t iteration{}; iteration < arguments->warmups; ++iteration) {
            const bool newFirst = iteration % 2U == 0U;
            const auto& first = newFirst ? arguments->current : arguments->legacy;
            const auto& second = newFirst ? arguments->legacy : arguments->current;
            auto firstRun = run(first, parameterFile, arguments->threads,
                                arguments->timeoutSeconds, false);
            auto secondRun = run(second, parameterFile, arguments->threads,
                                 arguments->timeoutSeconds, false);
            if (!firstRun || !secondRun || failure(*firstRun) || failure(*secondRun)) {
                std::cerr << parameterFile.filename().string() << ": warm-up failed\n";
                failed = true;
                warmupFailed = true;
                break;
            }
        }
        if (warmupFailed) continue;

        std::vector<std::chrono::nanoseconds> legacyTimes;
        std::vector<std::chrono::nanoseconds> newTimes;
        std::vector<long> legacyMemory;
        std::vector<long> newMemory;
        legacyTimes.reserve(arguments->repetitions);
        newTimes.reserve(arguments->repetitions);
        legacyMemory.reserve(arguments->repetitions);
        newMemory.reserve(arguments->repetitions);
        bool caseFailed{};
        for (std::size_t iteration{}; iteration < arguments->repetitions; ++iteration) {
            const bool newFirst = iteration % 2U == 0U;
            auto first = run(newFirst ? arguments->current : arguments->legacy,
                             parameterFile, arguments->threads,
                             arguments->timeoutSeconds, false);
            auto second = run(newFirst ? arguments->legacy : arguments->current,
                              parameterFile, arguments->threads,
                              arguments->timeoutSeconds, false);
            if (!first || !second || failure(*first) || failure(*second)) {
                caseFailed = true;
                break;
            }
            const auto& newRun = newFirst ? *first : *second;
            const auto& legacyRun = newFirst ? *second : *first;
            newTimes.push_back(newRun.elapsed);
            legacyTimes.push_back(legacyRun.elapsed);
            newMemory.push_back(newRun.peakRssKiB);
            legacyMemory.push_back(legacyRun.peakRssKiB);
        }
        if (caseFailed) {
            std::cerr << parameterFile.filename().string() << ": timed execution failed\n";
            failed = true;
            continue;
        }

        const auto legacyMedian = medianMilliseconds(legacyTimes);
        const auto newMedian = medianMilliseconds(newTimes);
        const auto speedup = legacyMedian / newMedian;
        speedups.push_back(speedup);
        std::cout << csvField(parameterFile.stem().string()) << ','
                  << std::fixed << std::setprecision(3) << legacyMedian << ','
                  << newMedian << ',' << speedup << ','
                  << static_cast<long long>(std::llround(median(legacyMemory))) << ','
                  << static_cast<long long>(std::llround(median(newMemory))) << '\n';
    }

    if (!arguments->verifyOnly && !speedups.empty()) {
        const auto logSum = std::accumulate(
            speedups.begin(), speedups.end(), 0.0,
            [](const double sum, const double value) { return sum + std::log(value); });
        std::cout << "GEOMEAN,,," << std::fixed << std::setprecision(3)
                  << std::exp(logSum / static_cast<double>(speedups.size())) << ",,\n";
    }
    if (arguments->verifyOnly) {
        std::cout << "COMPATIBLE " << speedups.size() << '/' << cases->size() << '\n';
    }
    return failed || speedups.size() != cases->size() ? 1 : 0;
}
