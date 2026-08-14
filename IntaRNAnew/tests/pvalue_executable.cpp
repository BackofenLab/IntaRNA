#include <algorithm>
#include <array>
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
#include <iostream>
#include <iterator>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

#include <fcntl.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

namespace {

struct Result {
    int exitCode{};
    std::string output;
    std::string diagnostics;
};

[[nodiscard]] auto temporaryFile(const std::string_view label)
    -> std::expected<std::pair<int, std::filesystem::path>, std::string> {
    std::error_code error;
    const auto directory = std::filesystem::temp_directory_path(error);
    if (error) return std::unexpected("cannot locate temporary directory: " + error.message());
    auto pattern = (directory / ("intarnanew-pvalue-" + std::string(label) + "-XXXXXX")).string();
    std::vector<char> bytes(pattern.begin(), pattern.end());
    bytes.push_back('\0');
    const int descriptor = ::mkstemp(bytes.data());
    if (descriptor < 0) {
        return std::unexpected("cannot create temporary capture: " +
                               std::error_code(errno, std::generic_category()).message());
    }
    return std::pair{descriptor, std::filesystem::path{bytes.data()}};
}

[[nodiscard]] auto readFile(const std::filesystem::path& path)
    -> std::expected<std::string, std::string> {
    std::ifstream input(path, std::ios::binary);
    if (!input) return std::unexpected("cannot open captured output");
    std::string bytes(std::istreambuf_iterator<char>{input}, std::istreambuf_iterator<char>{});
    if (input.bad()) return std::unexpected("cannot read captured output");
    return bytes;
}

[[nodiscard]] auto temporaryInput(
    const std::string_view label,
    const std::string_view contents) -> std::expected<std::filesystem::path, std::string> {
    auto file = temporaryFile(label);
    if (!file) return std::unexpected(file.error());
    ::close(file->first);
    std::ofstream output(file->second, std::ios::binary | std::ios::trunc);
    output.write(contents.data(), static_cast<std::streamsize>(contents.size()));
    output.close();
    if (!output) {
        std::error_code ignored;
        std::filesystem::remove(file->second, ignored);
        return std::unexpected("cannot write temporary p-value input");
    }
    return std::move(file->second);
}

[[nodiscard]] auto run(
    const std::filesystem::path& executable,
    std::vector<std::string> arguments) -> std::expected<Result, std::string> {
    auto output = temporaryFile("stdout");
    if (!output) return std::unexpected(output.error());
    auto diagnostics = temporaryFile("stderr");
    if (!diagnostics) {
        ::close(output->first);
        std::error_code ignored;
        std::filesystem::remove(output->second, ignored);
        return std::unexpected(diagnostics.error());
    }
    arguments.insert(arguments.begin(), executable.string());
    std::vector<char*> argv;
    argv.reserve(arguments.size() + 1U);
    for (auto& argument : arguments) argv.push_back(argument.data());
    argv.push_back(nullptr);

    const auto child = ::fork();
    if (child < 0) {
        ::close(output->first);
        ::close(diagnostics->first);
        std::error_code ignored;
        std::filesystem::remove(output->second, ignored);
        std::filesystem::remove(diagnostics->second, ignored);
        return std::unexpected("cannot fork p-value executable test");
    }
    if (child == 0) {
        if (::dup2(output->first, STDOUT_FILENO) < 0 ||
            ::dup2(diagnostics->first, STDERR_FILENO) < 0) {
            ::_exit(126);
        }
        ::close(output->first);
        ::close(diagnostics->first);
        ::execv(arguments.front().c_str(), argv.data());
        ::_exit(127);
    }
    ::close(output->first);
    ::close(diagnostics->first);
    int status{};
    pid_t waited{};
    bool timedOut{};
    const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(60);
    do {
        waited = ::waitpid(child, &status, WNOHANG);
        if (waited == 0 && std::chrono::steady_clock::now() >= deadline) {
            timedOut = true;
            static_cast<void>(::kill(child, SIGKILL));
            do {
                waited = ::waitpid(child, &status, 0);
            } while (waited < 0 && errno == EINTR);
            break;
        }
        if (waited == 0) ::usleep(1'000U);
    } while (waited == 0 || (waited < 0 && errno == EINTR));
    auto capturedOutput = readFile(output->second);
    auto capturedDiagnostics = readFile(diagnostics->second);
    std::error_code ignored;
    std::filesystem::remove(output->second, ignored);
    std::filesystem::remove(diagnostics->second, ignored);
    if (waited < 0) return std::unexpected("cannot wait for p-value executable test");
    if (timedOut) return std::unexpected("p-value executable timed out after 60 seconds");
    if (!capturedOutput) return std::unexpected(capturedOutput.error());
    if (!capturedDiagnostics) return std::unexpected(capturedDiagnostics.error());
    const int exitCode = WIFEXITED(status) ? WEXITSTATUS(status) : 128;
    return Result{exitCode, std::move(*capturedOutput), std::move(*capturedDiagnostics)};
}

[[nodiscard]] auto lines(const std::string_view text) -> std::size_t {
    return static_cast<std::size_t>(std::ranges::count(text, '\n'));
}

[[nodiscard]] auto finiteProbability(const std::string_view text) -> bool {
    double value{};
    const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
    const auto trailing = std::string_view{end, static_cast<std::size_t>(text.data() + text.size() - end)};
    return error == std::errc{} && std::isfinite(value) && value >= 0.0 && value <= 1.0 &&
           (trailing.empty() || trailing == "\n");
}

} // namespace

auto main(const int argc, char** argv) -> int {
    if (argc != 2) {
        std::cerr << "usage: IntaRNAnewPvalueExecutableTests PATH\n";
        return 2;
    }
    const std::filesystem::path executable{argv[1]};
    int failures{};
    const auto require = [&](const bool condition, const std::string_view description) {
        if (!condition) {
            std::cerr << "FAIL: " << description << '\n';
            ++failures;
        }
    };

    auto help = run(executable, {"--help"});
    require(help && help->exitCode == 0 && help->output.starts_with("intarnanew-pvalue") &&
                help->output.find("--cardinality") != std::string::npos,
            "help succeeds and identifies the executable");

    auto parameterFile = temporaryInput("parameters",
        "energy=B\nacc=N\nnoSeed=true\nmode=M\nmodel=S\n"
        "intLenMax=12\noutNumber=1\n");
    auto queryFasta = temporaryInput("query-fasta",
        ">query-first\nAACCGGUUAUCG\n>query-ignored\nAAAA\n");
    auto targetFasta = temporaryInput("target-fasta",
        ">target-first\nGCUAACCGGAUU\n>target-ignored\nUUUU\n");
    const auto removeInputs = [&] {
        std::error_code ignored;
        if (parameterFile) std::filesystem::remove(*parameterFile, ignored);
        if (queryFasta) std::filesystem::remove(*queryFasta, ignored);
        if (targetFasta) std::filesystem::remove(*targetFasta, ignored);
    };
    if (!parameterFile || !queryFasta || !targetFasta) {
        std::cerr << "cannot create p-value contract input: "
                  << (!parameterFile ? parameterFile.error()
                      : (!queryFasta ? queryFasta.error() : targetFasta.error())) << '\n';
        removeInputs();
        return 2;
    }
    const auto parameterOption = "--parameterFile=" + parameterFile->string();

    const std::vector<std::string> common{
        "--query=AACCGGUUAUCG", "--target=GCUAACCGGAUU", "--samples=32",
        "--shuffle-mode=b", "--randSeed=42", "--output=scores", parameterOption};
    auto oneArguments = common;
    oneArguments.push_back("--threads=1");
    auto fourArguments = common;
    fourArguments.push_back("--threads=4");
    auto oneThread = run(executable, std::move(oneArguments));
    auto fourThreads = run(executable, std::move(fourArguments));
    require(oneThread && fourThreads && oneThread->exitCode == 0 && fourThreads->exitCode == 0,
            "score-mode executions succeed");
    require(oneThread && fourThreads && oneThread->output == fourThreads->output,
            "score output is byte-identical across one and four threads");
    require(oneThread && lines(oneThread->output) == 32U,
            "score output contains the requested sample count");

    auto probability = run(executable, {
        "--query=AACCGGUUAUCG", "--target=GCUAACCGGAUU", "--samples=80",
        "--shuffle-mode=b", "--randSeed=42", "--threads=4",
        "--distribution=gauss", "--output=pvalue", parameterOption});
    require(probability && probability->exitCode == 0 && finiteProbability(probability->output),
            "p-value mode emits one finite probability in [0,1]");

    auto fasta = run(executable, {
        "--query=" + queryFasta->string(), "--target=" + targetFasta->string(),
        "--cardinality=8", "--shuffle-mode=b", "--randSeed=42", "--threads=2",
        "--output=scores", parameterOption});
    auto literal = run(executable, {
        "--query=AACCGGUUAUCG", "--target=GCUAACCGGAUU", "--samples=8",
        "--shuffle-mode=b", "--randSeed=42", "--threads=2",
        "--output=scores", parameterOption});
    require(fasta && literal && fasta->exitCode == 0 && literal->exitCode == 0 &&
                fasta->output == literal->output && lines(fasta->output) == 8U,
            "FASTA first-record input and --cardinality match literal --samples output");

    auto malformed = run(executable, {
        "--query=A", "--target=U", "--samples=0"});
    require(malformed && malformed->exitCode == 2 &&
                malformed->diagnostics.find("samples must be positive") != std::string::npos,
            "malformed sample count fails with a diagnostic");

    removeInputs();

    if (failures == 0) std::cout << "All p-value executable contract tests passed.\n";
    return failures == 0 ? 0 : 1;
}
