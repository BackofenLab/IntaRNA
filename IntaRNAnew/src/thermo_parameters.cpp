#include "thermo_parameters.hpp"

#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <mutex>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

namespace intarnanew::detail {
namespace {

namespace fs = std::filesystem;

// Keep input work bounded independently of caller-provided files.
constexpr double referenceTemperatureKelvin = 310.15;
constexpr std::uintmax_t parameterFileSizeLimit = 16U * 1024U * 1024U;
constexpr std::size_t parameterLineSizeLimit = 64U * 1024U;
constexpr std::size_t sectionTokenLimit = 20'000U;

using RawSections = std::map<std::string, std::vector<double>, std::less<>>;

constexpr std::array requiredSections{
    std::string_view{"stack"},
    std::string_view{"stack_enthalpies"},
    std::string_view{"mismatch_internal"},
    std::string_view{"mismatch_internal_enthalpies"},
    std::string_view{"mismatch_internal_1n"},
    std::string_view{"mismatch_internal_1n_enthalpies"},
    std::string_view{"mismatch_internal_23"},
    std::string_view{"mismatch_internal_23_enthalpies"},
    std::string_view{"mismatch_exterior"},
    std::string_view{"mismatch_exterior_enthalpies"},
    std::string_view{"dangle5"},
    std::string_view{"dangle5_enthalpies"},
    std::string_view{"dangle3"},
    std::string_view{"dangle3_enthalpies"},
    std::string_view{"int11"},
    std::string_view{"int11_enthalpies"},
    std::string_view{"int21"},
    std::string_view{"int21_enthalpies"},
    std::string_view{"int22"},
    std::string_view{"int22_enthalpies"},
    std::string_view{"bulge"},
    std::string_view{"bulge_enthalpies"},
    std::string_view{"internal"},
    std::string_view{"internal_enthalpies"},
    std::string_view{"NINIO"},
    std::string_view{"Misc"},
};

[[nodiscard]] auto trim(const std::string_view value) noexcept -> std::string_view {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string_view::npos) return {};
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1U);
}

[[nodiscard]] auto isRequiredSection(const std::string_view name) noexcept -> bool {
    return std::ranges::find(requiredSections, name) != requiredSections.end();
}

[[nodiscard]] auto stripBlockComments(const std::string_view line, bool& inComment) -> std::string {
    std::string result;
    result.reserve(line.size());
    for (std::size_t index = 0U; index < line.size();) {
        if (inComment) {
            const auto end = line.find("*/", index);
            if (end == std::string_view::npos) return result;
            inComment = false;
            index = end + 2U;
            continue;
        }
        const auto begin = line.find("/*", index);
        if (begin == std::string_view::npos) {
            result.append(line.substr(index));
            break;
        }
        result.append(line.substr(index, begin - index));
        inComment = true;
        index = begin + 2U;
    }
    return result;
}

[[nodiscard]] auto parseToken(const std::string_view token, const fs::path& path,
                              const std::size_t lineNumber) -> double {
    if (token == "INF") return std::numeric_limits<double>::infinity();
    if (token == "NST") return 0.0;
    if (token == "DEF") return -50.0;

    const std::string owned(token);
    char* end{};
    errno = 0;
    const auto value = std::strtod(owned.c_str(), &end);
    if (errno == ERANGE || end != owned.c_str() + owned.size() || !std::isfinite(value)) {
        throw std::invalid_argument("invalid parameter token '" + owned + "' at " +
                                    path.string() + ':' + std::to_string(lineNumber));
    }
    return value;
}

[[nodiscard]] auto parseParameterFile(const fs::path& path) -> RawSections {
    std::error_code error;
    const auto fileSize = fs::file_size(path, error);
    if (error) {
        throw std::invalid_argument("cannot determine ViennaRNA parameter file size: " + path.string());
    }
    if (fileSize > parameterFileSizeLimit) {
        throw std::invalid_argument("ViennaRNA parameter file exceeds the 16 MiB parser limit: " +
                                    path.string());
    }

    std::ifstream input(path);
    if (!input) throw std::invalid_argument("cannot open ViennaRNA parameter file: " + path.string());

    RawSections sections;
    std::string activeSection;
    std::string line;
    std::size_t lineNumber{};
    bool inComment{};
    bool versionSeen{};
    while (std::getline(input, line)) {
        ++lineNumber;
        if (line.size() > parameterLineSizeLimit) {
            throw std::invalid_argument("overlong line in ViennaRNA parameter file: " + path.string());
        }
        const auto cleaned = stripBlockComments(line, inComment);
        const auto view = trim(cleaned);
        if (view.empty()) continue;
        if (view.starts_with("##")) {
            if (!view.starts_with("## RNAfold parameter file v2")) {
                throw std::invalid_argument("unsupported ViennaRNA parameter file header in " + path.string());
            }
            versionSeen = true;
            activeSection.clear();
            continue;
        }
        if (view.front() == '#') {
            auto name = trim(view.substr(1U));
            if (const auto whitespace = name.find_first_of(" \t"); whitespace != std::string_view::npos) {
                name = name.substr(0U, whitespace);
            }
            activeSection.clear();
            if (isRequiredSection(name)) {
                auto [iterator, inserted] = sections.try_emplace(std::string(name));
                static_cast<void>(iterator);
                if (!inserted) {
                    throw std::invalid_argument("duplicate ViennaRNA parameter section # " +
                                                std::string(name) + " in " + path.string());
                }
                activeSection = std::string(name);
            }
            continue;
        }
        if (activeSection.empty()) continue;

        std::istringstream tokens{std::string(view)};
        std::string token;
        auto& values = sections.at(activeSection);
        while (tokens >> token) {
            if (values.size() >= sectionTokenLimit) {
                throw std::invalid_argument("too many values in ViennaRNA parameter section # " +
                                            activeSection + " in " + path.string());
            }
            values.push_back(parseToken(token, path, lineNumber));
        }
    }
    if (!input.eof()) throw std::invalid_argument("failed while reading ViennaRNA parameter file: " + path.string());
    if (inComment) throw std::invalid_argument("unterminated block comment in ViennaRNA parameter file: " + path.string());
    if (!versionSeen) throw std::invalid_argument("not a ViennaRNA v2 parameter file: " + path.string());
    return sections;
}

void appendEnvironmentDirectories(std::vector<fs::path>& directories, const char* variable) {
    const auto* value = std::getenv(variable);
    if (value == nullptr || *value == '\0') return;
#ifdef _WIN32
    constexpr char separator = ';';
#else
    constexpr char separator = ':';
#endif
    const std::string_view list(value);
    std::size_t begin{};
    while (begin <= list.size() && directories.size() < 64U) {
        const auto end = list.find(separator, begin);
        const auto item = trim(list.substr(begin, end == std::string_view::npos
            ? std::string_view::npos : end - begin));
        if (!item.empty()) directories.emplace_back(item);
        if (end == std::string_view::npos) break;
        begin = end + 1U;
    }
}

[[nodiscard]] auto canonicalRegularFile(const fs::path& path) -> std::optional<fs::path> {
    std::error_code error;
    if (!fs::is_regular_file(path, error) || error) return std::nullopt;
    auto canonical = fs::weakly_canonical(path, error);
    return error ? std::optional<fs::path>{path} : std::optional<fs::path>{std::move(canonical)};
}

[[nodiscard]] auto resolveParameterFile(const std::string_view parameterSet) -> fs::path {
    const auto requested = parameterSet.empty() ? std::string_view{"Turner04"} : parameterSet;
    if (auto direct = canonicalRegularFile(fs::path(requested))) return *direct;

    std::string_view fileName;
    if (requested == "Turner04" || requested == "rna_turner2004.par") {
        fileName = "rna_turner2004.par";
    } else if (requested == "Turner99" || requested == "rna_turner1999.par") {
        fileName = "rna_turner1999.par";
    } else if (requested == "Andronescu07" || requested == "rna_andronescu2007.par") {
        fileName = "rna_andronescu2007.par";
    } else if (requested.find('/') != std::string_view::npos ||
               requested.find('\\') != std::string_view::npos || requested.ends_with(".par")) {
        throw std::invalid_argument("ViennaRNA parameter file not found: " + std::string(requested));
    } else {
        throw std::invalid_argument("unknown ViennaRNA parameter set '" + std::string(requested) +
                                    "' (expected Turner04, Turner99, Andronescu07, or a .par path)");
    }

    std::vector<fs::path> directories;
    appendEnvironmentDirectories(directories, "INTARNANEW_PARAMETER_DIR");
    appendEnvironmentDirectories(directories, "VIENNA_RNA_DATAPATH");
    appendEnvironmentDirectories(directories, "VRNA_DATAPATH");
    if (const auto* prefix = std::getenv("CONDA_PREFIX"); prefix != nullptr && *prefix != '\0') {
        directories.emplace_back(fs::path(prefix) / "share" / "ViennaRNA");
    }
#ifdef _WIN32
    if (const auto* programData = std::getenv("PROGRAMDATA"); programData != nullptr && *programData != '\0') {
        directories.emplace_back(fs::path(programData) / "ViennaRNA");
    }
#else
    directories.emplace_back("/usr/local/share/ViennaRNA");
    directories.emplace_back("/usr/share/ViennaRNA");
#endif

    std::error_code error;
    auto current = fs::current_path(error);
    for (std::size_t depth = 0U; !error && !current.empty() && depth < 8U; ++depth) {
        directories.emplace_back(current / "share" / "ViennaRNA");
        directories.emplace_back(current / ".conda-env" / "share" / "ViennaRNA");
        directories.emplace_back(current / "intaRNA_legacy" / ".conda-env" / "share" / "ViennaRNA");
        const auto parent = current.parent_path();
        if (parent == current) break;
        current = parent;
    }

    std::set<fs::path> visited;
    for (const auto& directory : directories) {
        if (directory.empty() || !visited.insert(directory).second) continue;
        if (auto resolved = canonicalRegularFile(directory / fileName)) return *resolved;
    }
    throw std::invalid_argument("ViennaRNA parameter set '" + std::string(requested) +
                                "' was not found; set INTARNANEW_PARAMETER_DIR or pass an explicit .par path");
}

[[nodiscard]] auto section(const RawSections& sections, const std::string_view name,
                           const std::size_t expected, const fs::path& path) -> const std::vector<double>& {
    const auto iterator = sections.find(name);
    if (iterator == sections.end()) {
        throw std::invalid_argument("missing ViennaRNA parameter section # " + std::string(name) +
                                    " in " + path.string());
    }
    if (iterator->second.size() != expected) {
        throw std::invalid_argument("ViennaRNA parameter section # " + std::string(name) + " in " +
                                    path.string() + " has " + std::to_string(iterator->second.size()) +
                                    " values; expected " + std::to_string(expected));
    }
    return iterator->second;
}

[[nodiscard]] auto scaledValue(const double freeEnergy37, const double enthalpy,
                               const double temperatureKelvin, const bool quantizeToCentikcal = true) -> Energy {
    if (!std::isfinite(freeEnergy37) || !std::isfinite(enthalpy)) return infinity;
    if (std::abs(temperatureKelvin - referenceTemperatureKelvin) < 1e-9) {
        return freeEnergy37 / 100.0;
    }
    const auto value = enthalpy - (enthalpy - freeEnergy37) *
        temperatureKelvin / referenceTemperatureKelvin;
    return quantizeToCentikcal ? std::trunc(value) / 100.0 : value / 100.0;
}

[[nodiscard]] auto scaledSection(const RawSections& sections, const std::string_view freeEnergyName,
                                 const std::string_view enthalpyName, const std::size_t expected,
                                 const fs::path& path, const double temperatureKelvin) -> std::vector<Energy> {
    const auto& freeEnergies = section(sections, freeEnergyName, expected, path);
    const auto& enthalpies = section(sections, enthalpyName, expected, path);
    std::vector<Energy> result;
    result.reserve(expected);
    for (std::size_t index = 0U; index < expected; ++index) {
        result.push_back(scaledValue(freeEnergies[index], enthalpies[index], temperatureKelvin));
    }
    return result;
}

[[nodiscard]] auto buildParameters(const fs::path& path, const double temperatureCelsius)
    -> std::shared_ptr<const NearestNeighborParameters> {
    const auto temperatureKelvin = temperatureCelsius + 273.15;
    const auto raw = parseParameterFile(path);
    auto result = std::make_shared<NearestNeighborParameters>();
    result->stack = scaledSection(raw, "stack", "stack_enthalpies", 7U * 7U, path, temperatureKelvin);
    result->mismatchInternal = scaledSection(raw, "mismatch_internal",
        "mismatch_internal_enthalpies", 7U * 5U * 5U, path, temperatureKelvin);
    result->mismatchInternal1n = scaledSection(raw, "mismatch_internal_1n",
        "mismatch_internal_1n_enthalpies", 7U * 5U * 5U, path, temperatureKelvin);
    result->mismatchInternal23 = scaledSection(raw, "mismatch_internal_23",
        "mismatch_internal_23_enthalpies", 7U * 5U * 5U, path, temperatureKelvin);
    result->mismatchExterior = scaledSection(raw, "mismatch_exterior",
        "mismatch_exterior_enthalpies", 7U * 5U * 5U, path, temperatureKelvin);
    result->dangle5 = scaledSection(raw, "dangle5", "dangle5_enthalpies",
        7U * 5U, path, temperatureKelvin);
    result->dangle3 = scaledSection(raw, "dangle3", "dangle3_enthalpies",
        7U * 5U, path, temperatureKelvin);
    result->int11 = scaledSection(raw, "int11", "int11_enthalpies",
        7U * 7U * 5U * 5U, path, temperatureKelvin);
    result->int21 = scaledSection(raw, "int21", "int21_enthalpies",
        7U * 7U * 5U * 5U * 5U, path, temperatureKelvin);
    result->int22 = scaledSection(raw, "int22", "int22_enthalpies",
        6U * 6U * 4U * 4U * 4U * 4U, path, temperatureKelvin);
    result->bulge = scaledSection(raw, "bulge", "bulge_enthalpies", 31U, path, temperatureKelvin);
    result->internal = scaledSection(raw, "internal", "internal_enthalpies", 31U, path, temperatureKelvin);

    const auto& ninio = section(raw, "NINIO", 3U, path);
    result->ninioSlope = scaledValue(ninio[0], ninio[1], temperatureKelvin);
    result->ninioMaximum = ninio[2] / 100.0;

    const auto miscIterator = raw.find("Misc");
    if (miscIterator == raw.end() || (miscIterator->second.size() != 4U && miscIterator->second.size() != 6U)) {
        throw std::invalid_argument("ViennaRNA parameter section # Misc in " + path.string() +
                                    " must contain 4 or 6 values");
    }
    const auto& misc = miscIterator->second;
    result->duplexInit = scaledValue(misc[0], misc[1], temperatureKelvin);
    result->terminalAu = scaledValue(misc[2], misc[3], temperatureKelvin);
    // ViennaRNA v2 permits the historical LXC pair to be omitted; its specified default is 107.856/0.
    result->logarithmicLoopSlope = misc.size() == 6U
        ? scaledValue(misc[4], misc[5], temperatureKelvin, false)
        : scaledValue(107.856, 0.0, temperatureKelvin, false);
    return result;
}

struct CacheKey {
    fs::path path;
    double temperatureCelsius{};
    friend auto operator<=>(const CacheKey&, const CacheKey&) = default;
};

} // namespace

auto loadNearestNeighborParameters(const double temperatureCelsius,
                                   const std::string_view parameterSet)
    -> std::shared_ptr<const NearestNeighborParameters> {
    if (!std::isfinite(temperatureCelsius) || temperatureCelsius <= -273.15) {
        throw std::invalid_argument("energy temperature must be finite and above absolute zero");
    }
    const CacheKey key{resolveParameterFile(parameterSet), temperatureCelsius};
    static std::mutex cacheMutex;
    static std::map<CacheKey, std::weak_ptr<const NearestNeighborParameters>> cache;
    {
        const std::scoped_lock lock(cacheMutex);
        if (const auto found = cache.find(key); found != cache.end()) {
            if (auto existing = found->second.lock()) return existing;
            cache.erase(found);
        }
    }

    auto loaded = buildParameters(key.path, temperatureCelsius);
    {
        const std::scoped_lock lock(cacheMutex);
        std::erase_if(cache, [](const auto& entry) { return entry.second.expired(); });
        if (cache.size() < 16U) cache[key] = loaded;
    }
    return loaded;
}

} // namespace intarnanew::detail
