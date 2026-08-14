#include "folding_parameters.hpp"

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
#include <system_error>
#include <utility>

namespace intarnanew::folding_detail {
namespace {

namespace fs = std::filesystem;
constexpr double referenceTemperatureKelvin = 310.15;
constexpr std::uintmax_t maximumFileSize = 16U * 1024U * 1024U;
constexpr std::size_t maximumLineSize = 64U * 1024U;
constexpr std::size_t maximumSectionTokens = 20'000U;

using NumericSections = std::map<std::string, std::vector<double>, std::less<>>;
using StringSections = std::map<std::string, std::vector<std::string>, std::less<>>;

struct RawParameters {
    NumericSections numeric;
    StringSections textual;
};

constexpr std::array numericSectionNames{
    std::string_view{"stack"}, std::string_view{"stack_enthalpies"},
    std::string_view{"mismatch_hairpin"}, std::string_view{"mismatch_hairpin_enthalpies"},
    std::string_view{"mismatch_internal"}, std::string_view{"mismatch_internal_enthalpies"},
    std::string_view{"mismatch_internal_1n"}, std::string_view{"mismatch_internal_1n_enthalpies"},
    std::string_view{"mismatch_internal_23"}, std::string_view{"mismatch_internal_23_enthalpies"},
    std::string_view{"mismatch_multi"}, std::string_view{"mismatch_multi_enthalpies"},
    std::string_view{"mismatch_exterior"}, std::string_view{"mismatch_exterior_enthalpies"},
    std::string_view{"dangle5"}, std::string_view{"dangle5_enthalpies"},
    std::string_view{"dangle3"}, std::string_view{"dangle3_enthalpies"},
    std::string_view{"int11"}, std::string_view{"int11_enthalpies"},
    std::string_view{"int21"}, std::string_view{"int21_enthalpies"},
    std::string_view{"int22"}, std::string_view{"int22_enthalpies"},
    std::string_view{"hairpin"}, std::string_view{"hairpin_enthalpies"},
    std::string_view{"bulge"}, std::string_view{"bulge_enthalpies"},
    std::string_view{"internal"}, std::string_view{"internal_enthalpies"},
    std::string_view{"ML_params"}, std::string_view{"NINIO"}, std::string_view{"Misc"},
};

[[nodiscard]] auto trim(const std::string_view text) noexcept -> std::string_view {
    const auto first = text.find_first_not_of(" \t\r\n");
    if (first == std::string_view::npos) return {};
    const auto last = text.find_last_not_of(" \t\r\n");
    return text.substr(first, last - first + 1U);
}

[[nodiscard]] auto isNumericSection(const std::string_view name) noexcept -> bool {
    return std::ranges::find(numericSectionNames, name) != numericSectionNames.end();
}

[[nodiscard]] auto isSpecialLoopSection(const std::string_view name) noexcept -> bool {
    return name == "Triloops" || name == "Tetraloops" || name == "Hexaloops";
}

[[nodiscard]] auto stripBlockComments(const std::string_view line, bool& inComment) -> std::string {
    std::string result;
    result.reserve(line.size());
    for (std::size_t cursor{}; cursor < line.size();) {
        if (inComment) {
            const auto end = line.find("*/", cursor);
            if (end == std::string_view::npos) return result;
            inComment = false;
            cursor = end + 2U;
            continue;
        }
        const auto begin = line.find("/*", cursor);
        if (begin == std::string_view::npos) {
            result.append(line.substr(cursor));
            break;
        }
        result.append(line.substr(cursor, begin - cursor));
        inComment = true;
        cursor = begin + 2U;
    }
    return result;
}

[[nodiscard]] auto parseNumber(const std::string_view token, const fs::path& path,
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

[[nodiscard]] auto parseFile(const fs::path& path) -> RawParameters {
    std::error_code error;
    const auto fileSize = fs::file_size(path, error);
    if (error || fileSize > maximumFileSize) {
        throw std::invalid_argument(error ? "cannot determine parameter file size: " + path.string()
                                          : "parameter file exceeds the 16 MiB limit: " + path.string());
    }
    std::ifstream input(path);
    if (!input) throw std::invalid_argument("cannot open ViennaRNA parameter file: " + path.string());

    RawParameters result;
    std::string active;
    std::string line;
    std::size_t lineNumber{};
    bool inComment{};
    bool versionSeen{};
    while (std::getline(input, line)) {
        ++lineNumber;
        if (line.size() > maximumLineSize) {
            throw std::invalid_argument("overlong ViennaRNA parameter line in " + path.string());
        }
        const auto cleaned = stripBlockComments(line, inComment);
        const auto view = trim(cleaned);
        if (view.empty()) continue;
        if (view.starts_with("##")) {
            if (!view.starts_with("## RNAfold parameter file v2")) {
                throw std::invalid_argument("unsupported ViennaRNA parameter header in " + path.string());
            }
            versionSeen = true;
            active.clear();
            continue;
        }
        if (view.front() == '#') {
            auto name = trim(view.substr(1U));
            if (const auto blank = name.find_first_of(" \t"); blank != std::string_view::npos) {
                name = name.substr(0U, blank);
            }
            active.clear();
            if (isNumericSection(name)) {
                if (!result.numeric.try_emplace(std::string(name)).second) {
                    throw std::invalid_argument("duplicate ViennaRNA section # " + std::string(name));
                }
                active = std::string(name);
            } else if (isSpecialLoopSection(name)) {
                if (!result.textual.try_emplace(std::string(name)).second) {
                    throw std::invalid_argument("duplicate ViennaRNA section # " + std::string(name));
                }
                active = std::string(name);
            }
            continue;
        }
        if (active.empty()) continue;
        if (isSpecialLoopSection(active)) {
            result.textual.at(active).emplace_back(view);
            continue;
        }
        auto& values = result.numeric.at(active);
        std::istringstream tokens{std::string(view)};
        std::string token;
        while (tokens >> token) {
            if (values.size() == maximumSectionTokens) {
                throw std::invalid_argument("too many values in ViennaRNA section # " + active);
            }
            values.push_back(parseNumber(token, path, lineNumber));
        }
    }
    if (!input.eof() || inComment || !versionSeen) {
        throw std::invalid_argument("malformed ViennaRNA parameter file: " + path.string());
    }
    return result;
}

void appendEnvironmentPaths(std::vector<fs::path>& paths, const char* variable) {
    const auto* value = std::getenv(variable);
    if (value == nullptr || *value == '\0') return;
#ifdef _WIN32
    constexpr char separator = ';';
#else
    constexpr char separator = ':';
#endif
    const std::string_view encoded(value);
    std::size_t begin{};
    while (begin <= encoded.size() && paths.size() < 64U) {
        const auto end = encoded.find(separator, begin);
        const auto item = trim(encoded.substr(begin, end == std::string_view::npos
            ? std::string_view::npos : end - begin));
        if (!item.empty()) paths.emplace_back(item);
        if (end == std::string_view::npos) break;
        begin = end + 1U;
    }
}

[[nodiscard]] auto regularFile(const fs::path& path) -> std::optional<fs::path> {
    std::error_code error;
    if (!fs::is_regular_file(path, error) || error) return std::nullopt;
    auto canonical = fs::weakly_canonical(path, error);
    return error ? std::optional<fs::path>{path} : std::optional<fs::path>{std::move(canonical)};
}

[[nodiscard]] auto resolveParameterFile(const std::string_view parameterSet) -> fs::path {
    const auto requested = parameterSet.empty() ? std::string_view{"Turner04"} : parameterSet;
    if (auto direct = regularFile(fs::path(requested))) return *direct;
    std::string_view fileName;
    if (requested == "Turner04" || requested == "rna_turner2004.par") {
        fileName = "rna_turner2004.par";
    } else if (requested == "Turner99" || requested == "rna_turner1999.par") {
        fileName = "rna_turner1999.par";
    } else if (requested == "Andronescu07" || requested == "rna_andronescu2007.par") {
        fileName = "rna_andronescu2007.par";
    } else {
        throw std::invalid_argument("unknown or missing ViennaRNA parameter set '" +
                                    std::string(requested) + "'");
    }

    std::vector<fs::path> directories;
    appendEnvironmentPaths(directories, "INTARNANEW_PARAMETER_DIR");
    appendEnvironmentPaths(directories, "VIENNA_RNA_DATAPATH");
    appendEnvironmentPaths(directories, "VRNA_DATAPATH");
    if (const auto* prefix = std::getenv("CONDA_PREFIX"); prefix != nullptr && *prefix != '\0') {
        directories.emplace_back(fs::path(prefix) / "share" / "ViennaRNA");
    }
#ifdef _WIN32
    if (const auto* data = std::getenv("PROGRAMDATA"); data != nullptr && *data != '\0') {
        directories.emplace_back(fs::path(data) / "ViennaRNA");
    }
#else
    directories.emplace_back("/usr/local/share/ViennaRNA");
    directories.emplace_back("/usr/share/ViennaRNA");
#endif
    std::error_code error;
    auto current = fs::current_path(error);
    for (std::size_t depth{}; !error && !current.empty() && depth < 8U; ++depth) {
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
        if (auto file = regularFile(directory / fileName)) return *file;
    }
    throw std::invalid_argument("ViennaRNA parameter set '" + std::string(requested) +
                                "' was not found; set INTARNANEW_PARAMETER_DIR");
}

[[nodiscard]] auto section(const NumericSections& sections, const std::string_view name,
                           const std::size_t expected, const fs::path& path) -> const std::vector<double>& {
    const auto found = sections.find(name);
    if (found == sections.end() || found->second.size() != expected) {
        throw std::invalid_argument("ViennaRNA section # " + std::string(name) + " in " +
                                    path.string() + " must contain " + std::to_string(expected) + " values");
    }
    return found->second;
}

[[nodiscard]] auto scale(const double freeEnergy37, const double enthalpy,
                         const double temperatureKelvin,
                         const bool quantize = true) noexcept -> Energy {
    if (!std::isfinite(freeEnergy37) || !std::isfinite(enthalpy)) return infinity;
    if (std::abs(temperatureKelvin - referenceTemperatureKelvin) < 1e-9) return freeEnergy37 / 100.0;
    const auto value = enthalpy - (enthalpy - freeEnergy37) *
        temperatureKelvin / referenceTemperatureKelvin;
    return quantize ? std::trunc(value) / 100.0 : value / 100.0;
}

[[nodiscard]] auto scaledSection(const RawParameters& raw, const std::string_view energyName,
                                 const std::string_view enthalpyName, const std::size_t size,
                                 const fs::path& path, const double temperatureKelvin) -> std::vector<Energy> {
    const auto& energies = section(raw.numeric, energyName, size, path);
    const auto& enthalpies = section(raw.numeric, enthalpyName, size, path);
    std::vector<Energy> result;
    result.reserve(size);
    for (std::size_t index{}; index < size; ++index) {
        result.push_back(scale(energies[index], enthalpies[index], temperatureKelvin));
    }
    return result;
}

void parseSpecialLoops(FoldingParameters& destination, const RawParameters& raw,
                       const double temperatureKelvin, const fs::path& path) {
    for (const auto name : {std::string_view{"Triloops"}, std::string_view{"Tetraloops"},
                            std::string_view{"Hexaloops"}}) {
        const auto found = raw.textual.find(name);
        if (found == raw.textual.end()) continue;
        for (const auto& line : found->second) {
            std::istringstream input(line);
            std::string motif;
            std::string energyToken;
            std::string enthalpyToken;
            std::string extra;
            if (!(input >> motif >> energyToken >> enthalpyToken) || (input >> extra)) {
                throw std::invalid_argument("malformed special hairpin in " + path.string());
            }
            const auto energy = parseNumber(energyToken, path, 0U);
            const auto enthalpy = parseNumber(enthalpyToken, path, 0U);
            destination.specialHairpins.emplace(std::move(motif),
                                                scale(energy, enthalpy, temperatureKelvin));
        }
    }
}

[[nodiscard]] auto build(const fs::path& path, const double temperatureCelsius)
    -> std::shared_ptr<const FoldingParameters> {
    const auto temperatureKelvin = temperatureCelsius + 273.15;
    const auto raw = parseFile(path);
    auto result = std::make_shared<FoldingParameters>();
    result->stack = scaledSection(raw, "stack", "stack_enthalpies", 49U, path, temperatureKelvin);
    result->hairpin = scaledSection(raw, "hairpin", "hairpin_enthalpies", 31U, path, temperatureKelvin);
    result->bulge = scaledSection(raw, "bulge", "bulge_enthalpies", 31U, path, temperatureKelvin);
    result->internal = scaledSection(raw, "internal", "internal_enthalpies", 31U, path, temperatureKelvin);
    result->mismatchHairpin = scaledSection(raw, "mismatch_hairpin", "mismatch_hairpin_enthalpies",
                                            175U, path, temperatureKelvin);
    result->mismatchInternal = scaledSection(raw, "mismatch_internal", "mismatch_internal_enthalpies",
                                             175U, path, temperatureKelvin);
    result->mismatchInternal1n = scaledSection(raw, "mismatch_internal_1n", "mismatch_internal_1n_enthalpies",
                                               175U, path, temperatureKelvin);
    result->mismatchInternal23 = scaledSection(raw, "mismatch_internal_23", "mismatch_internal_23_enthalpies",
                                               175U, path, temperatureKelvin);
    result->mismatchMulti = scaledSection(raw, "mismatch_multi", "mismatch_multi_enthalpies",
                                          175U, path, temperatureKelvin);
    result->mismatchExterior = scaledSection(raw, "mismatch_exterior", "mismatch_exterior_enthalpies",
                                             175U, path, temperatureKelvin);
    result->dangle5 = scaledSection(raw, "dangle5", "dangle5_enthalpies", 35U, path, temperatureKelvin);
    result->dangle3 = scaledSection(raw, "dangle3", "dangle3_enthalpies", 35U, path, temperatureKelvin);
    result->int11 = scaledSection(raw, "int11", "int11_enthalpies", 1225U, path, temperatureKelvin);
    result->int21 = scaledSection(raw, "int21", "int21_enthalpies", 6125U, path, temperatureKelvin);
    result->int22 = scaledSection(raw, "int22", "int22_enthalpies", 9216U, path, temperatureKelvin);

    const auto& multiloop = section(raw.numeric, "ML_params", 6U, path);
    result->multiloopUnpaired = scale(multiloop[0], multiloop[1], temperatureKelvin);
    result->multiloopClosing = scale(multiloop[2], multiloop[3], temperatureKelvin);
    result->multiloopStem = scale(multiloop[4], multiloop[5], temperatureKelvin);
    const auto& ninio = section(raw.numeric, "NINIO", 3U, path);
    result->ninioSlope = scale(ninio[0], ninio[1], temperatureKelvin);
    result->ninioMaximum = ninio[2] / 100.0;
    const auto miscIterator = raw.numeric.find("Misc");
    if (miscIterator == raw.numeric.end() ||
        (miscIterator->second.size() != 4U && miscIterator->second.size() != 6U)) {
        throw std::invalid_argument("ViennaRNA section # Misc in " + path.string() +
                                    " must contain 4 or 6 values");
    }
    const auto& misc = miscIterator->second;
    result->terminalAu = scale(misc[2], misc[3], temperatureKelvin);
    result->logarithmicLoopSlope = misc.size() == 6U
        ? scale(misc[4], misc[5], temperatureKelvin, false)
        : scale(107.856, 0.0, temperatureKelvin, false);
    parseSpecialLoops(*result, raw, temperatureKelvin, path);
    return result;
}

struct CacheKey {
    fs::path path;
    double temperature{};
    friend auto operator<=>(const CacheKey&, const CacheKey&) = default;
};

} // namespace

auto loadFoldingParameters(const double temperatureCelsius, const std::string_view parameterSet)
    -> std::shared_ptr<const FoldingParameters> {
    if (!std::isfinite(temperatureCelsius) || temperatureCelsius <= -273.15) {
        throw std::invalid_argument("folding temperature must be finite and above absolute zero");
    }
    const CacheKey key{resolveParameterFile(parameterSet), temperatureCelsius};
    static std::mutex mutex;
    static std::map<CacheKey, std::weak_ptr<const FoldingParameters>> cache;
    {
        const std::scoped_lock lock(mutex);
        if (const auto found = cache.find(key); found != cache.end()) {
            if (auto existing = found->second.lock()) return existing;
            cache.erase(found);
        }
    }
    auto result = build(key.path, temperatureCelsius);
    {
        const std::scoped_lock lock(mutex);
        std::erase_if(cache, [](const auto& item) { return item.second.expired(); });
        if (cache.size() < 16U) cache.emplace(key, result);
    }
    return result;
}

} // namespace intarnanew::folding_detail
