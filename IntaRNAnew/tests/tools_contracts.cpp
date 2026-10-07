#include "intarnanew/sequence.hpp"
#include "intarnanew/tools/csv.hpp"
#include "intarnanew/tools/mutations.hpp"
#include "intarnanew/tools/pvalue.hpp"
#include "intarnanew/tools/statistics.hpp"
#include "intarnanew/tools/svg.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

namespace {

int failures{};

void require(const bool condition, const std::string_view description) {
    if (!condition) {
        std::cerr << "FAIL: " << description << '\n';
        ++failures;
    }
}

void close(const double observed, const double expected, const double tolerance,
           const std::string_view description) {
    require(std::abs(observed - expected) <= tolerance, description);
}

[[nodiscard]] auto dinucleotides(const std::string_view sequence)
    -> std::map<std::string, std::size_t> {
    std::map<std::string, std::size_t> result;
    for (std::size_t index = 0U; index + 1U < sequence.size(); ++index) {
        ++result[std::string{sequence.substr(index, 2U)}];
    }
    return result;
}

[[nodiscard]] auto sorted(std::string text) -> std::string {
    std::ranges::sort(text);
    return text;
}

void csvContracts() {
    std::istringstream input{
        "id;note;E\r\n"
        "a;\"semi;colon\";-2.5\r\n"
        "b;\"two\nlines and \"\"quote\"\"\";-1\r\n"};
    auto parsed = intarnanew::tools::readCsv(input);
    require(parsed.has_value(), "CSV parses RFC-style quoted fields and CRLF");
    if (!parsed) return;
    require(parsed->separator == ';', "CSV separator auto-detection");
    require(parsed->rows.size() == 2U, "CSV record count");
    require(parsed->rows[0][1] == "semi;colon", "CSV embedded separator");
    require(parsed->rows[1][1] == "two\nlines and \"quote\"", "CSV multiline/quote unescaping");
    auto serialized = intarnanew::tools::csvText(*parsed);
    require(serialized.has_value(), "CSV serialization succeeds");
    std::istringstream roundtrip{serialized.value_or("")};
    auto reparsed = intarnanew::tools::readCsv(roundtrip);
    require(reparsed.has_value() && reparsed->rows == parsed->rows,
            "CSV parse/write round trip preserves fields");

    intarnanew::tools::CsvTable second{
        .header = {"id", "probability"},
        .rows = {{"c", "0.1"}},
        .separator = ';'};
    std::array tables{*parsed, second};
    const std::array labels{std::string{"first"}, std::string{"second"}};
    auto fused = intarnanew::tools::fuseCsv(
        tables, labels, {.sourceColumn = "source", .deduplicate = false});
    require(fused.has_value(), "CSV schema-union fusion succeeds");
    if (fused) {
        require(fused->header == std::vector<std::string>{"id", "note", "E", "probability", "source"},
                "CSV schema union has stable first-occurrence order");
        require(fused->rows[2] == std::vector<std::string>{"c", "", "", "0.1", "second"},
                "CSV fusion fills absent fields and source label");
    }

    intarnanew::tools::CsvTable pathological{
        .header = {"value"},
        .rows = {{""}, {"comma,value"}, {"tab\tvalue"}, {"line\r\nvalue"},
                 {"quote\"value"}, {"UTF-8 π"}},
        .separator = ';'};
    auto pathologicalText = intarnanew::tools::csvText(pathological);
    std::istringstream pathologicalInput{pathologicalText.value_or("")};
    auto pathologicalRoundtrip = intarnanew::tools::readCsv(pathologicalInput);
    require(pathologicalRoundtrip.has_value() &&
                pathologicalRoundtrip->rows == pathological.rows,
            "CSV round trip preserves empty, control, quoted, and UTF-8 fields");
    std::istringstream duplicateHeader{"id;id\na;b\n"};
    require(!intarnanew::tools::readCsv(duplicateHeader),
            "CSV parser rejects duplicate column names");
    std::istringstream unterminated{"id;note\na;\"unfinished\n"};
    require(!intarnanew::tools::readCsv(unterminated),
            "CSV parser rejects unterminated quoted fields");
    intarnanew::tools::CsvTable malformed{
        .header = {"a", "a"}, .rows = {{"1", "2"}}, .separator = ';'};
    const std::array malformedTables{malformed};
    require(!intarnanew::tools::fuseCsv(malformedTables),
            "CSV fusion rejects malformed programmatic tables");
    std::ostringstream malformedOutput;
    require(!intarnanew::tools::writeCsv(malformed, malformedOutput),
            "CSV writer rejects malformed programmatic tables");
    require(!intarnanew::tools::parseFiniteNumber("1.2x") &&
                !intarnanew::tools::parseFiniteNumber("nan"),
            "numeric CSV conversion rejects trailing text and non-finite values");
}

void statisticContracts() {
    const std::array gaussian{-2.0, -1.0, 0.0, 1.0, 2.0};
    auto fit = intarnanew::tools::fitDistribution(
        gaussian, intarnanew::tools::DistributionKind::gaussian);
    require(fit.has_value(), "Gaussian fit succeeds");
    if (fit) {
        close(fit->location, 0.0, 1.0e-14, "Gaussian MLE location");
        close(fit->scale, std::sqrt(2.0), 1.0e-14, "Gaussian population MLE scale");
        auto probability = intarnanew::tools::interactionEnergyPValue(0.0, *fit);
        require(probability.has_value(), "Gaussian CDF evaluation succeeds");
        if (probability) close(*probability, 0.5, 1.0e-14, "Gaussian CDF at location");
    }
    const intarnanew::tools::DistributionFit gumbel{
        .kind = intarnanew::tools::DistributionKind::gumbel,
        .location = 3.0,
        .scale = 2.0};
    auto gumbelAtLocation = intarnanew::tools::cumulativeProbability(3.0, gumbel);
    require(gumbelAtLocation.has_value(), "Gumbel CDF evaluation succeeds");
    if (gumbelAtLocation) close(*gumbelAtLocation, std::exp(-1.0), 1.0e-14,
                                "Gumbel CDF at location");
    auto stableUpper = intarnanew::tools::tailProbability(
        203.0, gumbel, intarnanew::tools::ProbabilityTail::upper);
    require(stableUpper.has_value() && *stableUpper > 0.0,
            "Gumbel upper tail remains nonzero when 1-CDF would cancel");
    const intarnanew::tools::DistributionFit gev{
        .kind = intarnanew::tools::DistributionKind::gev,
        .location = -4.0,
        .scale = 1.5,
        .shape = 0.2};
    auto gevAtLocation = intarnanew::tools::cumulativeProbability(-4.0, gev);
    require(gevAtLocation.has_value(), "GEV CDF evaluation succeeds");
    if (gevAtLocation) close(*gevAtLocation, std::exp(-1.0), 1.0e-14,
                             "GEV CDF at location");

    std::vector<double> gumbelQuantiles;
    for (std::size_t index = 0U; index < 200U; ++index) {
        const double probability = (static_cast<double>(index) + 0.5) / 200.0;
        gumbelQuantiles.push_back(-3.0 - 2.0 * std::log(-std::log(probability)));
    }
    auto fittedGumbel = intarnanew::tools::fitDistribution(
        gumbelQuantiles, intarnanew::tools::DistributionKind::gumbel);
    require(fittedGumbel.has_value() && fittedGumbel->converged,
            "Gumbel MLE converges on deterministic population quantiles");
    if (fittedGumbel) {
        close(fittedGumbel->location, -3.0, 0.05, "Gumbel MLE recovers location");
        close(fittedGumbel->scale, 2.0, 0.05, "Gumbel MLE recovers scale");
    }
    std::vector<double> gevQuantiles;
    constexpr double gevLocation = -2.0;
    constexpr double gevScale = 1.2;
    constexpr double gevShape = 0.15;
    for (std::size_t index = 0U; index < 300U; ++index) {
        const double probability = (static_cast<double>(index) + 0.5) / 300.0;
        const double transformed = std::pow(-std::log(probability), -gevShape);
        gevQuantiles.push_back(gevLocation + gevScale * (transformed - 1.0) / gevShape);
    }
    auto fittedGev = intarnanew::tools::fitDistribution(
        gevQuantiles, intarnanew::tools::DistributionKind::gev);
    require(fittedGev.has_value() && fittedGev->converged,
            "GEV MLE converges on deterministic population quantiles");
    if (fittedGev) {
        close(fittedGev->location, gevLocation, 0.06, "GEV MLE recovers location");
        close(fittedGev->scale, gevScale, 0.06, "GEV MLE recovers scale");
        close(fittedGev->shape, gevShape, 0.04, "GEV MLE recovers shape");
    }

    const std::array background{-7.0, -6.0, -5.0, -4.0};
    auto empirical = intarnanew::tools::empiricalInteractionPValue(-5.5, background);
    require(empirical.has_value(), "empirical p-value succeeds");
    if (empirical) close(*empirical, 3.0 / 5.0, 1.0e-14, "plus-one empirical lower tail");

    const std::array probabilities{0.01, 0.04, 0.03, 0.002};
    auto bonferroni = intarnanew::tools::adjustPValues(
        probabilities, intarnanew::tools::AdjustmentMethod::bonferroni);
    require(bonferroni.has_value(), "Bonferroni adjustment succeeds");
    if (bonferroni) require(*bonferroni == std::vector<double>{0.04, 0.16, 0.12, 0.008},
                            "Bonferroni exact values");
    auto holm = intarnanew::tools::adjustPValues(
        probabilities, intarnanew::tools::AdjustmentMethod::holm);
    require(holm.has_value(), "Holm adjustment succeeds");
    if (holm) require(*holm == std::vector<double>{0.03, 0.06, 0.06, 0.008},
                      "Holm exact values");
    auto bh = intarnanew::tools::adjustPValues(
        probabilities, intarnanew::tools::AdjustmentMethod::benjaminiHochberg);
    require(bh.has_value(), "Benjamini-Hochberg adjustment succeeds");
    if (bh) require(*bh == std::vector<double>{0.02, 0.04, 0.04, 0.008},
                    "Benjamini-Hochberg exact values");
    auto hochberg = intarnanew::tools::adjustPValues(
        probabilities, intarnanew::tools::AdjustmentMethod::hochberg);
    require(hochberg.has_value() &&
                *hochberg == std::vector<double>{0.03, 0.04, 0.04, 0.008},
            "Hochberg exact values");
    auto by = intarnanew::tools::adjustPValues(
        probabilities, intarnanew::tools::AdjustmentMethod::benjaminiYekutieli);
    require(by.has_value(), "Benjamini-Yekutieli adjustment succeeds");
    if (by) {
        close((*by)[0], 1.0 / 24.0, 1.0e-14, "BY first adjusted value");
        close((*by)[1], 1.0 / 12.0, 1.0e-14, "BY second adjusted value");
        close((*by)[2], 1.0 / 12.0, 1.0e-14, "BY third adjusted value");
        close((*by)[3], 1.0 / 60.0, 1.0e-14, "BY fourth adjusted value");
    }
    const std::array invalid{0.1, -0.1};
    require(!intarnanew::tools::adjustPValues(
                 invalid, intarnanew::tools::AdjustmentMethod::holm),
            "invalid p-values are rejected");
}

void mutationContracts() {
    const intarnanew::Sequence query{"q", "GCAU", 1};
    const intarnanew::Sequence target{"t", "CGUA", 10};
    const std::array pairs{intarnanew::BasePair{.target = 0U, .query = 0U}};
    auto flipped = intarnanew::tools::enumerateMutations(
        query, target, pairs, intarnanew::tools::MutationGenerator::flip);
    require(flipped.has_value() && flipped->size() == 1U, "flip mutation enumeration");
    if (flipped && !flipped->empty()) {
        require(flipped->front().encoding(query, target) == "G1C&C10G",
                "mutation external-coordinate encoding");
        auto sequences = intarnanew::tools::applyMutation(query, target, flipped->front());
        require(sequences.has_value() && sequences->mutatedQuery == "CCAU" &&
                    sequences->mutatedTarget == "GGUA",
                "mutation application yields four-combination inputs");
    }
    auto parsed = intarnanew::tools::parseMutationEncoding("G1C&C10G", query, target);
    require(parsed.has_value(), "explicit mutation parsing succeeds");
    if (parsed && flipped && !flipped->empty()) {
        require(*parsed == flipped->front(), "explicit mutation parsing round trip");
    }
    auto any = intarnanew::tools::enumerateMutations(
        query, target, pairs, intarnanew::tools::MutationGenerator::any);
    require(any.has_value() && any->size() == 4U,
            "any generator emits every pairing combination mutating both bases");
    const std::array filter{intarnanew::tools::CandidateFilter::removeCg};
    auto filtered = intarnanew::tools::enumerateMutations(
        query, target, pairs, intarnanew::tools::MutationGenerator::any, filter);
    require(filtered.has_value() && filtered->empty(), "pair-class candidate filter");
    if (flipped && !flipped->empty()) {
        auto invalidMutation = flipped->front();
        invalidMutation.mutatedTarget = invalidMutation.wildTarget;
        require(!intarnanew::tools::applyMutation(query, target, invalidMutation),
                "manually constructed non-compensatory mutation is rejected");
    }

    const std::string input{"AACCGGUUAACCGG"};
    auto mono = intarnanew::tools::shuffleMononucleotides(input, 42U);
    require(mono.has_value() && sorted(*mono) == sorted(input),
            "mononucleotide shuffle preserves exact composition");
    require(mono == intarnanew::tools::shuffleMononucleotides(input, 42U),
            "mononucleotide shuffle is seed-deterministic");
    auto di = intarnanew::tools::shuffleDinucleotides(input, 42U);
    require(di.has_value() && dinucleotides(*di) == dinucleotides(input),
            "dinucleotide shuffle preserves directed dinucleotide counts");
    require(di.has_value() && di->front() == input.front() && di->back() == input.back(),
            "dinucleotide shuffle preserves trail endpoints");

    constexpr std::array<char, 4U> bases{'A', 'C', 'G', 'U'};
    bool exhaustiveInvariant = true;
    std::size_t sequenceCount = 1U;
    for (std::size_t length = 1U; length <= 6U && exhaustiveInvariant; ++length) {
        sequenceCount *= bases.size();
        for (std::size_t encoded = 0U; encoded < sequenceCount && exhaustiveInvariant; ++encoded) {
            auto remaining = encoded;
            std::string sequence(length, 'A');
            for (std::size_t position = 0U; position < length; ++position) {
                sequence[position] = bases[remaining % bases.size()];
                remaining /= bases.size();
            }
            for (std::uint64_t seed = 0U; seed < 3U; ++seed) {
                auto shuffled = intarnanew::tools::shuffleDinucleotides(sequence, seed);
                if (!shuffled || shuffled->size() != sequence.size() ||
                    shuffled->front() != sequence.front() || shuffled->back() != sequence.back() ||
                    dinucleotides(*shuffled) != dinucleotides(sequence)) {
                    exhaustiveInvariant = false;
                    break;
                }
            }
        }
    }
    require(exhaustiveInvariant,
            "dinucleotide shuffle invariants hold exhaustively through length six");
}

void pvalueOrchestrationContracts() {
    const intarnanew::tools::RandomScoreOptions single{
        .cardinality = 16U,
        .mode = intarnanew::tools::ShuffleMode::both,
        .preservation = intarnanew::tools::ShufflePreservation::mononucleotide,
        .randomSeed = 901U,
        .threads = 1U};
    auto evaluator = [](const std::string_view query, const std::string_view target)
        -> std::expected<double, std::string> {
        double score{};
        for (std::size_t index = 0; index < std::min(query.size(), target.size()); ++index) {
            score -= query[index] == target[index] ? 1.0 : 0.0;
        }
        return score;
    };
    auto oneThread = intarnanew::tools::sampleRandomInteractionScores(
        "AACCGGUU", "GGUUCCAA", single, evaluator);
    auto parallel = single;
    parallel.threads = 4U;
    auto fourThreads = intarnanew::tools::sampleRandomInteractionScores(
        "AACCGGUU", "GGUUCCAA", parallel, evaluator);
    require(oneThread.has_value() && fourThreads.has_value() && *oneThread == *fourThreads,
            "random-score samples are invariant to thread count");
}

void svgContracts() {
    const std::array profile{
        intarnanew::tools::ProfilePoint{1.0, -1.0},
        intarnanew::tools::ProfilePoint{2.0, std::nullopt},
        intarnanew::tools::ProfilePoint{3.0, -2.0}};
    intarnanew::tools::ProfileSvgOptions profileOptions;
    profileOptions.title = "A&B <profile>";
    auto first = intarnanew::tools::profileSvg(profile, profileOptions);
    auto second = intarnanew::tools::profileSvg(profile, profileOptions);
    require(first.has_value() && first == second, "profile SVG is deterministic");
    if (first) {
        require(first->find("A&amp;B &lt;profile&gt;") != std::string::npos,
                "profile SVG escapes XML text");
        const auto firstPolyline = first->find("<polyline");
        const auto secondPolyline = first->find("<polyline", firstPolyline + 1U);
        require(firstPolyline != std::string::npos && secondPolyline != std::string::npos,
                "missing profile point splits line segments");
        require(first->find("<circle") != std::string::npos,
                "single-point profile segments remain visible");
    }
    intarnanew::tools::HeatmapData heatmap{
        .xLabels = {"q1", "q2"},
        .yLabels = {"t1", "t2"},
        .values = {-3.0, 2.0, std::nullopt, -1.0}};
    auto svg = intarnanew::tools::heatmapSvg(heatmap);
    require(svg.has_value() && svg->find("q1 × t1: -3") != std::string::npos &&
                svg->find("q2 × t1: 0") != std::string::npos &&
                svg->find("q1 × t2: NA") != std::string::npos,
            "heatmap values, positive clamp, and missing cells are explicit");
    auto invalidOptions = profileOptions;
    invalidOptions.width = 100U;
    require(!intarnanew::tools::profileSvg(profile, invalidOptions),
            "SVG rejects unusably small dimensions");
    const std::array regions{
        intarnanew::tools::RegionSpan{"target&A", -5, 3},
        intarnanew::tools::RegionSpan{"targetB", 10, 14}};
    auto regionSvg = intarnanew::tools::regionsSvg(regions);
    require(regionSvg.has_value() &&
                regionSvg->find("target&amp;A: -5–3") != std::string::npos,
            "region SVG preserves coordinates and escapes identifiers");
    const std::array invalidRegions{intarnanew::tools::RegionSpan{"bad", 4, 2}};
    require(!intarnanew::tools::regionsSvg(invalidRegions),
            "region SVG rejects reversed spans");
}

} // namespace

int main() {
    csvContracts();
    statisticContracts();
    mutationContracts();
    pvalueOrchestrationContracts();
    svgContracts();
    if (failures == 0) {
        std::cout << "All native companion-tool contract tests passed.\n";
    }
    return failures == 0 ? 0 : 1;
}
