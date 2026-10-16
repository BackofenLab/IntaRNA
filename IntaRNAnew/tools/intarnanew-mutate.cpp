#include "intarnanew/sequence.hpp"
#include "intarnanew/tools/mutations.hpp"

#include <charconv>
#include <cstdlib>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

namespace {

[[nodiscard]] auto help() -> std::string {
    return R"(intarnanew-mutate - reusable compensatory-mutation enumeration

Usage:
  intarnanew-mutate -q SEQUENCE -t SEQUENCE --pairs qIndex&tIndex[,..] [options]
  intarnanew-mutate -q SEQUENCE -t SEQUENCE --mutation-encoding G1C&U7G

Options:
  -q, --query SEQUENCE       raw query RNA
  -t, --target SEQUENCE      raw target RNA
      --q-index-origin N     first query coordinate (default: 1)
      --t-index-origin N     first target coordinate (default: 1)
      --pairs LIST           external query&target base-pair coordinates
  -g, --generator MODE      flip or any (default: flip)
  -f, --filter CLASS        GU, AU, or CG; repeatable
      --mutation-encoding M select one explicit mutation
  -h, --help                show this help

Output is deterministic TSV with mutation encoding and the four ww/wm/mw/mm
sequence combinations. Ranking measures require a prediction-library runner and
are intentionally outside this enumeration layer.
)";
}

[[nodiscard]] auto integer(const std::string_view text) -> std::expected<long long, std::string> {
    long long result{};
    const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), result);
    if (error != std::errc{} || end != text.data() + text.size()) return std::unexpected("invalid integer '" + std::string{text} + "'");
    return result;
}

[[nodiscard]] auto pairs(
    const std::string_view text,
    const intarnanew::Sequence& query,
    const intarnanew::Sequence& target) -> std::expected<std::vector<intarnanew::BasePair>, std::string> {
    std::vector<intarnanew::BasePair> result;
    std::size_t start{};
    while (start < text.size()) {
        const auto comma = text.find(',', start);
        const auto token = text.substr(start, comma == std::string_view::npos ? text.size() - start : comma - start);
        const auto amp = token.find('&');
        if (amp == std::string_view::npos) return std::unexpected("each pair must be queryIndex&targetIndex");
        auto qExternal = integer(token.substr(0U, amp));
        auto tExternal = integer(token.substr(amp + 1U));
        if (!qExternal) return std::unexpected(qExternal.error());
        if (!tExternal) return std::unexpected(tExternal.error());
        auto q = query.internalIndex(*qExternal);
        auto t = target.internalIndex(*tExternal);
        if (!q) return std::unexpected(q.error());
        if (!t) return std::unexpected(t.error());
        result.push_back({.target = *t, .query = *q});
        if (comma == std::string_view::npos) break;
        start = comma + 1U;
    }
    if (result.empty()) return std::unexpected("pair list is empty");
    return result;
}

} // namespace

int main(const int argc, const char* const* argv) {
    std::string queryText, targetText, pairText, encoding;
    long long queryOrigin = 1, targetOrigin = 1;
    auto generator = intarnanew::tools::MutationGenerator::flip;
    std::vector<intarnanew::tools::CandidateFilter> filters;
    for (int index = 1; index < argc; ++index) {
        const std::string_view option{argv[index]};
        const auto value = [&]() -> std::string_view {
            if (++index >= argc) { std::cerr << "intarnanew-mutate: missing value for " << option << '\n'; std::exit(2); }
            return argv[index];
        };
        if (option == "-h" || option == "--help") { std::cout << help(); return 0; }
        if (option == "-q" || option == "--query") queryText = value();
        else if (option == "-t" || option == "--target") targetText = value();
        else if (option == "--pairs") pairText = value();
        else if (option == "--mutation-encoding") encoding = value();
        else if (option == "--q-index-origin" || option == "--t-index-origin") {
            auto parsed = integer(value()); if (!parsed) { std::cerr << parsed.error() << '\n'; return 2; }
            if (option == "--q-index-origin") queryOrigin = *parsed; else targetOrigin = *parsed;
        } else if (option == "-g" || option == "--generator") {
            auto mode = value();
            if (mode == "flip") generator = intarnanew::tools::MutationGenerator::flip;
            else if (mode == "any") generator = intarnanew::tools::MutationGenerator::any;
            else { std::cerr << "intarnanew-mutate: generator must be flip or any\n"; return 2; }
        } else if (option == "-f" || option == "--filter") {
            auto filter = value();
            if (filter == "GU") filters.push_back(intarnanew::tools::CandidateFilter::removeGu);
            else if (filter == "AU") filters.push_back(intarnanew::tools::CandidateFilter::removeAu);
            else if (filter == "CG" || filter == "GC") filters.push_back(intarnanew::tools::CandidateFilter::removeCg);
            else { std::cerr << "intarnanew-mutate: this layer supports GU, AU, and CG filters\n"; return 2; }
        } else { std::cerr << "intarnanew-mutate: unknown option " << option << '\n'; return 2; }
    }
    if (queryText.empty() || targetText.empty()) { std::cerr << "intarnanew-mutate: query and target are required\n"; return 2; }
    try {
        intarnanew::Sequence query("query", queryText, queryOrigin), target("target", targetText, targetOrigin);
        std::vector<intarnanew::tools::MutationCandidate> mutations;
        if (!encoding.empty()) {
            auto mutation = intarnanew::tools::parseMutationEncoding(encoding, query, target);
            if (!mutation) { std::cerr << "intarnanew-mutate: " << mutation.error() << '\n'; return 1; }
            mutations.push_back(*mutation);
        } else {
            auto interactionPairs = pairs(pairText, query, target);
            if (!interactionPairs) { std::cerr << "intarnanew-mutate: " << interactionPairs.error() << '\n'; return 1; }
            auto generated = intarnanew::tools::enumerateMutations(query, target, *interactionPairs, generator, filters);
            if (!generated) { std::cerr << "intarnanew-mutate: " << generated.error() << '\n'; return 1; }
            mutations = std::move(*generated);
        }
        std::cout << "mutation\tqueryIndex\ttargetIndex\tbpWildtype\tbpMutated"
                     "\twwQuery\twwTarget\twmQuery\twmTarget"
                     "\tmwQuery\tmwTarget\tmmQuery\tmmTarget\n";
        for (const auto& mutation : mutations) {
            auto sequences = intarnanew::tools::applyMutation(query, target, mutation);
            if (!sequences) { std::cerr << "intarnanew-mutate: " << sequences.error() << '\n'; return 1; }
            std::cout << mutation.encoding(query, target) << '\t' << query.externalIndex(mutation.queryIndex)
                      << '\t' << target.externalIndex(mutation.targetIndex) << '\t'
                      << mutation.wildQuery << mutation.wildTarget << '\t'
                      << mutation.mutatedQuery << mutation.mutatedTarget << '\t'
                      << sequences->wildQuery << '\t' << sequences->wildTarget << '\t'
                      << sequences->wildQuery << '\t' << sequences->mutatedTarget << '\t'
                      << sequences->mutatedQuery << '\t' << sequences->wildTarget << '\t'
                      << sequences->mutatedQuery << '\t' << sequences->mutatedTarget << '\n';
        }
    } catch (const std::exception& error) {
        std::cerr << "intarnanew-mutate: " << error.what() << '\n'; return 1;
    }
    return 0;
}
