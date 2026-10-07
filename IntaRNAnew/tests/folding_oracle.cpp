#include "intarnanew/accessibility.hpp"
#include "intarnanew/config.hpp"
#include "intarnanew/folding.hpp"
#include "intarnanew/sequence.hpp"
#include "intarnanew/types.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <string_view>

namespace {

auto failures = 0;

void check(const bool condition, const std::string_view message) {
    if (!condition) {
        std::cerr << "FAILED: " << message << '\n';
        ++failures;
    }
}

[[nodiscard]] auto close(const double observed, const double expected,
                         const double tolerance = 1e-10) noexcept -> bool {
    return std::abs(observed - expected) <= tolerance *
        std::max({1.0, std::abs(observed), std::abs(expected)});
}

[[nodiscard]] auto options() -> intarnanew::FoldingOptions {
    intarnanew::FoldingOptions result;
    result.parameterSet = "Turner04";
    result.maximumPairSpan = 100U;
    result.includeDangles = false;
    return result;
}

} // namespace

auto main() -> int {
    using namespace intarnanew;
    const auto rt = gasConstantKcal * 310.15;

    {
        const Sequence sequence("unpairable", "AAAAAA");
        const auto ensemble = makeTurnerFoldingEnsemble(sequence, options());
        check(close(ensemble->logPartition(), 0.0), "unpairable RNA has Z=1");
        check(close(ensemble->ensembleFreeEnergy(), 0.0), "unpairable RNA has ensemble energy zero");
        check(close(ensemble->jointUnpairedProbability({0U, 5U}), 1.0),
              "all positions of an unpairable RNA are jointly unpaired");
    }

    {
        const Sequence sequence("base-pair", "GAAAC");
        auto baseOptions = options();
        baseOptions.maximumPairSpan = sequence.size();
        const auto ensemble = makeBasePairFoldingEnsemble(sequence, baseOptions);
        const auto expectedLogZ = std::log1p(std::exp(1.0));
        check(close(ensemble->logPartition(), expectedLogZ),
              "base-pair model sums the empty and one-pair structures exactly");
        check(close(ensemble->ensembleFreeEnergy(), -expectedLogZ),
              "base-pair model uses its documented RT=1");
        check(close(ensemble->jointUnpairedProbability({0U, 4U}),
                    1.0 / (1.0 + std::exp(1.0))),
              "base-pair model computes exact joint interval probability");
        baseOptions.noLonelyPairs = true;
        const auto noLp = makeBasePairFoldingEnsemble(sequence, baseOptions);
        check(close(noLp->logPartition(), 0.0),
              "base-pair noLP removes the isolated pair");
    }

    {
        // Exactly two structures exist: all-unpaired and one GC hairpin with
        // three unpaired nucleotides. Turner04 assigns 5.4 kcal/mol.
        const Sequence sequence("single-hairpin", "GAAAC");
        const auto ensemble = makeTurnerFoldingEnsemble(sequence, options());
        const auto weight = std::exp(-5.4 / rt);
        const auto expectedLogZ = std::log1p(weight);
        check(close(ensemble->logPartition(), expectedLogZ),
              "single-hairpin partition equals its exhaustive two-structure sum");
        check(close(ensemble->jointUnpairedProbability({0U, 4U}), 1.0 / (1.0 + weight)),
              "whole-interval Pu is the all-unpaired structure probability");

        auto constrained = options();
        constrained.constraint = "p...p";
        const auto forcedPair = makeTurnerFoldingEnsemble(sequence, constrained);
        check(close(forcedPair->logPartition(), -5.4 / rt),
              "p constraints condition the ensemble on the required pair");
        check(forcedPair->jointUnpairedProbability({0U, 0U}) == 0.0,
              "a p-constrained position cannot be queried as unpaired");
    }

    {
        // Official ViennaRNA 2.7.2 public-executable differential fixture
        // (Turner04, 37 C, dangles=0): Eall=-0.17 kcal/mol and the global
        // all-unpaired structure probability is 0.755597. Legacy IntaRNA
        // 3.4.1 reports 0.7589416 because its accessibility adapter applies a
        // distinct exterior-loop convention; the provider model is explicit.
        const Sequence sequence("legacy-differential", "GCAAAUGC");
        const auto ensemble = makeTurnerFoldingEnsemble(sequence, options());
        const auto observedProbability = ensemble->jointUnpairedProbability({0U, 7U});
        check(close(ensemble->ensembleFreeEnergy(), -0.17, 0.012),
              "Turner ensemble free energy matches the black-box fixture");
        check(close(observedProbability, 0.755597, 2e-5),
              "Turner joint Pu matches the black-box fixture");

        auto dangleOptions = options();
        dangleOptions.includeDangles = true;
        const Sequence dangleSequence("dangle-two", "AGGGAAACCC");
        const auto dangles = makeTurnerFoldingEnsemble(dangleSequence, dangleOptions);
        check(close(dangles->ensembleFreeEnergy(), -1.80, 0.012),
              "Turner dangle-2 ensemble matches official RNAfold");
    }

    {
        // Public RNAeval 2.7.2 Turner04/dangles=0 forced-structure oracles.
        // These distinguish the three nucleotide axes of the asymmetric int21
        // table; symmetric A-only fixtures cannot catch a permutation.
        auto first = options();
        first.constraint = "pxxpxxxpxp";
        const auto int21TwoByOne = makeTurnerFoldingEnsemble(
            Sequence("int21-2x1", "AACGCAGCGU"), first);
        check(close(int21TwoByOne->ensembleFreeEnergy(), 8.90, 1e-10),
              "Turner int21 2x1 nucleotide order matches RNAeval");

        auto second = options();
        second.constraint = "pxpxxxpxxp";
        const auto int21OneByTwo = makeTurnerFoldingEnsemble(
            Sequence("int21-1x2", "AAUCGUAAGU"), second);
        check(close(int21OneByTwo->ensembleFreeEnergy(), 10.10, 1e-10),
              "Turner int21 1x2 nucleotide order matches RNAeval");
    }

    {
        // Turner99 and Andronescu07 use the valid six-value ViennaRNA Misc
        // section (DuplexInit, TerminalAU, LXC pairs).
        for (const auto& [name, hairpin] :
             std::initializer_list<std::pair<const char*, Energy>>{
                 {"Turner99", 5.70}, {"Andronescu07", 4.75}}) {
            auto named = options();
            named.parameterSet = name;
            const auto ensemble = makeTurnerFoldingEnsemble(
                Sequence("named-parameters", "GAAAC"), named);
            const auto weight = std::exp(-hairpin / rt);
            check(close(ensemble->logPartition(), std::log1p(weight), 1e-10),
                  "six-value Misc parameter set loads into the folding ensemble");
            check(close(ensemble->ensembleFreeEnergy(), -rt * std::log1p(weight), 1e-10),
                  "named folding parameter set preserves its own hairpin partition");
        }
    }

    {
        const Sequence sequence("lonely", "GAAAC");
        auto noLonelyPairs = options();
        noLonelyPairs.noLonelyPairs = true;
        const auto ensemble = makeTurnerFoldingEnsemble(sequence, noLonelyPairs);
        check(close(ensemble->logPartition(), 0.0), "noLP removes an isolated hairpin pair");

        const Sequence guSequence("gu-end", "GAAAU");
        auto noGuEnds = options();
        noGuEnds.noGuHelixEnds = true;
        const auto guEnsemble = makeTurnerFoldingEnsemble(guSequence, noGuEnds);
        check(close(guEnsemble->logPartition(), 0.0), "noGUend removes an isolated GU helix end");

        auto internalGu = options();
        internalGu.noGuHelixEnds = true;
        internalGu.constraint = "ppp....ppp";
        const Sequence internalGuSequence("internal-gu", "GGGAAAACUC");
        const auto accepted = makeTurnerFoldingEnsemble(internalGuSequence, internalGu);
        check(std::isfinite(accepted->logPartition()),
              "noGUend permits a GU pair stacked on both helix sides");

        auto terminalGu = options();
        terminalGu.noGuHelixEnds = true;
        terminalGu.constraint = "ppxxxxxxpp";
        bool rejectedTerminalGu{};
        try {
            static_cast<void>(makeTurnerFoldingEnsemble(
                Sequence("terminal-gu", "GGAAAAAACU"), terminalGu));
        } catch (const std::invalid_argument&) {
            rejectedTerminalGu = true;
        }
        check(rejectedTerminalGu, "noGUend rejects a GU pair at the outer helix end");
    }

    {
        // Public RNAfold/RNAeval 2.7.2 oracle for the single structure
        // ((...)(...)): dangles=0 gives 17.40 kcal/mol and dangles=2 gives
        // 13.20 kcal/mol. This pins dangle decoration between directly
        // adjacent multiloop stems.
        const Sequence sequence("forced-multiloop", "GGAAACGAAACC");
        auto noDangles = options();
        noDangles.constraint = "ppxxxppxxxpp";
        const auto dangleZero = makeTurnerFoldingEnsemble(sequence, noDangles);
        check(close(dangleZero->ensembleFreeEnergy(), 17.40, 0.012),
              "forced multiloop matches the public dangle-0 oracle");

        auto dangleTwo = noDangles;
        dangleTwo.includeDangles = true;
        const auto decorated = makeTurnerFoldingEnsemble(sequence, dangleTwo);
        check(close(decorated->ensembleFreeEnergy(), 13.20, 0.012),
              "forced multiloop matches the public dangle-2 oracle");
    }

    {
        // Eall/Zall are full-monomer quantities even when Pu/ED comes from a
        // local fold, a table, or is disabled. This is also a black-box legacy
        // contract: GAAAC under energy=B always has Z=1+e and Eall=-ln(Z).
        const Sequence sequence("summary", "GAAAC");
        Config config;
        config.energy = EnergyKind::basePair;
        config.query.accessibilityWindow = 4U;
        config.query.accessibilitySpan = 4U;
        const auto expectedLogZ = std::log1p(std::exp(1.0));

        config.query.accessibility = AccessibilityKind::compute;
        auto computedResult = makeAccessibility(sequence, config.query, config);
        check(computedResult.has_value(), "local computed accessibility constructs");
        if (computedResult) {
            check(computedResult.value()->ensembleLogPartition().has_value() &&
                  close(*computedResult.value()->ensembleLogPartition(), expectedLogZ),
                  "local computed accessibility exposes the global monomer partition");
            const auto physicalRt = gasConstantKcal * 310.15;
            const auto openingEnergy = std::trunc(expectedLogZ * 100.0) / 100.0;
            check(close(computedResult.value()->unpairedProbability({0U, 0U}),
                        std::exp(-openingEnergy / physicalRt)),
                  "base-pair Pu converts its RT=1 ensemble ED using physical RT");
            check(close(computedResult.value()->openingEnergy({0U, 0U}), openingEnergy),
                  "base-pair opening energy preserves public centikcal semantics");
        }

        config.query.accessibility = AccessibilityKind::disabled;
        auto disabledResult = makeAccessibility(sequence, config.query, config);
        check(disabledResult.has_value(), "disabled accessibility constructs with monomer summary");
        if (disabledResult) {
            check(disabledResult.value()->ensembleFreeEnergy().has_value() &&
                  close(*disabledResult.value()->ensembleFreeEnergy(), -expectedLogZ),
                  "disabled accessibility retains the global monomer free energy");
            check(close(disabledResult.value()->openingEnergy({0U, 4U}), 0.0),
                  "disabled interval opening energy remains zero");
        }

        const auto tablePath = std::filesystem::temp_directory_path() /
                               "intarnanew-folding-summary-table.txt";
        {
            std::ofstream table(tablePath);
            table << "1 1\n2 1\n3 1\n4 1\n5 1\n";
        }
        config.query.accessibility = AccessibilityKind::probabilitiesFile;
        config.query.accessibilityFile = tablePath.string();
        auto tableResult = makeAccessibility(sequence, config.query, config);
        check(tableResult.has_value(), "table accessibility constructs with monomer summary");
        if (tableResult) {
            check(tableResult.value()->ensembleLogPartition().has_value() &&
                  close(*tableResult.value()->ensembleLogPartition(), expectedLogZ),
                  "table accessibility retains the global monomer partition");
            check(close(tableResult.value()->positionUnpairedProbability(2U), 1.0),
                  "table values still govern interval accessibility");
        }
        std::error_code error;
        std::filesystem::remove(tablePath, error);
        check(!error, "temporary summary table fixture is removed");
    }

    {
        const Sequence sequence("constraints", "GAAAC", -2);
        auto constrained = options();
        constrained.constraint = "x:-2--2,b:2-2";
        const auto ensemble = makeTurnerFoldingEnsemble(sequence, constrained);
        check(close(ensemble->logPartition(), 0.0), "x/b endpoints are unpaired in the folding ensemble");
        check(close(ensemble->jointUnpairedProbability({0U, 4U}), 1.0),
              "x/b constraints preserve their conditional unpaired probability");
    }

    {
        const auto shapePath = std::filesystem::temp_directory_path() /
                               "intarnanew-folding-shape-oracle.txt";
        {
            std::ofstream shape(shapePath);
            shape << "1 G 0\n5 C 0\n";
        }
        const Sequence sequence("shape", "GAAAC");
        auto shaped = options();
        shaped.shapeFile = shapePath.string();
        shaped.shapeMethod = "Zb0.89";
        shaped.shapeConversion = "S";
        const auto ensemble = makeTurnerFoldingEnsemble(sequence, shaped);
        const auto weight = std::exp(-(5.4 + 2.0 * 0.89) / rt);
        check(close(ensemble->logPartition(), std::log1p(weight)),
              "Zarringhalam S conversion adds the documented paired-state penalties");
        std::error_code error;
        std::filesystem::remove(shapePath, error);
        check(!error, "temporary SHAPE fixture is removed");
    }

    {
        bool rejected{};
        try {
            validateShapeEncoding("Dm1.8m2", "");
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        check(rejected, "duplicate SHAPE method tags are rejected");
        rejected = false;
        try {
            validateShapeEncoding("D", "S");
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        check(rejected, "conversion without Zarringhalam method is rejected");
    }

    if (failures == 0) std::cout << "All Turner folding oracle tests passed.\n";
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
