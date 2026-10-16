#pragma once

#include "intarnanew/sequence.hpp"
#include "intarnanew/types.hpp"

#include <memory>
#include <optional>
#include <string>
#include <string_view>

namespace intarnanew {

// Thermodynamic options for a single-RNA secondary-structure ensemble.
// A zero pair span means the whole sequence.  The implementation uses
// log-space arithmetic, therefore partitionScale is validated for command-line
// compatibility but is not needed to keep the recurrence finite.
struct FoldingOptions {
    double temperatureCelsius{37.0};
    std::string parameterSet{"Turner04"};
    Index maximumPairSpan{};
    Index maximumInternalLoop{30U};
    bool noLonelyPairs{};
    bool noGuHelixEnds{};
    bool includeDangles{true};
    double partitionScale{1.07};
    std::string constraint;
    std::string shapeFile;
    std::string shapeMethod;
    std::string shapeConversion;
};

// Immutable exact global ensemble for the documented Turner nearest-neighbour
// model.  Interval probabilities are joint probabilities, not products of
// single-position marginals.
class FoldingEnsemble {
public:
    virtual ~FoldingEnsemble() = default;

    [[nodiscard]] virtual auto logPartition() const noexcept -> double = 0;
    [[nodiscard]] virtual auto ensembleFreeEnergy() const noexcept -> Energy = 0;
    [[nodiscard]] virtual auto jointUnpairedProbability(Interval interval) const -> double = 0;
};

[[nodiscard]] auto makeTurnerFoldingEnsemble(
    const Sequence& sequence,
    const FoldingOptions& options) -> std::unique_ptr<FoldingEnsemble>;

// Nussinov-style single-RNA ensemble used with the public base-pair energy
// model. Each intramolecular pair contributes -1 and RT is exactly 1.
[[nodiscard]] auto makeBasePairFoldingEnsemble(
    const Sequence& sequence,
    const FoldingOptions& options) -> std::unique_ptr<FoldingEnsemble>;

// Public validation helpers are shared by the folding and accessibility APIs.
// They throw std::invalid_argument with a stable diagnostic on malformed data.
void validateShapeEncoding(std::string_view method, std::string_view conversion);

} // namespace intarnanew
