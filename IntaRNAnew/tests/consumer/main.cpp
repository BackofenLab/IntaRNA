#include "intarnanew/accessibility.hpp"
#include "intarnanew/cli.hpp"
#include "intarnanew/compression.hpp"
#include "intarnanew/config.hpp"
#include "intarnanew/energy.hpp"
#include "intarnanew/folding.hpp"
#include "intarnanew/helix_blocks.hpp"
#include "intarnanew/output.hpp"
#include "intarnanew/output_plan.hpp"
#include "intarnanew/parallel.hpp"
#include "intarnanew/predictor.hpp"
#include "intarnanew/runner.hpp"
#include "intarnanew/sequence.hpp"
#include "intarnanew/tools/csv.hpp"
#include "intarnanew/tools/mutations.hpp"
#include "intarnanew/tools/pvalue.hpp"
#include "intarnanew/tools/statistics.hpp"
#include "intarnanew/tools/svg.hpp"
#include "intarnanew/types.hpp"

#include <array>
#include <cmath>
#include <span>

auto main() -> int {
    const intarnanew::Sequence sequence{"consumer", "ACGU"};
    if (sequence.size() != 4U || intarnanew::reverseComplement(sequence.str()) != "ACGU") {
        return 1;
    }

    const std::array<double, 3U> pValues{0.01, 0.04, 0.03};
    const auto adjusted = intarnanew::tools::adjustPValues(
        std::span<const double>{pValues},
        intarnanew::tools::AdjustmentMethod::benjaminiHochberg);
    if (!adjusted || adjusted->size() != pValues.size() ||
        std::abs((*adjusted)[0] - 0.03) > 1e-12) {
        return 2;
    }

    const auto compressed = intarnanew::gzipCompress("installed consumer");
    if (!compressed) return 3;
    const auto decompressed = intarnanew::gzipDecompress(*compressed);
    if (!decompressed || *decompressed != "installed consumer") return 4;

    intarnanew::Config config;
    config.energy = intarnanew::EnergyKind::basePair;
    config.mode = intarnanew::PredictionMode::exact;
    config.model = intarnanew::InteractionModel::singleSite;
    config.seed.required = false;
    config.target.accessibility = intarnanew::AccessibilityKind::disabled;
    config.query.accessibility = intarnanew::AccessibilityKind::disabled;
    config.output.number = 1U;
    const intarnanew::Sequence target{"target", "CCAACACC"};
    const intarnanew::Sequence query{"query", "GG"};
    const auto prediction = intarnanew::predictPair(config, target, query);
    if (!prediction || prediction->prediction.interactions.empty()) return 5;
    return 0;
}
