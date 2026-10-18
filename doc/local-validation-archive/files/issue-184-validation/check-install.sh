#!/usr/bin/env bash
set -euo pipefail
validation=/home/alex/CodexProjects/IntaRNA/issue-184-validation
dependencies=/home/alex/CodexProjects/IntaRNA/issue-246-validation/toolchain
compiler=/home/alex/CodexProjects/IntaRNA/issue-246-validation/toolchain-gcc14/bin/x86_64-conda-linux-gnu-g++
prefix="$validation/install-release"
if test "${1:-all}" != consumer; then
for header in "$prefix/include/IntaRNA/"*.h; do
    printf '#include <IntaRNA/%s>\nint main() { return 0; }\n' "$(basename "$header")" > "$validation/header.cpp"
    "$compiler" -std=c++23 -I"$prefix/include" -I"$dependencies/include" \
        -fsyntax-only "$validation/header.cpp"
done
echo 'All installed public headers compiled independently'
fi
cat > "$validation/consumer.cpp" <<'CPP'
#include <IntaRNA/PredictorEvalOnly.h>
#include <IntaRNA/AccessibilityDisabled.h>
#include <IntaRNA/InteractionEnergyBasePair.h>
#include <IntaRNA/OutputHandlerInteractionList.h>
INITIALIZE_EASYLOGGINGPP
int main() {
    using namespace IntaRNA;
    RnaSequence t("t", "GG"), q("q", "CC");
    AccessibilityDisabled a(t, 0, nullptr), b(q, 0, nullptr);
    ReverseAccessibility reversed(b);
    InteractionEnergyBasePair energy(a, reversed);
    OutputHandlerInteractionList output(OutputConstraint(), 10);
    PredictorEvalOnly eval(energy, output, nullptr,
        PredictorEvalOnly::parseInteractions("1||&1||", t, q));
    eval.predict();
    return output.reported() == 1 && (*output.begin())->energy == -200 ? 0 : 1;
}
CPP
export PKG_CONFIG_PATH="$prefix/lib/pkgconfig:$dependencies/lib/pkgconfig"
"$compiler" -std=c++23 -DELPP_FEATURE_PERFORMANCE_TRACKING -DELPP_NO_DEFAULT_LOG_FILE "$validation/consumer.cpp" -I"$dependencies/include" \
    -L"$dependencies/lib" -Wl,-rpath,"$dependencies/lib" \
    $(pkg-config --cflags --libs IntaRNA) -leasylogging -o "$validation/consumer"
"$validation/consumer"
echo 'Installed headers and pkg-config evaluation consumer passed'
