#!/usr/bin/env bash
set -euo pipefail
task_root=/home/alex/CodexProjects/IntaRNA
task_validation="$task_root/issue-267-validation"
task_deps="$task_root/issue-246-validation/toolchain"
task_cxx="$task_root/issue-246-validation/toolchain-gcc14/bin/x86_64-conda-linux-gnu-g++"
cd "$task_root/pr-worktrees/issue-267"
make install > "$task_validation/logs/install-release.log" 2>&1
cat > "$task_validation/installed-consumer.cpp" <<'CPP'
#include <IntaRNA/VrnaHandler.h>
#include <IntaRNA/RnaSequence.h>
#include <IntaRNA/AccessibilityVrna.h>
INITIALIZE_EASYLOGGINGPP
int main() {
    IntaRNA::VrnaHandler handler;
    IntaRNA::RnaSequence sequence("test", "GGGGAAAACCCC");
    IntaRNA::AccessibilityVrna accessibility(sequence, 6, nullptr, handler, 10);
    return accessibility.getED(0, 2) < IntaRNA::E_INF ? 0 : 1;
}
CPP
export PKG_CONFIG_PATH="$task_validation/install/lib/pkgconfig:$task_deps/lib/pkgconfig"
"$task_cxx" -std=c++23 -fopenmp -DELPP_FEATURE_PERFORMANCE_TRACKING -DELPP_NO_DEFAULT_LOG_FILE \
  "$task_validation/installed-consumer.cpp" \
  -L"$task_deps/lib" -Wl,-rpath,"$task_deps/lib" \
  $(pkg-config --cflags --libs IntaRNA) -leasylogging -o "$task_validation/installed-consumer"
"$task_validation/installed-consumer"
for task_header in "$task_validation/install/include/IntaRNA/"*.h; do
  printf '#include <IntaRNA/%s>\nint main() { return 0; }\n' "$(basename "$task_header")" > "$task_validation/header-probe.cpp"
  "$task_cxx" -std=c++23 -I"$task_validation/install/include" -I"$task_deps/include" \
    -fsyntax-only "$task_validation/header-probe.cpp"
done
printf 'Installed consumer and %s standalone public headers passed.\n' "$(find "$task_validation/install/include/IntaRNA" -maxdepth 1 -name '*.h' | wc -l)"
