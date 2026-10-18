#!/bin/bash
set -euo pipefail
validation=/home/alex/CodexProjects/IntaRNA/issue-245-validation
toolchain=/home/alex/CodexProjects/IntaRNA/issue-246-validation/toolchain
for header in Accessibility AccessibilityFromStream; do
  printf '#include <IntaRNA/%s.h>\nint main() { return 0; }\n' "$header" > "$validation/header.cpp"
  "$toolchain/bin/x86_64-conda-linux-gnu-g++" -std=c++23 -Wno-deprecated-declarations \
    -I"$validation/install-release/include" -I"$toolchain/include" -fsyntax-only "$validation/header.cpp"
done
