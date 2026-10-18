#!/bin/bash
set -euo pipefail
validation=/home/alex/CodexProjects/IntaRNA/issue-245-validation
toolchain=/home/alex/CodexProjects/IntaRNA/issue-246-validation/toolchain
export LD_LIBRARY_PATH="$toolchain/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PKG_CONFIG_PATH="$validation/install-release/lib/pkgconfig"
"$toolchain/bin/x86_64-conda-linux-gnu-g++" -std=c++23 -O3 \
  "$validation/benchmark.cpp" -L"$toolchain/lib" -Wl,-rpath,"$toolchain/lib" \
  $(pkg-config --cflags --libs IntaRNA) -L"$validation/install-release/lib" -leasylogging \
  -o "$validation/benchmark"
"$validation/benchmark" 100000 "$validation" > "$validation/benchmark.csv"
