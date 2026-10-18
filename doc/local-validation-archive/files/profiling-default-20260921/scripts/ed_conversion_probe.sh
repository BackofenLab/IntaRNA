#!/usr/bin/env bash
set -euo pipefail
profile_root=$(cd "$(dirname "$0")/.." && pwd)
dependency_root=${INTARNA_DEPS:-/home/alex/CursorProjects/IntaRNA_optimization/intaRNA_legacy/.conda-env}
"$dependency_root/bin/x86_64-conda-linux-gnu-g++" -std=c++23 -O2 \
  -I"$profile_root/source/src" -I"$dependency_root/include" \
  "$profile_root/scripts/ed_conversion_probe.cpp" -o "$profile_root/diagnostics/ed_conversion_probe"
"$profile_root/diagnostics/ed_conversion_probe" > "$profile_root/diagnostics/ed-conversion.tsv"
cat "$profile_root/diagnostics/ed-conversion.tsv"
