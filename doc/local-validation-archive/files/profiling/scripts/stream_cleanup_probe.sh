#!/usr/bin/env bash
set -euo pipefail
profile_root=$(cd "$(dirname "$0")/.." && pwd)
dependency_root=${INTARNA_DEPS:-/home/alex/CursorProjects/IntaRNA_optimization/intaRNA_legacy/.conda-env}
"$dependency_root/bin/x86_64-conda-linux-gnu-g++" -std=c++23 -O2 \
    "$profile_root/scripts/stream_cleanup_probe.cpp" \
    "$profile_root/source/src/IntaRNA/libIntaRNA.a" \
    -L"$dependency_root/lib" -Wl,-rpath,"$dependency_root/lib" \
    -lboost_regex -lboost_iostreams -lz -fopenmp \
    -o "$profile_root/results/stream_cleanup_probe"
"$profile_root/results/stream_cleanup_probe" "$profile_root/inputs/t100.fa" \
    > "$profile_root/results/stream-cleanup-probe.json"
cat "$profile_root/results/stream-cleanup-probe.json"
