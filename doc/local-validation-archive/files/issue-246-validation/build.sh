#!/bin/bash
set -euo pipefail
validation=/home/alex/CodexProjects/IntaRNA/issue-246-validation
toolchain=$validation/toolchain
source_dir=$1
label=$2
shift 2
export PATH="$toolchain/bin:$PATH"
export LD_LIBRARY_PATH="$toolchain/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PKG_CONFIG_PATH="$toolchain/lib/pkgconfig"
export OPENBLAS_NUM_THREADS=1
cd "$source_dir"
bash autotools-init.sh > "$validation/logs/$label-autotools.log" 2>&1
./configure CXX="$toolchain/bin/x86_64-conda-linux-gnu-g++" CC="$toolchain/bin/x86_64-conda-linux-gnu-gcc" CPPFLAGS="-I$toolchain/include" LDFLAGS="-L$toolchain/lib -Wl,-rpath,$toolchain/lib" --with-boost="$toolchain" --with-boost-libdir="$toolchain/lib" --with-vrna="$toolchain" --prefix="$validation/install-$label" "$@" > "$validation/logs/$label-configure.log" 2>&1
make -j4 > "$validation/logs/$label-build.log" 2>&1
make -j4 check > "$validation/logs/$label-check.log" 2>&1
