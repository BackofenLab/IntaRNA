#!/usr/bin/env bash
set -euo pipefail
validation=/home/alex/CodexProjects/IntaRNA/issue-184-validation
dependencies=/home/alex/CodexProjects/IntaRNA/issue-246-validation/toolchain
compiler=/home/alex/CodexProjects/IntaRNA/issue-246-validation/toolchain-gcc14
source_dir=$1
mode=$2
export PATH="$compiler/bin:$PATH"
export PKG_CONFIG_PATH="$dependencies/lib/pkgconfig"
export OPENBLAS_NUM_THREADS=1
cd "$source_dir"
bash autotools-init.sh > "$validation/logs/$mode-bootstrap.log" 2>&1
./configure CXX="$compiler/bin/x86_64-conda-linux-gnu-g++" \
    CC="$compiler/bin/x86_64-conda-linux-gnu-gcc" CPPFLAGS="-I$dependencies/include" \
    LDFLAGS="-L$dependencies/lib -Wl,-rpath,$dependencies/lib" \
    --with-boost="$dependencies" --with-boost-libdir="$dependencies/lib" \
    --with-vrna="$dependencies" --prefix="$validation/install-$mode" \
    --enable-debug > "$validation/logs/$mode-configure.log" 2>&1
make -j2 > "$validation/logs/$mode-build.log" 2>&1
make tests -j2 > "$validation/logs/$mode-tests.log" 2>&1
