#!/usr/bin/env bash
set -euo pipefail
validation=/home/alex/CodexProjects/IntaRNA/issue-212-validation
dependencies=/home/alex/CodexProjects/IntaRNA/issue-246-validation/toolchain
compiler=/home/alex/CodexProjects/IntaRNA/issue-246-validation/toolchain-gcc14
export PATH="$compiler/bin:$PATH"
export PKG_CONFIG_PATH="$dependencies/lib/pkgconfig"
export OPENBLAS_NUM_THREADS=1
cd "$validation/review-debug"
bash autotools-init.sh > "$validation/logs/review-debug-bootstrap.log" 2>&1
./configure CXX="$compiler/bin/x86_64-conda-linux-gnu-g++" \
    CC="$compiler/bin/x86_64-conda-linux-gnu-gcc" CPPFLAGS="-I$dependencies/include" \
    CXXFLAGS="-O0 -g -fsanitize=undefined,float-cast-overflow -fno-sanitize-recover=all" \
    LDFLAGS="-L$dependencies/lib -Wl,-rpath,$dependencies/lib -fsanitize=undefined,float-cast-overflow" \
    --with-boost="$dependencies" --with-boost-libdir="$dependencies/lib" \
    --with-vrna="$dependencies" --prefix="$validation/install-review-debug" \
    --enable-debug > "$validation/logs/review-debug-configure.log" 2>&1
make tests -j2 > "$validation/logs/review-debug-tests.log" 2>&1
