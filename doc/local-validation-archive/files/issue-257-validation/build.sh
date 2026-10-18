#!/usr/bin/env bash
set -euo pipefail
v=/home/alex/CodexProjects/IntaRNA/issue-257-validation
deps=/home/alex/CodexProjects/IntaRNA/issue-246-validation/toolchain
cc=/home/alex/CodexProjects/IntaRNA/issue-246-validation/toolchain-gcc14
export PATH="$cc/bin:$PATH" PKG_CONFIG_PATH="$deps/lib/pkgconfig" OPENBLAS_NUM_THREADS=1
cd "$1"
bash autotools-init.sh > "$v/logs/$2-bootstrap.log" 2>&1
./configure CC="$cc/bin/x86_64-conda-linux-gnu-gcc" CXX="$cc/bin/x86_64-conda-linux-gnu-g++" CPPFLAGS="-I$deps/include" LDFLAGS="-L$deps/lib -Wl,-rpath,$deps/lib" --with-vrna="$deps" --with-boost="$deps" --with-boost-libdir="$deps/lib" --prefix="$v/install-$2" "${@:3}" > "$v/logs/$2-configure.log" 2>&1
make -j2 > "$v/logs/$2-build.log" 2>&1
make tests -j2 > "$v/logs/$2-tests.log" 2>&1
