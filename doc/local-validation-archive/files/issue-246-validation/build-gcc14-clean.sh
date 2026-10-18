#!/bin/bash
set -euo pipefail
validation=/home/alex/CodexProjects/IntaRNA/issue-246-validation
compiler=$validation/toolchain-gcc14
deps=$validation/toolchain
source_dir=$1
label=$2
shift 2
export PATH="$compiler/bin:$deps/bin:$PATH"
export LD_LIBRARY_PATH="$deps/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PKG_CONFIG_PATH="$deps/lib/pkgconfig"
export OPENBLAS_NUM_THREADS=1
cd "$source_dir"
bash autotools-init.sh > "$validation/logs/$label-autotools.log" 2>&1
./configure CXX="$compiler/bin/x86_64-conda-linux-gnu-g++" CC="$compiler/bin/x86_64-conda-linux-gnu-gcc" CPPFLAGS="-I$deps/include" LDFLAGS="-L$deps/lib -Wl,-rpath,$deps/lib" --with-boost="$deps" --with-boost-libdir="$deps/lib" --with-vrna="$deps" --prefix="$validation/install-$label" "$@" > "$validation/logs/$label-configure.log" 2>&1
make clean > "$validation/logs/$label-clean.log" 2>&1
make -j4 > "$validation/logs/$label-build.log" 2>&1
make -j4 check > "$validation/logs/$label-check.log" 2>&1
make install > "$validation/logs/$label-install.log" 2>&1
