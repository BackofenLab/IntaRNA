#!/bin/bash
set -euo pipefail
root=/home/alex/CodexProjects/IntaRNA
validation="$root/issue-267-validation"
minimum="$validation/vrna-2.7.0"
deps="$root/issue-246-validation/toolchain"
compiler="$root/issue-246-validation/toolchain-gcc14"
export PATH="$compiler/bin:$deps/bin:$PATH"
export LD_LIBRARY_PATH="$deps/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PKG_CONFIG_PATH="$minimum/lib/pkgconfig:$deps/lib/pkgconfig"
export OPENBLAS_NUM_THREADS=1
cd "$validation/minimum"
pkg-config --modversion RNAlib2 > "$validation/logs/minimum-dependency-version.log"
printf '#include <ViennaRNA/vrna_config.h>\nVRNA_VERSION\n' | \
  "$compiler/bin/x86_64-conda-linux-gnu-g++" -x c++ -E \
  -isystem "$minimum/include" -isystem "$deps/include" - \
  >> "$validation/logs/minimum-dependency-version.log"
bash autotools-init.sh > "$validation/logs/bootstrap-minimum.log" 2>&1
./configure \
  CC="$compiler/bin/x86_64-conda-linux-gnu-gcc" \
  CXX="$compiler/bin/x86_64-conda-linux-gnu-g++" \
  CPPFLAGS="-isystem $minimum/include -isystem $deps/include" \
  LDFLAGS="-L$minimum/lib -L$deps/lib -Wl,-rpath,$deps/lib" \
  --with-vrna="$minimum" --with-boost="$deps" --with-boost-libdir="$deps/lib" \
  --with-zlib="$deps" --prefix="$validation/install-minimum" \
  > "$validation/logs/configure-minimum.log" 2>&1
# Put the selected ViennaRNA archive ahead of the other dependency prefix.
# configure otherwise prepends the Boost library directory to AM_LDFLAGS.
link_flags="-L$minimum/lib -L$deps/lib -Wl,-rpath,$deps/lib"
make -j4 AM_LDFLAGS="$link_flags" > "$validation/logs/build-minimum.log" 2>&1
make tests -j4 AM_LDFLAGS="$link_flags" > "$validation/logs/tests-minimum.log" 2>&1
"$compiler/bin/x86_64-conda-linux-gnu-g++" -std=c++23 -O2 -fopenmp \
  -isystem "$minimum/include" -isystem "$deps/include" \
  doc/analysis/vrna-window-threads.cpp -o "$validation/vrna-window-threads-2.7.0" \
  -L"$minimum/lib" -L"$deps/lib" -Wl,-rpath,"$deps/lib" \
  $(pkg-config --libs RNAlib2) \
  > "$validation/logs/window-threads-minimum-build.log" 2>&1
"$validation/vrna-window-threads-2.7.0" setup-locked 4 3 \
  > "$validation/logs/window-threads-minimum.log" 2>&1
