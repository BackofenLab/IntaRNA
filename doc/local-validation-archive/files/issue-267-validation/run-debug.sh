#!/usr/bin/env bash
set -euo pipefail
root=/home/alex/CodexProjects/IntaRNA
validation="$root/issue-267-validation"
deps="$root/issue-246-validation/toolchain"
gcc14="$root/issue-246-validation/toolchain-gcc14/bin"
cd "$validation/debug"
bash autotools-init.sh > "$validation/logs/bootstrap-debug.log" 2>&1
./configure \
  CC="$gcc14/x86_64-conda-linux-gnu-gcc" \
  CXX="$gcc14/x86_64-conda-linux-gnu-g++" \
  CPPFLAGS="-I$deps/include" \
  LDFLAGS="-L$deps/lib -Wl,-rpath,$deps/lib" \
  --with-vrna="$deps" \
  --with-boost="$deps" \
  --with-boost-libdir="$deps/lib" \
  --with-zlib="$deps" \
  --enable-debug \
  --prefix="$validation/install-debug" > "$validation/logs/configure-debug.log" 2>&1
make -j2 > "$validation/logs/build-debug.log" 2>&1
make tests -j2 > "$validation/logs/tests-debug.log" 2>&1
printf 'Debug configure, build, and full tests passed\n'
tail -18 "$validation/logs/tests-debug.log"
