#!/bin/bash
set -euo pipefail
validation=/home/alex/CodexProjects/IntaRNA/issue-267-validation
compiler=/home/alex/CodexProjects/IntaRNA/issue-246-validation/toolchain-gcc14
deps=/home/alex/CodexProjects/IntaRNA/issue-246-validation/toolchain
export PATH="$compiler/bin:$deps/bin:$PATH"
export LD_LIBRARY_PATH="$deps/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PKG_CONFIG_PATH="$deps/lib/pkgconfig"
export OPENBLAS_NUM_THREADS=1
cd "$validation/upstream/ViennaRNA-2.7.0"
./configure CC="$compiler/bin/x86_64-conda-linux-gnu-gcc" \
  CXX="$compiler/bin/x86_64-conda-linux-gnu-g++" \
  CPPFLAGS="-I$deps/include" LDFLAGS="-L$deps/lib -Wl,-rpath,$deps/lib" \
  --prefix="$validation/vrna-2.7.0" --disable-shared --enable-static \
  --disable-lto --without-svm --without-gsl --without-swig \
  --without-perl --without-python --without-doc --without-cla \
  --without-check --disable-unittests --without-kinfold --without-forester \
  --without-rnalocmin --without-rnaxplorer \
  > "$validation/logs/vrna-2.7.0-configure.log" 2>&1
# ViennaRNA 2.7.0 omits utils CPPFLAGS when SVM support is disabled.
# Supply its ordinary project flags without changing the upstream sources.
make -C src/ViennaRNA -j2 'libRNA_utils_la_CPPFLAGS=$(AM_CPPFLAGS)' > "$validation/logs/vrna-2.7.0-build.log" 2>&1
make -C src/ViennaRNA install 'libRNA_utils_la_CPPFLAGS=$(AM_CPPFLAGS)' > "$validation/logs/vrna-2.7.0-install.log" 2>&1
make install-pkgconfigDATA >> "$validation/logs/vrna-2.7.0-install.log" 2>&1
