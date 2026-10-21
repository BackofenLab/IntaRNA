#!/usr/bin/env bash
# Build every linked dependency for wasm32; never reuse native libraries.
set -euo pipefail
repo=$(cd "$(dirname "$0")/../.." && pwd)
build=${WASM_BUILD_DIR:-"$repo/build-wasm"}
mkdir -p "$build"
build=$(cd "$build" && pwd)
prefix="$build/deps"
jobs=${JOBS:-2}
for tool in emcc em++ emar emranlib emconfigure emmake curl pkg-config node; do
  command -v "$tool" >/dev/null || { echo "Missing tool: $tool" >&2; exit 1; }
done
mkdir -p "$build/downloads" "$build/sources" "$prefix" "$build/site/wasm"
fetch() {
  local url=$1 name=$2 checksum=$3
  if [[ ! -f "$build/downloads/$name" ]]; then
    curl --fail --location --retry 3 "$url" -o "$build/downloads/$name.part"
    mv "$build/downloads/$name.part" "$build/downloads/$name"
  fi
  echo "$checksum  $build/downloads/$name" | sha256sum --check -
}
fetch https://zlib.net/fossils/zlib-1.3.1.tar.gz zlib-1.3.1.tar.gz \
  9a93b2b7dfdac77ceba5a558a580e74667dd6fede4585b91eefb60f03b72df23
fetch https://archives.boost.io/release/1.87.0/source/boost_1_87_0.tar.bz2 boost_1_87_0.tar.bz2 \
  af57be25cb4c4f4b413ed692fe378affb4352ea50fbe294a11ef548f4d527d89
fetch https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_7_x/ViennaRNA-2.7.0.tar.gz ViennaRNA-2.7.0.tar.gz \
  9a99fd68ed380894defb4d5e6a8a2871629270028cdf7f16f0a05da6e8c71473

# Keep pkg-config searches inside the WASM prefix, including configure probes.
export PKG_CONFIG_LIBDIR="$prefix/lib/pkgconfig"
export PKG_CONFIG_PATH="$PKG_CONFIG_LIBDIR"
export CFLAGS='-O3 -fexceptions'
export CXXFLAGS='-O3 -fexceptions'
export CPPFLAGS="-I$prefix/include"
export LDFLAGS="-L$prefix/lib -fexceptions"
if [[ ! -f "$prefix/lib/libz.a" ]]; then
  tar -xzf "$build/downloads/zlib-1.3.1.tar.gz" -C "$build/sources"
  (cd "$build/sources/zlib-1.3.1"
   emconfigure ./configure --static --prefix="$prefix"
   emmake make -j "$jobs"
   emmake make install)
fi
if [[ ! -f "$prefix/.boost-1.87.0-js-exceptions-complete" ]]; then
  tar -xjf "$build/downloads/boost_1_87_0.tar.bz2" -C "$build/sources"
  (cd "$build/sources/boost_1_87_0"
   ./bootstrap.sh --with-libraries=filesystem,program_options,regex,iostreams,serialization,system
   ./b2 -j "$jobs" toolset=emscripten variant=release link=static threading=single \
     exception-handling=on exception-handling-method=js \
     --with-filesystem --with-program_options --with-regex --with-iostreams \
     --with-serialization --with-system cxxflags="-std=c++23 -fexceptions" \
     -sNO_BZIP2=1 -sNO_LZMA=1 -sNO_ZSTD=1 -sZLIB_INCLUDE="$prefix/include" \
     -sZLIB_LIBPATH="$prefix/lib" --prefix="$prefix" install)
  touch "$prefix/.boost-1.87.0-js-exceptions-complete"
fi
if [[ ! -f "$prefix/lib/libRNA.a" ]]; then
  tar -xzf "$build/downloads/ViennaRNA-2.7.0.tar.gz" -C "$build/sources"
  (cd "$build/sources/ViennaRNA-2.7.0"
   CPPFLAGS="$CPPFLAGS -I$PWD/src" emconfigure ./configure --host=wasm32-unknown-emscripten --prefix="$prefix" \
     --disable-shared --enable-static --disable-openmp --disable-pthreads \
     --disable-lto --disable-simd --disable-vectorize --disable-mpfr \
     --without-svm --without-gsl --without-swig --without-perl --without-python \
     --without-doc --without-cla --without-check --without-kinfold \
     --without-forester --without-rnalocmin --without-rnaxplorer
   # Backport dlib's upstream char32_t typedef for modern libc++ char_traits.
   # Only the downloaded dependency is patched; RNAlib algorithms are unchanged.
   python3 "$repo/scripts/wasm/patch-vienna.py" "$PWD"
   emmake make -C src/ViennaRNA -j "$jobs"
   emmake make -C src/ViennaRNA install
   emmake make install-pkgconfigDATA)
fi

# Legacy Autoconf probes declare exit()/inflateEnd() with the wrong signatures,
# which wasm-ld rejects. Verify real dependency APIs before caching those probes.
em++ -std=c++23 -O2 -fexceptions -I"$prefix/include" \
  "$repo/scripts/wasm/dependency-probe.cpp" -L"$prefix/lib" \
  -lboost_filesystem -lboost_system -lboost_program_options -lboost_regex \
  -lboost_iostreams -lboost_serialization -lRNA -lz -o "$build/dependency-probe.js"
node "$build/dependency-probe.js"
export ac_cv_lib_boost_system_exit=yes
export ac_cv_lib_boost_filesystem_exit=yes
export ac_cv_lib_boost_program_options_exit=yes
export ac_cv_lib_boost_regex_exit=yes
export ac_cv_lib_boost_iostreams_exit=yes
export ac_cv_lib_z_inflateEnd=yes
# RNAlib's .pc embeds a native-only -Wl,-fno-lto; this is a static WASM build.
export VRNA_LIBS="-L$prefix/lib -lRNA -lm"

# Use the native Autotools source lists and CLI entry point, with a JS output.
(cd "$repo" && bash autotools-init.sh)
mkdir -p "$build/intarna"
(cd "$build/intarna"
 CXXFLAGS="$CXXFLAGS -Werror=inconsistent-missing-override" CPPFLAGS="-I$build/intarna/src -I$repo/src $CPPFLAGS" emconfigure "$repo/configure" --host=wasm32-unknown-emscripten \
   --disable-multithreading --disable-log-coloring --with-mdspan=kokkos \
   --with-vrna="$prefix" --with-boost="$prefix" \
   --with-boost-libdir="$prefix/lib" --with-zlib="$prefix"
 emmake make -C src -j "$jobs" EXEEXT=.js \
   LDFLAGS="$LDFLAGS -sMODULARIZE=1 -sEXPORT_ES6=1 -sEXPORT_NAME=createIntaRNA -sENVIRONMENT=web,worker,node -sINVOKE_RUN=0 -sEXIT_RUNTIME=1 -sEXPORTED_RUNTIME_METHODS=FS,callMain -sALLOW_MEMORY_GROWTH=1 -sMAXIMUM_MEMORY=1073741824 -sSTACK_SIZE=8388608"
 cp src/bin/IntaRNA.js "$build/site/wasm/intarna.js"
 cp src/bin/IntaRNA.wasm "$build/site/wasm/intarna.wasm")
# ES modules also work in Node smoke tests without changing package scope.
printf '{"type":"module"}\n' > "$build/site/wasm/package.json"
node "$repo/scripts/wasm/write-build-info.mjs" "$repo" "$build/site" \
  "$(emcc -dumpversion)" "$(git -C "$repo" rev-parse HEAD)" \
  "$(git -C "$repo" status --porcelain --untracked-files=no)"
echo "WASM assets: $build/site/wasm"
