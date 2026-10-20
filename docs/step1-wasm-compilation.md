# Step 1: WebAssembly compilation (#286)

`scripts/wasm/build.sh` cross-compiles the unmodified CLI entry point and the
Autotools library source lists. `.github/workflows/build-wasm.yml` uses the official
[setup-emsdk action](https://github.com/emscripten-core/setup-emsdk) with Emscripten
4.0.23. It runs on pull requests, pushes to this repository's `master` branch,
and manual dispatches. Build tools (Autotools, pkg-config, a native C++ compiler,
curl, tar and bzip2) must be available.

## Reproduce locally

Activate Emscripten 4.0.23 with `source /path/to/emsdk/emsdk_env.sh`, then run:

```sh
bash scripts/wasm/build.sh
node scripts/wasm/smoke.mjs
```

`WASM_BUILD_DIR` changes the build directory; `JOBS` defaults to 2. Source
archives are downloaded into `build-wasm/downloads` and verified by pinned
SHA-256 hashes on every invocation. The build directory is ignored. Delete
`build-wasm/deps` and `build-wasm/sources` for a fresh dependency build after
changing dependency versions or compiler flags. Do not reuse this prefix with a
different SDK. Native libraries are never accepted from system pkg-config paths.

The script supplies absolute source and generated-header include paths for the
out-of-tree IntaRNA build. ViennaRNA 2.7.0 bundles dlib 19.24, whose Unicode
`uint32` string requires a generic `std::char_traits` removed in modern libc++.
The script backports [dlib's upstream `char32_t` typedef](https://github.com/davisking/dlib/blob/master/dlib/unicode/unicode.h)
in the downloaded source only, adds its matching serialization overload and
supplies a missing dependent template argument rejected by current Clang.
`scripts/wasm/patch-vienna.py` checks the expected source before each replacement.
It does not change the RNA algorithms.

## Dependencies and runtime

- zlib 1.3.1: static library for compressed accessibility and Boost.Iostreams.
- Boost 1.87.0: headers plus static filesystem, program_options, regex,
  iostreams, serialization and system libraries, built with Boost.Build's
  Emscripten toolset. Optional bzip2, lzma and zstd backends are disabled.
- ViennaRNA 2.7.0: static RNAlib, without OpenMP, pthreads, SIMD, external
  solvers, bindings, documentation or companion applications. Compiling only
  RNAlib avoids executing cross-compiled helper applications on the host.
- IntaRNA: C++23, bundled Kokkos mdspan, standard CLI, default energy models
  and accessibility calculations. Single-thread Emscripten builds omit
  `-pthread`; native builds retain their existing threading flags.

All C++ compilation and linking enables exceptions (`-fexceptions`), including
Boost and the CLI's input/error handling. The SDK's single-thread mutex support
retains ViennaRNA's normal setup/logger locking without requiring browser
`SharedArrayBuffer` or COOP/COEP headers (unavailable on ordinary GitHub Pages).

The output is an ES module factory (`createIntaRNA`) in
`build-wasm/site/wasm/intarna.js`, paired with `intarna.wasm`. The factory does
not run main automatically; it exports `FS` and `callMain`. The runtime exits
and flushes C++ streams after one call. Each prediction must use a new instance.
Memory grows on demand to a maximum of 1 GiB; stack size is 8 MiB. Allocation
failure is reported as a failed prediction; devices may impose lower limits.

See [Emscripten's Autotools guidance](https://emscripten.org/docs/compiling/Building-Projects.html)
for `emconfigure`/`emmake`. Generated configure files, dependency sources and WASM
binaries are build outputs, not tracked source. The artifact must keep JS and
WASM from the same build together. The Node smoke test runs a real
Turner/default-accessibility prediction and checks its energy, coordinates and
suboptimal count against the native CLI.
