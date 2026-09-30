# Vector/mdspan storage evaluation (issue #246)

The migration preserves scientific output and compact storage while retaining
GCC 14 support through bundled Kokkos headers when native mdspan is unavailable.
In the GCC 14 experiment, median wall-time changes range from 5.0% faster to
4.4% slower across workloads, with similar peak memory use. Performance is
mixed; these measurements do not establish a general speedup.

The earlier native-mdspan experiment with GCC 16 was effectively neutral
(1.7% faster to 1.2% slower). Both experiments, toolchains, and raw sample sets
are kept separately below. Treat this as a storage modernization.

## Scope and representation

All production `boost::numeric::ublas` storage aliases have been replaced by
owning containers in `IntaRNA/Matrix.h`. The dynamic programming recurrences,
traversal order, energy types and floating-point operations are unchanged.
`boost::multi_array` seed/helix recurrence tensors are outside this issue's
uBLAS scope and remain unchanged.

| Container | Logical shape | Stored elements | mdspan access |
| --- | --- | --- | --- |
| `Matrix<T>` | rows by columns | rows * columns | 2D row-major `[i,j]` |
| `UpperBandedMatrix<T>` | upper band including diagonal | rows * min(columns, upper+1) | 2D row-major `[i,j-i]` |
| `UpperTriangularMatrix<T>` | square upper triangle | n * (n+1) / 2 | 1D packed offset |

The triangular container uses a one-dimensional view because its packed row
lengths vary; a standard affine 2D layout would require quadratic rectangular
storage or a custom mapping. The band contains row-end padding, as did uBLAS,
and clamps superdiagonals that exceed the number of columns. In particular,
10,000 rows with ten superdiagonals take 110,000 elements, not 100 million.
The accessibility constructors still request their existing extra
superdiagonals for dangling-end probabilities.

Each vector owns its data. A short-lived mdspan is created at access time,
so there is no view pointer to repair after copying, moving or resizing.
Default resize preserves the overlapping logical cells; preserving resize
initializes new cells to zero/default values. Callers of non-preserving resize
must overwrite cells before reading. Const access outside a triangle/band
returns zero. Mutable access must address a stored cell, checked by assertions
in debug builds. Allocation-size overflow is rejected before allocating.
The intentionally small interface supports the operations used by IntaRNA,
not the full uBLAS expression/iterator interface. Triangular matrices must be
square and banded matrices must have lower bandwidth zero.

## Build requirements

The existing C++23/OpenMP requirement is unchanged: GCC 14 and Apple Clang
remain supported. Configure first compiles actual mutable 2D and const 1D
`std::mdspan` accesses over vector storage. If the native API is unavailable,
it compiles the same accesses using the bundled Kokkos implementation.
`INTARNA_USE_STD_MDSPAN` is set to `1` for native mdspan or `0` for Kokkos in
both the build configuration and the installed `IntaRNA/intarna_config.h`.
`Matrix.h` includes that public configuration and selects the corresponding
header and namespace; Kokkos types are not injected into `std`.

The default `--with-mdspan=auto` chooses native mdspan when supported.
`--with-mdspan=std` requires the native API and reports an actionable error
if it is unavailable. `--with-mdspan=kokkos` selects the bundled headers even
on a toolchain with native support, allowing reproducible backend comparisons.
Both implementations use the same owning containers and storage tests.

The vendored headers are copied unchanged from the latest
[Kokkos mdspan stable revision `8989f70749e28f337e6f7aa210db88659dba6f2f`](https://github.com/kokkos/mdspan/tree/8989f70749e28f337e6f7aa210db88659dba6f2f)
(as of 2026-09-29). Only upstream `include/mdspan` and `include/experimental`
headers are stored in `src/mdspan` and `src/experimental`. Both folders contain
an exact copy of the upstream `LICENSE` (Apache-2.0 WITH LLVM-exception).
No upstream build system, tests, benchmarks, or generated artifacts are
imported. Header-only use needs no additional upstream resources.

Source distributions contain both complete header trees and licenses.
Installation preserves their relative include layout under
`include/IntaRNA/mdspan` and `include/IntaRNA/experimental`, avoiding global
installation into another project's `mdspan` or `experimental` include tree.
Installed consumers use the configured backend and the existing IntaRNA
include path; no additional mdspan package or include flag is needed.

Linux CI retains its GCC 14 release/debug layout and checks automatic Kokkos
selection. Apple Clang CI checks native selection. Both compile every installed
public header independently and link/run an installed matrix consumer.
Boost remains a dependency for other project components and for the independent
reference containers in the storage regression tests.

## Reproduction

Build the parent revision `c5823e0` and this change in separate directories
with the same compiler, dependencies and configure flags. For each checkout:

```sh
bash autotools-init.sh
./configure CXX=/path/to/g++ CC=/path/to/gcc \
  --with-boost=/path/to/deps --with-boost-libdir=/path/to/deps/lib \
  --with-vrna=/path/to/deps \
  CPPFLAGS=-I/path/to/deps/include \
  LDFLAGS='-L/path/to/deps/lib -Wl,-rpath,/path/to/deps/lib'
make -j4
make -j4 check
```

Then, with no build or other benchmark running:

```sh
python3 tests/benchmark/compare-matrix-storage.py \
  /path/to/parent/src/bin/IntaRNA /path/to/candidate/src/bin/IntaRNA \
  /path/to/new-results-directory --repetitions=7 --cpu=2
```

The script generates seeded synthetic inputs, copies the repository's biological
fhlA/OxyS inputs, randomizes execution order, excludes one warmup per binary and
case, and retains stdout, stderr and each timing sample. It compares stdout
bytes for every execution, failing immediately on a mismatch or process error.
Each execution has a 180-second timeout that terminates both GNU time and its
IntaRNA child, retaining any captured stdout and stderr. Interrupting the
benchmark also terminates the process group. The biological input
files are included in source distributions so the script also runs from an
extracted release archive.
It records binary/input hashes, exact arguments, versions, affinity and controlled
environment. Wall time includes startup, loading and output; peak RSS comes from
GNU time. The test cases cover accessibility bands, dense predictors, seed
bulges, helix blocks, exact MFE, exact ensemble, seed extension and Nussinov
triangular storage. This is a bounded, single-host performance experiment.

## Backend-selection validation on 2026-09-29

- GCC 14.4 selects bundled Kokkos automatically; the release build passes
  33,264 assertions in 41 API cases and all 20 CLI golden cases.
- Native GCC 16.2 and bundled Kokkos/GCC 14.4 each pass the 28,923 storage
  assertions with AddressSanitizer and UndefinedBehaviorSanitizer. LeakSanitizer
  is disabled because of the sandbox's process-tracing limit.
- Configure checks cover GCC 14 automatic fallback, GCC 16 native selection,
  forced Kokkos on GCC 16, and actionable failures for forced native on GCC 14,
  an invalid backend option, and a missing bundled header tree.
- An installed GCC 14 consumer includes only installed IntaRNA/dependency
  headers, exercises all three matrix types, and links/runs successfully.
- All 66 installed public headers compile independently with GCC 14.
- All 30 imported headers and both license copies match pinned upstream bytes.
- An extracted source archive contains the complete bundled headers, licenses,
  and GCC 14 samples; configure selects Kokkos and its matrix consumer runs.

## Original native-backend validation on 2026-09-28

- Parent release: 4,341 assertions in 37 API cases; all 20 CLI golden cases pass.
- Candidate release and debug (`--enable-debug`): 33,264 assertions in 41 API
  cases; all 20 CLI golden cases pass in each build.
- Standalone storage tests: 28,923 assertions in four cases, passing with
  AddressSanitizer and UndefinedBehaviorSanitizer (`-O0 -g1
  -fsanitize=address,undefined -fno-omit-frame-pointer`). LeakSanitizer was
  disabled because the execution sandbox does not support its process tracing;
  no claim of leak-sanitizer coverage is made.
- The original native-only configure probe passed on GCC 16.2 and rejected
  a simulated missing `<mdspan>`. The current automatic mode instead selects
  Kokkos when native support is unavailable, as checked above.
- `make install` installs `Matrix.h`; a separate C++23 consumer including the
  installed accessibility, Nussinov and seed-extension headers compiles and runs.
- Benchmark process regression checks pass, including timeout cleanup of GNU
  time's child process and cleanup on interruption. Both checks reject runners
  without the corresponding cleanup.
- The benchmark runs from an extracted `make dist` archive: all ten workloads
  complete with matching outputs (one warmup and one measured run per binary,
  40 executions total). This is a packaging smoke check; the performance table
  below retains the original seven-repetition measurements.

The storage cases cover empty and rectangular shapes, logical preservation on
resize, zero-initialization of new cells, const structural zeros, the last
stored superdiagonal, packed triangular boundaries, non-scalar cells,
copy/move/swap independence, non-preserving resize, and allocation-size overflow.
Reference values come from uBLAS. New uBLAS primitive cells must be initialized
explicitly before comparison; its banded preserving resize is also unsuitable
as an oracle for some rectangular/empty transitions, so that reference is
rebuilt explicitly from overlapping cells.

## Original native-backend performance results (GCC 16)

Baseline: upstream `c5823e0`. Candidate: the storage migration in this change.

Both are optimized builds using conda-forge GCC/libstdc++ 16.2.0, Boost 1.85.0
and ViennaRNA 2.7.2 on Linux x86-64, AMD Ryzen 5 7530U. Configure adds
`-O3 -fno-strict-aliasing`; OpenMP is enabled with `--threads=1`.
Runs are pinned to CPU 2 with `LC_ALL=C`, `OPENBLAS_NUM_THREADS=1`,
`OMP_DYNAMIC=FALSE`. No build or sanitizer runs overlap these timings.

One warmup plus seven measured executions per binary/case: 160 executions
total, all stdout results byte-identical within and between variants.
The table shows median wall seconds and median peak RSS in MiB. A negative
time change means faster. RSS is a process-level measurement and includes
the executable, shared libraries and ViennaRNA allocations.

| Workload | uBLAS time (s) | mdspan time (s) | Change | uBLAS / mdspan RSS (MiB) |
| --- | ---: | ---: | ---: | ---: |
| biological-default | 0.0538 | 0.0542 | +0.8% | 16.80 / 16.67 |
| banded-default | 1.0801 | 1.0784 | -0.2% | 19.89 / 20.01 |
| banded-narrow | 0.7524 | 0.7592 | +0.9% | 17.71 / 17.59 |
| dense-no-accessibility | 3.9499 | 3.9981 | +1.2% | 15.05 / 15.05 |
| seed-bulges | 1.0204 | 1.0173 | -0.3% | 23.75 / 23.75 |
| helix-block | 0.1953 | 0.1975 | +1.1% | 18.38 / 18.26 |
| exact-mfe | 0.6855 | 0.6920 | +0.9% | 16.25 / 16.25 |
| exact-ensemble | 0.9845 | 0.9869 | +0.3% | 16.25 / 16.12 |
| seed-extension | 0.1306 | 0.1313 | +0.5% | 18.10 / 18.09 |
| triangular-base-pair | 2.4305 | 2.3897 | -1.7% | 14.25 / 14.38 |

Measured wall-time ranges (min–max across the seven measured samples):

| Workload | uBLAS (s) | mdspan (s) |
| --- | ---: | ---: |
| biological-default | 0.0534–0.0721 | 0.0529–0.0569 |
| banded-default | 1.0723–1.0850 | 1.0621–1.1033 |
| banded-narrow | 0.7463–0.7617 | 0.7476–0.7638 |
| dense-no-accessibility | 3.9071–4.1057 | 3.9527–4.0275 |
| seed-bulges | 1.0123–1.0264 | 1.0078–1.0262 |
| helix-block | 0.1949–0.1980 | 0.1934–0.2185 |
| exact-mfe | 0.6802–0.7005 | 0.6849–0.7092 |
| exact-ensemble | 0.9778–0.9980 | 0.9760–1.0010 |
| seed-extension | 0.1297–0.1322 | 0.1295–0.1320 |
| triangular-base-pair | 2.4121–2.5025 | 2.3744–2.4138 |

All samples (including explicitly marked warmups) and output hashes are in
[mdspan-storage-samples.tsv](mdspan-storage-samples.tsv). Full invocation
arguments and deterministic input generation are in the benchmark script.

Binary SHA-256 values:

- baseline: `e82b2031d59d2a5c07620cb7361ca1cb42eb5d0c052752762ee131c35d2ef0bb`
- candidate: `76799b89844581511e1010b72af10bcf2e28c47fe1111c21171733491e66e1c8`

These measurements describe the selected workloads on one host and one
toolchain. They do not establish a universal speedup. Short biological
runs include substantial startup overhead; differences near the observed
spread should be treated as inconclusive. Clang/macOS and other CPU
architectures still need their own measurements.

## Bundled Kokkos performance results (GCC 14, 2026-09-29)

The parent `c5823e0` and current candidate were rebuilt cleanly with the same
GCC 14.4.0 compiler and libstdc++ headers, Boost 1.85.0, ViennaRNA 2.7.2,
`-O3 -fno-strict-aliasing`, and OpenMP. Both link to the same dependency-prefix
libstdc++ 16.2 runtime; mdspan availability is determined by the GCC 14 headers.
The candidate automatically selects the pinned Kokkos backend. The host, CPU 2
affinity, workload arguments, random seed, environment controls, one warmup and
seven measured runs are as described above. No builds, header compilation or
sanitizer checks overlap these timings.

All 160 executions produce matching stdout for their respective workloads.
Median wall-time changes range from -5.0% to +4.4%.
These observations apply to this host and workload set; they do not establish
a universal performance improvement.

| Workload | uBLAS time (s) | Kokkos time (s) | Change | uBLAS / Kokkos RSS (MiB) |
| --- | ---: | ---: | ---: | ---: |
| biological-default | 0.0719 | 0.0728 | +1.2% | 16.42 / 16.30 |
| banded-default | 1.5545 | 1.5099 | -2.9% | 19.63 / 19.63 |
| banded-narrow | 1.1298 | 1.1367 | +0.6% | 17.33 / 17.21 |
| dense-no-accessibility | 5.8946 | 6.1527 | +4.4% | 14.68 / 14.68 |
| seed-bulges | 1.5445 | 1.4732 | -4.6% | 23.37 / 23.50 |
| helix-block | 0.2665 | 0.2666 | +0.0% | 18.00 / 18.00 |
| exact-mfe | 1.0318 | 1.0234 | -0.8% | 15.87 / 15.75 |
| exact-ensemble | 1.4179 | 1.4298 | +0.8% | 15.75 / 15.75 |
| seed-extension | 0.1764 | 0.1782 | +1.0% | 17.84 / 17.72 |
| triangular-base-pair | 3.5270 | 3.3494 | -5.0% | 14.00 / 13.88 |

Measured wall-time ranges across the seven samples:

| Workload | uBLAS (s) | Kokkos (s) |
| --- | ---: | ---: |
| biological-default | 0.0662–0.1053 | 0.0721–0.0936 |
| banded-default | 1.4561–2.2961 | 1.4573–1.5641 |
| banded-narrow | 1.0549–1.2632 | 1.0636–1.2372 |
| dense-no-accessibility | 5.6979–8.0424 | 5.9377–6.3798 |
| seed-bulges | 1.4834–2.1358 | 1.3990–1.5722 |
| helix-block | 0.2456–0.2855 | 0.2433–0.2901 |
| exact-mfe | 1.0026–1.1652 | 0.9963–1.0786 |
| exact-ensemble | 1.3865–1.5064 | 1.3904–1.4736 |
| seed-extension | 0.1670–0.1934 | 0.1621–0.1911 |
| triangular-base-pair | 3.3628–4.6398 | 3.2633–5.2531 |

All samples, warmups, and output hashes are in
[mdspan-storage-gcc14-samples.tsv](mdspan-storage-gcc14-samples.tsv).

Binary SHA-256 values:

- baseline: `f2d3293f7598980f4e539e925e1cac4bba35b08af688735d5d67250733329571`
- candidate: `b4842a322fa67d61990a3e4a623e19e3c5a049f2bf4d4db3213615c18b900b17`
