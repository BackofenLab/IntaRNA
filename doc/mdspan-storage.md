# Vector/mdspan storage evaluation (issue #246)

The migration preserves scientific output and compact storage. In the measured
single-host experiment, performance is effectively neutral: median wall-time
changes range from 1.7% faster to 1.2% slower, with similar peak memory use.
There is no substantial performance gain to justify the stricter toolchain
requirement by itself. Treat this as a storage modernization; retain the
measurements below when deciding whether to adopt the new build requirement.

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

Configure compiles actual mutable 2D and const 1D `std::mdspan` accesses over
vector storage. An unsupported standard library fails with an actionable error,
including when its compiler accepts `-std=c++23`. Public-header consumers also
need `<mdspan>` and C++23.

Native support is available in [GCC 16/libstdc++](https://gcc.gnu.org/onlinedocs/libstdc++/manual/status.html)
and [LLVM 18/libc++](https://releases.llvm.org/18.1.8/projects/libcxx/docs/Status/Cxx23.html).
The local validation below uses GCC. GitHub Actions also passed the Apple Clang
17 release build on macOS 15, including the full test suite, independent public
header compilation, and the installed pkg-config consumer ([CI run](https://github.com/BackofenLab/IntaRNA/actions/runs/36453446578/job/109033658550)).
The Linux CI jobs select GCC 16 from conda-forge: GCC 14's standard library
does not provide mdspan and correctly fails the configure check. Compiler
implementation packages avoid activation hooks overriding configure's release
and debug optimization flags.
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

## Validation on 2026-09-28

- Parent release: 4,341 assertions in 37 API cases; all 20 CLI golden cases pass.
- Candidate release and debug (`--enable-debug`): 33,264 assertions in 41 API
  cases; all 20 CLI golden cases pass in each build.
- Standalone storage tests: 28,923 assertions in four cases, passing with
  AddressSanitizer and UndefinedBehaviorSanitizer (`-O0 -g1
  -fsanitize=address,undefined -fno-omit-frame-pointer`). LeakSanitizer was
  disabled because the execution sandbox does not support its process tracing;
  no claim of leak-sanitizer coverage is made.
- Native GCC 16.2 configure probe passes. A GCC 16.2 configure run with an
  intentionally unavailable `<mdspan>` fails at the new probe with the expected
  diagnostic, after the other C++23 checks pass.
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

## Same-host performance results

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
