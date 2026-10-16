# IntaRNAnew

IntaRNAnew is an independent, standalone C++23 implementation of RNA–RNA
interaction prediction. Its scientific core, command-line interface, FASTA
handling, accessibility calculation, thermodynamic scoring, ensemble arithmetic,
tracking output, and tests are all native C++. It has no Python, Boost, ViennaRNA,
zlib, or other runtime/library dependency beyond the C++ standard library. The
nearest-neighbor model reads public ViennaRNA v2 `.par` data files but does not
link to the ViennaRNA library.

The implementation was designed from the public IntaRNA model and verified with
end-to-end RNA–RNA interaction behavior from IntaRNA Legacy. No source code was
copied from either historical implementation.

## Build

GCC 13 or Clang 17 and GNU Make are sufficient:

```sh
make -j
make check
```

The project also provides a CMake build for systems with CMake 3.25 or newer.
Both build systems require a C++23 standard library and a working thread
implementation. The installed libraries are static, so the public ABI does not
depend on platform-specific shared-library export annotations.

The Make build tracks included headers with compiler-generated dependency files.
Its debug and sanitizer targets run clean and build sequentially, avoiding a
parallel clean/build race:

```sh
make debug
make sanitize
```

To stage an installation, smoke every installed executable, and compile a
standalone C++ consumer against only the installed headers and archives:

```sh
make install-check
```

For a CMake build and install:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
cmake --install build --prefix /desired/prefix
```

An installed CMake consumer uses `find_package(IntaRNAnew CONFIG REQUIRED)` and
links `IntaRNAnew::core` or `IntaRNAnew::tools`. A complete consumer smoke project
is provided in `tests/consumer`.

## Quick start

```sh
./IntaRNAnew --target=CCAACCCACC --query=GGGG \
  --energy=B --acc=N --seedBP=3 --mode=M --model=S \
  --outMode=C --outNumber=10
```

Inputs may be literal IUPAC RNA sequences, FASTA files, multi-FASTA files, or
`STDIN`. Run `./IntaRNAnew --fullhelp` for the supported compatibility surface.

## Architecture

- `Sequence` and `SequenceReader`: validated IUPAC RNA and FASTA parsing.
- `AccessibilityProvider`: disabled, table-backed, and native log-space
  pseudoknot-free accessibility models.
- `HybridEnergyModel`: Nussinov-like and independent nearest-neighbor models.
- `Predictor`: endpoint-aware antiparallel DAG dynamic programming with seed
  automata, exact/heuristic/seed-only execution, helix filtering, k-best overlap
  policies, and log-space ensemble aggregation.
- `predictPair`: reusable in-process evaluation with the same region/window
  planning and global reduction contract as one command-line sequence pair.
- `OutputFormatter`: normal, detailed, CSV, ensemble, accessibility, profile,
  matrix, and spot-probability outputs.
- `app/main.cpp`: immutable configuration, deterministic multi-sequence task
  execution using `std::jthread`, and transactional result reduction.

The recurrence definitions, partition semantics, correctness findings, proof of
the cubic all-interval accessibility algorithm, performance measurements, and
the separate Legacy/Next audit are documented in
[MATHEMATICAL_AUDIT.md](MATHEMATICAL_AUDIT.md).

## Model scope

The base-pair (`--energy=B`) model is the strongest cross-version compatibility
oracle because it is independent of thermodynamic tables and uses `RT=1`. The
native nearest-neighbor model parses the public ViennaRNA v2 stack, loop,
mismatch, dangle, terminal, and enthalpy sections and applies temperature
interpolation without linking ViennaRNA. `--energyVRNA` accepts `Turner04`,
`Turner99`, `Andronescu07`, or an explicit `.par` path. Named sets are searched
in `INTARNANEW_PARAMETER_DIR`, `VIENNA_RNA_DATAPATH`, `VRNA_DATAPATH`, the active
Conda prefix, common system data directories, and bounded project-relative data
directories. Missing or incomplete tables are reported; they are never replaced
with approximate built-ins. All energy components remain visible in CSV output.
For builds outside this source tree, provide public `.par` files through one of
the search locations above. CMake test builds can set
`-DINTARNANEW_TEST_PARAMETER_DIR=/path/to/ViennaRNA`; a sibling legacy test
environment is detected only as a local development convenience.

Turner partition functions use the ordinary unsmoothed Boltzmann sum, equivalent
to ViennaRNA with `pf_smooth=0`. ViennaRNA's default smoothed dangle-2
partition can therefore differ slightly even when MFE energies agree.

Model P is an exact interaction partition only with exact prediction mode
(`--model=P --mode=M`) for the configured constraints. Heuristic mode H can
prune positive-weight paths and consequently underestimates Z. The
`IntaRNAens` personality selects model P but retains heuristic mode unless
`--mode=M` is also requested; disable the seed as well when an unconstrained
all-interaction ensemble is intended.

## Validation

`make check` runs real RNA–RNA interaction scenarios derived from the legacy
black-box integration suite—not tests of legacy internal classes. The suite covers
exact base-pair predictions, seed constraints, no-lonely-pair behavior, overlap
policies, coordinates, IUPAC input, and nearest-neighbor legacy-oracle motifs for
stacks, bulges, internal loops, exterior mismatches, temperature, and named sets.

On POSIX systems, `make check` also contains an executable-level 17-case parity
gate. It is skipped when no legacy executable was supplied. Set
`INTARNA_LEGACY_BIN` to a runnable IntaRNA Legacy binary to run both programs on
the same public `.parameter` corpus and compare stdout byte-for-byte; legacy logs
and diagnostics are routed to `/dev/null`. The standalone `make legacy-check`
target runs the same gate. CMake captures this environment variable at configure
time; it can also be set explicitly with
`-DINTARNA_LEGACY_BIN=/path/to/IntaRNA`. No legacy source or internal classes are
used.
