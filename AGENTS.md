# Working on IntaRNA

These instructions apply throughout this repository. Follow more specific
`AGENTS.md` instructions in a subdirectory if present. Keep this guide aligned
with the code and build files when their conventions change.

## Where to work

- `src/IntaRNA/`: the C++ library, including accessibility and energy models,
  predictors, seed/helix handlers, output handlers and prediction trackers.
- `src/bin/`: the executable and `CommandLineParsing` for CLI options and wiring.
- `tests/`: Catch-based API tests; `tests/data/`: CLI parameters and reference
  results exercised by `tests/runIntaRNA.sh`.
- `python/`, `perl/`, `R/`: companion tools; consult their local documentation
  when changing them.
- [README.md](README.md): user-facing installation, options and library usage.
  [doc/refactor/](doc/refactor/) records earlier analyses and measurements;
  its phase-specific branch instructions are historical, not the current
  contribution workflow.
- use '\n' (LF) as sole line ending in all text-based files.
  Convert CRLF to LF in editors or with `dos2unix` when needed.
- Use UTF-8 encoding for source and text files.

## Build and validation

Use the Autotools build. [configure.ac](configure.ac) defines requirements and
options; [.github/workflows/build.yml](.github/workflows/build.yml) is the
reference for the current platform checks: GCC 14 release/debug on Linux and
Apple Clang release on macOS.

C++23 support is required in both the compiler and its standard library.
`configure` checks `std::string::contains`, `std::string::resize_and_overwrite`
and `std::stringstream::view`. Do not lower the language standard to work around
an old toolchain. Dependencies include Boost, ViennaRNA (>= 2.4.14), zlib,
pkg-config and Autotools; OpenMP is used by the default multithreaded build.
[conda-build-env.yml](conda-build-env.yml) supplies library/build dependencies,
but does not select the C++ compiler. Use the platform setup in CI when needed.

Configure prefers native `std::mdspan` and otherwise uses the bundled Kokkos
headers. `INTARNA_USE_STD_MDSPAN` in the installed public configuration records
the choice; `--with-mdspan=std|kokkos` can select a backend explicitly. Preserve
the upstream headers and licenses in `src/mdspan` and `src/experimental` when
editing project code; their provenance is recorded in `doc/mdspan-storage.md`.

From the repository root, with dependencies available in standard locations:

```bash
bash autotools-init.sh
./configure
make -j2
make tests -j2
```

For dependencies in an activated conda environment, replace `./configure` with:

```bash
./configure --with-vrna="$CONDA_PREFIX" \
  --with-boost="$CONDA_PREFIX" --with-boost-libdir="$CONDA_PREFIX/lib" \
  --with-zlib="$CONDA_PREFIX"
```

Select the compiler with `CC`/`CXX` before configuring. Use `--enable-debug` for
debug checks and `--prefix` with a writable installation directory when testing
installation. CI includes platform-specific linker/OpenMP flags, installation,
standalone public-header compilation and an installed pkg-config consumer.
Consult those checks for build or public API changes.

- `make tests` (also `make test`) builds the program and runs both the API and
  CLI regression suites. Failure details are in `tests/test-suite.log` and the
  individual `tests/*.log` files. `make V=1` shows compiler/linker commands.
- For focused API iteration after building the library, use
  `make -C tests runApiTests`, then e.g. `./tests/runApiTests '[RNAsequence]'`.
  Run `make tests -j2` before submitting implementation changes; include debug
  validation when changing assertions, bounds or ownership.
- Add behavior/regression tests to the relevant `tests/*_test.cpp`; register
  new test sources in `tests/Makefile.am`. CLI regressions use paired
  `tests/data/*.parameter` and `*.testresult` files, automatically included by
  `tests/data/Makefile.am`; list other fixture types there as needed. Re-run
  the bootstrap after changing `Makefile.am`.
- Explain intentional reference-output changes and verify them independently;
  do not regenerate expected results merely to make a failure pass.
- For documentation-only changes, check paths, commands and the diff, and
  compile any C++ examples. A full C++ rebuild is unnecessary when library,
  executable and build behavior are unchanged. Report precisely which checks
  ran and any environmental blockers.

## C++ and header conventions

- Preserve the existing library families and virtual extension points. Use
  C++23 facilities where they help, while retaining support for the compiler
  and standard-library combinations covered by CI. Adding a newer library
  facility may require a configure feature check, not just `-std=c++23`.
- Match nearby formatting: tabs for C++ indentation, `IntaRNA` namespace,
  `INTARNA_..._H_` include guards, and existing class/member naming. Keep edits
  focused; do not reformat unrelated code or modernize vendored
  `src/easylogging++.*` or `tests/catch.hpp` as part of an ordinary fix.
- Keep class bodies focused on declarations. Define very short functions
  **explicitly `inline` after the class definition in the same header**, inside
  the namespace and include guard. Put nontrivial non-template implementations
  in the corresponding `.cpp`; keep template definitions visible in headers.
  [RnaSequence.h](src/IntaRNA/RnaSequence.h) illustrates the layout. For example:

```cpp
#ifndef INTARNA_EXAMPLE_H_
#define INTARNA_EXAMPLE_H_

namespace IntaRNA {

class Example {
public:
	/**
	 * Whether computation is complete.
	 * @return true if the result is complete
	 */
	bool
	isComplete() const;

private:
	//! whether computation is complete
	bool complete = false;
};

inline
bool
Example::isComplete() const
{
	return complete;
}

} // namespace IntaRNA

#endif
```

- Public headers must include what they use and compile independently when
  installed. Register new library headers/sources in `src/IntaRNA/Makefile.am`.
  Edit `configure.ac`, `Makefile.am` and source templates, not generated
  `configure`, `Makefile.in`, `Makefile` or configuration headers. Keep generated
  files, binaries and test output out of commits.

## Scientific correctness

- Preserve prediction semantics, CLI defaults and output formats unless the
  task explicitly changes them. Test energy values, interaction coordinates,
  traceback and relevant seed/helix/accessibility constraints, not just crashes.
- Internal sequence positions are zero-based; input/output indices can be
  shifted. Sequence 2 uses reversed accessibility in energy calculations.
  Reuse `RnaSequence`, `ReverseAccessibility` and offset-wrapper conversions.
- Use the energy/partition types, unit conversions, infinity sentinels and
  comparison helpers from [general.h](src/IntaRNA/general.h). Internal energies
  use hundredths of kcal/mol; do not mix them with kcal/mol or Boltzmann values.
- Preserve ownership and lifetimes when introducing RAII or non-owning views.
  Retain synchronization around shared state and ViennaRNA calls; existing
  accessibility calculations require serialization for thread safety.
- For recurrence or predictor changes, extend the relevant regression tests
  and short-sequence oracles (`PredictorTinyOracle_test.cpp` and
  `PredictorSeedOracle_test.cpp`). For optimizations, compare correctness first,
  then measure against the base with identical inputs, options, toolchain and
  thread count; record runtime/memory and any numerical differences.

## Documentation and ChangeLog

Document classes and public methods at their declarations using the existing
Doxygen `/** ... */` style, with a short purpose, `@param` and `@return` where
applicable. Explain non-obvious index direction, inclusive bounds, units,
constraints, ownership and exceptions. Use `//!` for brief member documentation
and implementation comments to explain recurrences and assumptions. Update
comments with the code; preserve existing author attribution. Doxygen setup is
in `doc/doxygen.cfg`; update `README.md` and CLI help for user-visible changes.

Update **[ChangeLog](ChangeLog)** (exact capitalization) for every introduced
change, including documentation, tests and build changes, in the same PR:

1. Add a concise bullet to the appropriate component/category in the opening
   "changes in development version since last release" summary.
2. Prepend a detailed entry below that summary's separator, before the previous
   dated entry, using `YYMMDD Contributor Name`. List affected files/classes with
   ` *`, and explain additions (`+`), modifications (`*`) or removals (`-`) below
   them. Describe the behavior/reason and link the issue or PR when available.
3. Retain existing summary bullets and dated history. Do not create a release
   section or change the package version unless that is part of the task.

## Branches and pull requests

Implement each independent change on its own descriptive branch from the
intended current base (normally `master`); use a separate worktree when another
task already occupies the checkout. Preserve unrelated local work. For dependent
changes, state the dependency and target the appropriate parent branch.

Open a PR for the implementation with a clear title and description covering
the target problem/issue, final changes, validation commands and results, and
any remaining limitations or compatibility effects. Use `Fixes #<issue>` only
when the PR resolves it fully. Include reproducible measurements for performance
claims. Review the final diff and run `git diff --check` before submitting.
