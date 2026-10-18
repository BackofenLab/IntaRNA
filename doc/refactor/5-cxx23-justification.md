# Is the C++23 requirement justified? — issue #244

Assessment of upstream `0a4568ad6e3da52f935221ff0da98a4a171e9833`,
21 September 2026. This is an assessment of the requirement's benefits,
not a proposal to downgrade the implementation.

**Conclusion:** C++23 has a concrete implementation rationale, principally
efficient construction of ViennaRNA constraint strings. The current source
really uses C++23, so changing only `configure.ac` would be incorrect.
However, the available evidence does not establish that mandatory C++23
delivers a material performance or correctness advantage over an earlier
baseline. Requiring it is defensible as a deliberate maintenance/toolchain
policy; it should not be presented as necessary for the scientific algorithms
or as a demonstrated speed improvement.

## Historical scope and recovery

This document preserves the September 2026 assessment of the revision named
above. Later changes, including the mdspan migration in PR #248, are outside
its scope. Source line numbers refer to that historical revision.

The shorter decision was published on [issue #244](https://github.com/BackofenLab/IntaRNA/issues/244),
but this full audit remained an untracked file in the local audit worktree.
The accompanying [archive](cxx23-audit/README.md) preserves the standalone
compiler probes and an abandoned C++17 compatibility experiment. The patch is
historical material for review and is not applied to the current build.

## Which benefits actually require C++23?

| Facility and locations | Standard | Benefit | Weight in the justification |
| --- | --- | --- | --- |
| `string::resize_and_overwrite` in `AccessibilityVrna.cpp:59` and `InteractionEnergyVrna.cpp:48` | C++23 | Fills owned string storage directly without first initializing bytes that will immediately be overwritten; the resulting string is null-terminated for ViennaRNA | The strongest C++23-specific reason: avoids redundant work with a standard API |
| `string::contains` in `RnaSequence.h:319`, `:503`, `:573` | C++23 | Expresses alphabet membership and ambiguity checks directly | Readability benefit; no inherent algorithmic advantage over `find(...) != npos` |
| String-stream `view()` in `OutputHandlerCsv.cpp`, `OutputHandlerEnsemble.cpp`, `OutputHandlerText.cpp` | C++20 | Inserts completed output without making an intermediate owning string | A real copy reduction, but cannot justify **23** specifically |
| Structured bindings and `unordered_map::try_emplace` in `PredictorMfeEns.cpp:121` | C++17 | Updates a partition entry with one lookup while preserving addition order | Useful optimization, independent of C++23 |
| `unique_ptr`, custom deleters, `make_unique`, direct sequence `c_str()` access | C++11/14 or earlier | Automatic cleanup, exception safety and removal of redundant sequence copies | Valuable modernization, but independent of C++23 |

The language/library classification is documented in the
[libstdc++ implementation tables](https://gcc.gnu.org/onlinedocs/libstdc++/manual/status.html).
The purpose and tradeoffs of direct string construction are described in
[WG21 P1072R10](https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2021/p1072r10.html).

The constraint helpers build buffers around ViennaRNA calls, rather than
changing predictor recurrences. Automatic ownership of those buffers is
useful, but ownership and exception safety do not themselves require
`resize_and_overwrite`. Its specific additional benefit is avoiding redundant
initialization while keeping a normal `std::string`.

No C++23-specific language construct was identified in the production source.
The requirement is driven by library APIs. The migration explicitly deferred
facilities such as `mdspan`, `format` and `print`; possible future uses do not
establish a present requirement. See the original
[migration record](2-c++23.md#implemented-c23-modernization).

## What the existing measurements establish

The [recorded modernization benchmark](2-c++23.md#substantive-c23-follow-up-benchmark)
compares `321364c56e914e9eff56ddba4846852297fe5c55` with
`7a9517b0587dd308659b61f81b26373ac7a5d0b5` on the same pinned
GCC 14.4.0 / Boost 1.85.0 / ViennaRNA 2.7.2 environment.
The evidence is also identified in
[PR #234](https://github.com/BackofenLab/IntaRNA/pull/234) and
[benchmark evidence commit 69174eb](https://github.com/BackofenLab/IntaRNA-refactor-benchmark/commit/69174eb00354e24427a6ee37c98d4363930114ca).

- Reported core timing ratios span 0.987x–1.022x (previous/modernized).
- The output-heavy case is 0.994x: 0.163162 s versus 0.164213 s.
- Canonical scientific outputs match, and peak RSS is near-stable.
- The recorded analysis explicitly finds no demonstrated general speedup.

These measurements support preserving behavior across the modernization.
They do **not** isolate a C++23 benefit: both revisions already request
C++23, and the candidate bundles string changes with ownership, sequence-copy
and partition-lookup changes. Nor do small mixed timings prove that the
C++23 facilities have zero benefit. There is no isolated end-to-end result
showing that `resize_and_overwrite` materially improves an IntaRNA workload.

The most useful missing evidence would compare the current implementation
with an otherwise identical constraint-buffer implementation, using the same
compiler and C++23 mode. That isolates the reason for the requirement without
confounding it with compiler versions or unrelated optimizations. It should
cover constraint-heavy workloads and representative normal use, and report
dispersion and output equality. A microbenchmark alone would not establish
an application-level benefit.

## Maintenance and compiler policy

Keeping one C++23 baseline allows the implementation to use the standard
buffer API directly without compatibility branches. The GCC 14 and Apple
Clang CI results establish a tested support matrix. Those are legitimate
maintenance considerations, but the same migration tests do not establish
that older standards are intrinsically unsafe or unsuitable.

There is also a distinction between using C++23 APIs and requiring a compiler
to advertise the final C++23 language value. The vendored Autoconf macro
rejects `__cplusplus < 202302L`. Locally, GCC 13.3 advertises `202100L` in
`-std=c++23` mode, but compiles and runs the **unchanged, exact library probe
from current `configure.ac`** with `-pedantic-errors`. The existing configure
check rejects that compiler before reaching its library probe. This is a
narrow diagnostic of the gate, not certification of GCC 13 for the complete
application. It shows why a compiler-version policy needs its own rationale.

## Lead-developer recommendation

Keep the current source and configure requirement unchanged for this audit.
Answer issue #244 with the actual dependency and its limited justification:

> C++23 is currently used intentionally: `resize_and_overwrite` avoids
> redundant initialization when constructing ViennaRNA constraint strings,
> and `contains` simplifies sequence checks. The stream-copy optimization
> requires only C++20, while ownership and partition-map improvements use
> older standards. Existing benchmarks preserve scientific output but do not
> demonstrate a material speedup attributable to C++23. Consequently, C++23
> is a chosen implementation/toolchain baseline, not a scientific necessity.
> If that baseline excludes users we intend to support, the present
> performance evidence is insufficient to insist on it.

I would defend the code's use of C++23 as purposeful. I would not defend the
hard requirement on correctness or demonstrated performance grounds with
the evidence currently available. No downgrade is included in this audit.
