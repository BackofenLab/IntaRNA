# C++17 compiler baseline (issue #244)

The audit of upstream `0a4568ad6e3da52f935221ff0da98a4a171e9833` found that
lowering only `configure.ac` would fail: the source used two C++23 library
facilities and one C++20 facility. None is required by the prediction
algorithms. This change replaces them and makes strict C++17 the supported
baseline for IntaRNA and its library consumers.

| Source | Previous requirement | C++17-compatible replacement |
| --- | --- | --- |
| `RnaSequence.h` | C++23 `string::contains` | `find` compared with `string::npos` |
| `AccessibilityVrna.cpp`, `InteractionEnergyVrna.cpp` | C++23 `string::resize_and_overwrite` | Sized, owned strings filled by index |
| `OutputHandlerCsv.cpp`, `OutputHandlerEnsemble.cpp`, `OutputHandlerText.cpp` | C++20 string-stream `view` | String-stream `str` |
| `PredictorMfeEns.cpp` | C++17 structured bindings and `unordered_map::try_emplace` | Retained unchanged |

C++14 is insufficient for the retained partition accumulation. Keeping
`try_emplace` preserves its single lookup and floating-point addition order.
The scoped ownership, direct `c_str()` calls into ViennaRNA, and other
prediction optimizations also remain in place.

The compatibility tradeoff is an initialization pass over constraint buffers
and temporary string copies when inserting completed output. Constraints
remain owned and null-terminated; every position is filled before ViennaRNA
uses them. Output uses the same string insertion semantics and the same
OpenMP critical sections. This change makes no performance-equivalence claim.

`configure.ac` requires strict C++17 and separately probes the standard
library's `unordered_map::try_emplace`, so a mismatched compiler/library pair
fails at configure time. The vendored Autoconf macro already supports C++17
and needs no changes. Installed pkg-config metadata continues to leave the
language mode to the consumer.

CI builds with GCC 11, GCC 14 (release and debug), and Apple Clang. Each job
uses C++17 with `-pedantic-errors`, runs the API and CLI suites, installs the
library, compiles every installed public header independently, and links and
runs an external pkg-config consumer.
