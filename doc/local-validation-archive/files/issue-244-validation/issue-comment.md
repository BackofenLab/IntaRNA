**Decision: retain C++23 as the current implementation baseline; no change to `configure.ac`.**

The code uses `std::string::resize_and_overwrite` to construct ViennaRNA constraint strings directly in owned storage without redundant initialization, and `std::string::contains` for sequence checks. Changing only the configure requirement would not make the current source compatible with earlier standards.

This is an implementation and maintenance choice consistent with the [GCC 14 / Apple Clang support policy](https://github.com/BackofenLab/IntaRNA/discussions/232). It is not a scientific requirement: stream views need only C++20, and the ownership and partition-map improvements use earlier standards. The [recorded benchmarks](https://github.com/BackofenLab/IntaRNA/blob/0a4568ad6e3da52f935221ff0da98a4a171e9833/doc/refactor/2-c%2B%2B23.md#substantive-c23-follow-up-benchmark) preserve canonical outputs but do not demonstrate a material speedup attributable specifically to C++23.

A lower language baseline would require a separate compatibility change and validation.
