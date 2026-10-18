# Accessibility I/O benchmark

This benchmark isolates writing and reading the same accessibility matrix as
RNAplfold-style gzip-compressed ED text and as a binary `.agz` archive. It folds
one deterministic pseudorandom RNA with ViennaRNA (37 C, Turner04, folding window
150, interaction length 100), then alternates formats for three trials. Folding
and verification are outside the timed sections. Every loaded ED cell, including
the dangling-end band, is compared to the computed matrix. Binary input must
match exactly; text conversion differences are reported.

Build and install IntaRNA, then compile against that installation:

```bash
export PKG_CONFIG_PATH=/path/to/intarna-install/lib/pkgconfig
c++ -std=c++23 -O3 doc/benchmarks/accessibility.cpp \
  $(pkg-config --cflags --libs IntaRNA) -leasylogging \
  -o accessibility-benchmark
mkdir -p accessibility-benchmark-output
./accessibility-benchmark 100000 accessibility-benchmark-output
```

Use the same compiler, dependency include/library paths and runtime library
paths as the IntaRNA build. With dependencies outside system paths, additional
`-L` and `-Wl,-rpath` flags may be necessary. Output columns are
`trial,format,write_seconds,read_seconds,compressed_bytes,different_cells,max_ED_difference`.
ED differences are measured in internal units (hundredths of kcal/mol). Files are overwritten
between trials; the output directory must already exist. These timings include
compression and stream closing, but not an `fsync` to durable storage. Results
depend on sequence, bandwidth, compressor, filesystem and machine load.
