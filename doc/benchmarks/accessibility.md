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

## Example measurement (2026-09-29)

Linux x86-64, Ryzen 5 7530U, GCC 16.2.0 release build, Boost 1.85.0,
ViennaRNA 2.7.2, 100,000 bases, one folding/I/O thread. Medians of the three
alternating trials above (other build work was running on the machine):

| Format | Write (s) | Read (s) | Compressed bytes | ED cells differing from source | Maximum ED difference |
| --- | ---: | ---: | ---: | ---: | ---: |
| ED text + gzip | 8.726 | 3.365 | 20,891,786 | 551,861 | 1 (0.01 kcal/mol) |
| Binary `.agz` | 3.582 | 0.349 | 16,620,982 | 0 | 0 |

In this local run, binary output was about 2.4 times faster to write, 9.6 times
faster to read, and 20% smaller. Text differences arise in existing decimal
conversion; the benchmark does not modify that behavior. These figures describe
I/O of the same computed data, not an end-to-end genome-screen speedup.
