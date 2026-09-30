# Accessibility binary I/O benchmark

This benchmark measures two independent choices on identical ED data:

- Export stored matrix rows directly, or gather each row through `getED()` using
  the generic `Accessibility::writeBinary()` implementation.
- Write/read an uncompressed Boost binary archive, or add the production gzip
  filter used for `.agz` files.

It tests both `AccessibilityVrna` and a previously loaded
`AccessibilityFromStream`. Direct export passes non-owning matrix row views to
Boost without copying ED rows. Constrained accessibility data and other
implementations retain the generic path to preserve their `getED()` semantics.
All reads use the optimized loader, which fills retained rows in place and
validates discarded tails when the requested interaction length is smaller.
The version 1 archive layout is unchanged.

## Reproduction

Build and install IntaRNA, then compile against that installation:

```bash
export PKG_CONFIG_PATH=/path/to/intarna-install/lib/pkgconfig
c++ -std=c++23 -O3 doc/benchmarks/accessibility.cpp \
  $(pkg-config --cflags --libs IntaRNA) -leasylogging \
  -o accessibility-benchmark
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mkdir -p accessibility-benchmark-output/random accessibility-benchmark-output/ecoli
./accessibility-benchmark 100000 accessibility-benchmark-output/random > random.csv
curl -fL --retry 2 \
  'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000913.3&rettype=fasta&retmode=text&seq_start=1&seq_stop=100000' \
  -o ecoli-100k.fa
./accessibility-benchmark 100000 accessibility-benchmark-output/ecoli ecoli-100k.fa > ecoli.csv
```

Use the same compiler, dependency include/library paths and runtime library
paths as the IntaRNA build. With dependencies outside system paths, additional
`-L` and `-Wl,-rpath` flags may be necessary. The output directory must exist;
archive files are overwritten between trials.

Without a FASTA argument, the sequence is deterministic pseudorandom RNA (the
32-bit generator and seed 245 are in the source). With a FASTA argument, the
benchmark uses the requested number of bases from its first record. The second
measurement uses bases 1-100,000 of
[E. coli K-12 MG1655, NC_000913.3](https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3).
The downloaded FASTA has SHA-256
`7ad0766d4eb40e247ad5c2a9de2518738a39c8d261d27efb3bf1218c7525d0bc`.

Each dataset is folded once with ViennaRNA at 37 C, Turner04, folding window
150 and maximum interaction length 100. Three trials rotate through the eight
producer/method/format combinations. Output columns are
`trial,producer,serialization,format,write_seconds,read_seconds,bytes,different_cells`.
Every reloaded ED cell, including the dangling-end band, is compared exactly to
the original folded matrix; any difference aborts the run. Folding, the initial
stream-source load and verification are outside the timed sections.

Both formats use Boost file-descriptor streams and the production buffering
(512 KiB for output, default buffering for input); only the gzip filter differs.
Timings include compressor finalization and stream close, without `fsync` or
cache eviction. Reads therefore benefit from the filesystem cache. These are
local I/O wall times, not durable-storage throughput or end-to-end genome-screen
speedups. The short raw writes vary across trials; the CSV files retain that
variation rather than reporting only a best case.

## Measurements (2026-09-30)

Linux x86-64, Ryzen 5 7530U, GCC 16.2.0 release build, native `std::mdspan`, Boost
1.85.0 and ViennaRNA 2.7.2. One folding/I/O thread, datasets run sequentially,
with no concurrent IntaRNA builds/tests. Each input has 100,000 bases. The
following values are medians of three trials; all 48 reloads matched every ED
cell exactly. Per-trial results: [synthetic RNA](accessibility-random.csv) and
[E. coli](accessibility-ecoli.csv).

### Compression after direct matrix export

| Dataset | Producer | Format | Write (s) | Read (s) | Bytes |
| --- | --- | --- | ---: | ---: | ---: |
| Synthetic | Vrna | raw | 0.105 | 0.022 | 40,479,873 |
| Synthetic | Vrna | gzip | 2.796 | 0.248 | 16,620,982 |
| Synthetic | FromStream | raw | 0.102 | 0.021 | 40,479,873 |
| Synthetic | FromStream | gzip | 2.762 | 0.257 | 16,620,982 |
| E. coli | Vrna | raw | 0.108 | 0.024 | 40,479,873 |
| E. coli | Vrna | gzip | 3.019 | 0.271 | 16,850,736 |
| E. coli | FromStream | raw | 0.109 | 0.024 | 40,479,873 |
| E. coli | FromStream | gzip | 2.997 | 0.268 | 16,850,736 |

Gzip reduces the archive size by 59% on synthetic RNA and 58% on E. coli; raw
files are about 2.4 times as large. It adds roughly 2.7-2.9 seconds per write and
0.23-0.25 seconds per read in these measurements. Compression is therefore
retained for `.agz`: its CPU cost is substantial, but the uncompressed output is
also substantially larger. This follows the review criterion of removing gzip
only if the size saving is small. C++ callers can already use uncompressed
binary streams; this comparison does not introduce a new CLI format.

### Direct versus generic export

All times below are write medians. Both paths generate the same version 1
payload. Reads do not depend on which writer produced that payload; their
per-trial timings are included in the CSV files.

| Dataset | Producer | Raw direct (s) | Raw generic (s) | Gzip direct (s) | Gzip generic (s) |
| --- | --- | ---: | ---: | ---: | ---: |
| Synthetic | Vrna | 0.105 | 0.157 | 2.796 | 2.853 |
| Synthetic | FromStream | 0.102 | 0.150 | 2.762 | 3.054 |
| E. coli | Vrna | 0.108 | 0.384 | 3.019 | 3.016 |
| E. coli | FromStream | 0.109 | 0.162 | 2.997 | 3.052 |

Direct export removes the row scratch buffer and per-cell virtual `getED()`
lookups for unconstrained matrix-backed producers. The raw measurements expose
that benefit; gzip dominates compressed write time, where some differences are
within run-to-run variation. API tests independently check zero `getED()` calls
on direct export, byte-identical generic/direct archives, constraint masking,
and exact round trips. An archive from the original version 1 implementation
was also loaded and re-exported byte-for-byte identically.
