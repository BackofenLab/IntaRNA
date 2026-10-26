# Companion-tool audit (PR #298)

The branch builds on `feat/294-python-intarna-fuse` (PR #296). The R-to-Python
migration, shared infrastructure, CopomuS fixes, IntaRNApvalue fixes, and
packaging/validation are separate commits.

## Survey and common conventions

| Area | Review outcome |
| --- | --- |
| IntaRNA-fuse | Retain its fusion model, decimal energies and literal-sequence interface; share checked subprocess calls and executable lookup. |
| CSV p-values / region plots | Retain positional interfaces and the R-compatible formulas; share strict CSV header/row validation and atomic CSV publication. |
| CopomuS | Share RNA/FASTA validation and subprocess handling; fix explicit mutation dispatch, CG/GC filtering, signed indexing, zero/NA handling, ranking, and profile ownership. |
| IntaRNApvalue | Correct both-strand sampling, consistent constraints and failure handling; share FASTA and executable handling; expose empirical estimates and preserve fitted models. |
| Dinucleotide shuffling | Normalize at entry, handle length one, use an optional RNG, and consume adjacency lists with deques. |
| Packaging / executable builder | Declare Python 3.10 and real scientific dependencies; install shared modules with both pip and Autotools; replace shell-based PyInstaller calls and cleanup tied to the caller's directory. |

`intarna_tools.py` contains only standard-library I/O and validation. Scientific
calculations stay in their tools: the genome-wide CSV fit uses negated energies,
whereas IntaRNApvalue retains its historical raw-energy fitted distributions.
Combining these fits as one helper would change one tool's statistical model.
Likewise, fusion requires literal sequences for constraints, while CopomuS and
IntaRNApvalue can read a single FASTA record.

All entry points have side-effect-free imports and `main(argv)` interfaces.
Existing options remain; `GC` is an alias of CopomuS's advertised `CG` filter.
IntaRNApvalue adds `--max-attempts` and `-d empirical` (`none` alias). Sampled
p-values remain conditional on a reported interaction. Errors never enter
numeric/CSV stdout, and CopomuS/output-analysis CSV files are replaced only after
successful computation. See each tool's README for interpretation and examples.

## Correctness and performance checks

The Python unittest suite covers malformed/quoted CSV, zero/NA values, signed
coordinate round trips, explicit mutation bounds, all plot formats, independent
shuffle pairing, composition preservation, retry limits, subprocess errors,
parameter forwarding, and old-integration/new-CDF equivalence. Native checks
exercise both/single-strand sampling, deterministic output across thread counts,
profile combinations, explicit mutations and failure limits. The fusion suite
retains its independent synthetic energy oracle and optional native check.

CopomuS reuses the one wildtype prediction and immutable profile contents.
For N mutation candidates, measurement calls decrease from 4N to 3N+1, with
candidate-selection calls unchanged. A two-candidate regression verifies seven
calls and complete temporary-file cleanup. Filter membership uses sets without
changing candidate order. No native predictor code or defaults are changed.

Local comparison on Linux, Python 3.12.7, the same GCC 14 release IntaRNA binary,
ViennaRNA 2.7.2, and one IntaRNA thread (three runs each):

```bash
python3 CopomuS.py \
  -q GCUGAAAAACAUAACCCAUAAAAUGCUAGCUGUACCAGGAACCA \
  -t GGUUUCUUCGCCUCUGCGUUCACCAAAGUGUUCACCC --threads 1 -m E
```

Comparing pre-review commit `b4a9c99` with the revised sources, all columns,
mutation order, ranks and energy values match after normalizing line endings.
Wall times measured by `/usr/bin/time -f '%e %M'` were 2.90/2.94/2.97 seconds
before and 2.22/2.29/2.23 seconds after. Peak RSS ranges were 17,704–17,712 KiB
and 17,584–17,712 KiB. These are workload-specific measurements made during a
concurrent build, not a general throughput guarantee.

For the shuffle loop, use `s = 'AGGAUGGGGGACCUUCUCGUUCCU' * 4000`, reset
`random.seed(7)` before each call, and time `dinucl_shuffle(s)` with
`time.perf_counter()`. Three runs took 0.239–0.243 seconds before and 0.093–0.095
seconds after, with byte-identical shuffles. Deque consumption removes repeated
front deletion from Python lists; mononucleotide/dinucleotide counts and endpoint
preservation are tested separately. Both versions retain adjacency storage
linear in sequence length.

## Reproducing validation

```bash
python3 -m pip install -r python/requirements-output-tools.txt
python3 -m unittest discover -s python/tests -v
python3 tests/testIntaRNAFuse.py
bash autotools-init.sh
./configure # supply your compiler and dependency prefixes as needed
make tests -j2
make -C python check-python-tools
make install
make dist
```

Place a native IntaRNA on PATH to enable the optional native Python tests and
RNAfold on PATH for the fusion integration test. CI covers Python 3.10 with
minimum dependency versions, Python 3.13 with current dependencies, installed
entry points, and native execution in the C++ build matrix. PyInstaller itself
is optional; tests check its invocation and failure propagation without building
a platform-specific frozen executable.
