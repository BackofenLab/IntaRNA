# IntaRNApvalue

Estimate the significance of an observed IntaRNA energy using a background of
**dinucleotide-preserving shuffled sequences**. This differs from
`IntaRNA_CSV_p-value.py`, which fits the energies of a genome-wide target screen.

## Installation and invocation

Use Python **3.10+**, IntaRNA on PATH (or beside the companion scripts / in the
checkout's `src/bin`), NumPy >= 1.23 and SciPy >= 1.9 for fitted distributions.
Scores-only and empirical modes need only the standard library and IntaRNA.
From the repository root:

```bash
python3 -m pip install ./python
python3 -m intarnapvalue --help
# Alternatively, without installing the package:
PYTHONPATH=python python3 -m intarnapvalue --help
```

`pip` installs NumPy and SciPy automatically. IntaRNA itself must be installed
separately. `--help` works without IntaRNA or scientific Python dependencies.
For an optional standalone executable, install PyInstaller and run
`python3 python/bin/build.py`; see `python/bin/README.md`.

## Sampling and interpretation

Each shuffle preserves length, endpoints, and mono/dinucleotide counts. Lowercase
and DNA input are normalized to uppercase RNA. FASTA inputs may wrap across lines
but must contain exactly one nonempty record; multiple records are rejected.
In `b` mode, each target and each query are independently shuffled, and IntaRNA
compares corresponding pairs in bounded batches, not all combinations.

`--samples` is the number of **reported interactions** collected. Pairs with no
reported interaction under the selected constraints are discarded and replaced.
The resulting p-value is therefore **conditional on a reported interaction**; it
does not estimate the unconditional probability over all shuffled pairs. The
number of discarded pairs is reported on stderr. The default attempt limit is
100 times `--samples`; reaching it is an error, with no partial score output.
Use `--max-attempts` to change that limit. An IntaRNA command failure stops
immediately. A parameter file is applied equally to the observed and shuffled
predictions; output settings enforce one energy per pair.

The fitted modes retain the historical models on **raw energies**: SciPy
`genextreme` (GEV), `gumbel_l` (left-tail Gumbel), or `norm` (Gaussian). The p-value
is the fitted CDF at the observed energy, replacing numerical PDF integration.
This preserves the former model, including its GEV orientation; the separate
CSV tool fits **negated** energies and has different statistical semantics.
At least three distinct finite energies are required for fitting. This is a
numerical minimum, not assurance of statistical adequacy. Inspect a sufficiently
large background for the intended application. Degenerate fits are errors.

`-d empirical` (alias `none`) reports the fraction of sampled energies less than
or equal to the observed energy, including ties. It accepts constant samples and
may return zero; its resolution is `1 / samples`, with no pseudocount correction.
Zero interaction energies remain valid. If the original pair has no reported
interaction, p-value mode fails; scores mode can still sample the background.

## Arguments

| Option | Default | Meaning |
| --- | --- | --- |
| `-q`, `--query` | Required | Literal RNA/DNA or a single-record FASTA path |
| `-t`, `--target` | Required | Literal RNA/DNA or a single-record FASTA path |
| `-s`, `--samples` | Required | Positive number of favorable shuffled pairs |
| `-m`, `--shuffle-mode` | Required | `q` (query), `t` (target), or `b` (both) |
| `-d`, `--distribution` | `gev` | `gev`, `gumbel`, `gauss`, `empirical`, or `none` |
| `-o`, `--output` | `pvalue` | `pvalue` or `scores` (one energy per line) |
| `-p`, `--parameterFile` | None | IntaRNA parameters used for every prediction |
| `--threads` | `0` | Nonnegative IntaRNA thread count; 0 uses available CPUs |
| `--randSeed` | None | Seed string; repeat with the same inputs/options/software |
| `--max-attempts` | `100 * samples` | Maximum pairs tried; must be at least samples |

```bash
python3 -m intarnapvalue \
  -q GCUGAAAAACAUAACCCAUAAAAUGCUAGCUGUACCAGGAACCA \
  -t GGUUUCUUCGCCUCUGCGUUCACCAAAGUGUUCACCC \
  -s 1000 -m b --threads 1 --randSeed example -d empirical
python3 -m intarnapvalue -q query.fasta -t target.fasta -s 1000 -m t -o scores > background.txt
```

A seed uses an instance-local RNG and does not reseed Python's global generator.
Corrected both-strand sampling changes results from earlier releases. Batch
results are restored to input order so thread scheduling does not reorder scores.
Diagnostics use stderr; successful stdout contains only the requested numbers.

## Python API and tests

Passing an argument list configures an instance without running predictions:

```python
from intarnapvalue.intarna_pvalue import IntaRNApvalue

analysis = IntaRNApvalue(['-q', 'GGGGGGGG', '-t', 'CCCCCCCC',
                         '-s', '10', '-m', 't', '--randSeed', 'example'])
observed = analysis.get_single_score(analysis.query, analysis.target)  # float or None
scores, discarded = analysis.get_scores()
if observed is not None:
    pvalue = analysis.calculate_pvalue_empirical(observed, scores)
```

The Python API raises exceptions for failed predictions and invalid data. CLI
argument errors use exit code 2; runtime errors use exit code 1.
From the repository root, run `python3 -m unittest discover -s python/tests -v`.
Tests cover both/single-strand sampling, sample counts, bounded failures,
parameter forwarding, RNG isolation, FASTA validation, zero energies, shuffle
composition, and agreement between CDFs and the former numerical integration.
