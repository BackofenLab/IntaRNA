# IntaRNAnew native companion tools

These utilities are C++23 implementations with no Python, R, Perl, shell, or
external-process runtime. Their reusable APIs live below
`include/intarnanew/tools/`. They fail on malformed data instead of silently
dropping rows or coercing values.

## `intarnanew-stats`

`fit` estimates Gaussian, Gumbel, or generalized-extreme-value parameters from
a numeric CSV column. Gaussian fitting is the closed-form population maximum
likelihood estimate. Gumbel and GEV fitting minimize negative log likelihood
with a deterministic multi-start Nelder-Mead search; GEV shape is constrained
to `[-1,1]`. The fit reports whether the numerical convergence tolerance was
reached.

`pvalue` appends the lower-tail probability `P(X <= E)`. This is the relevant
tail for IntaRNA energy because a more-negative energy is a better score.
`--empirical` instead uses the explicitly specified plus-one estimator

```
(number of background scores <= E + 1) / (background size + 1).
```

`adjust` implements Bonferroni, Holm, Hochberg, Benjamini-Hochberg, and
Benjamini-Yekutieli adjustment. Results remain in input order, and original row
order breaks equal-p-value ties deterministically.

Examples:

```
intarnanew-stats fit -i interactions.csv -c E -d gev
intarnanew-stats pvalue -i interactions.csv -o with-p.csv -c E -d gev
intarnanew-stats adjust -i with-p.csv -c p-value -m bh -o adjusted.csv
```

The API in `pvalue.hpp` adds deterministic mono- or exact directed-
dinucleotide-preserving shuffling and parallel score sampling. It takes a C++
score callback, so a prediction-library runner can be connected without
forking an `IntaRNA` executable. Sampling order and generated sequences do not
depend on thread count.

## `intarnanew-pvalue`

`intarnanew-pvalue` connects this API to the native predictor in-process. It
accepts literal RNAs or the first record of FASTA/gzip inputs, a positive sample
count (`--samples`, `--cardinality`, or `--scores`), shuffle side, optional
prediction parameter file, deterministic seed, distribution, and thread count.
`--output scores` prints the sampled energies in stable sample order;
`--output pvalue` prints the fitted lower-tail probability.

```
intarnanew-pvalue -q AACCGGUU -t GGUUCCAA -s 100 -m b \
  --randSeed 42 --threads 4 --output pvalue
```

## `intarnanew-csv`

This is a robust interaction-table parser and fusion utility. It supports CSV,
TSV, quoted separators, quoted newlines, and escaped quotes. Fusion produces a
stable union of the schemas in first-occurrence order, preserves row order, and
leaves fields empty when an input did not provide a column.

```
intarnanew-csv --source-column source -o all.csv first.csv second.csv
intarnanew-csv --deduplicate --separator tab a.tsv b.tsv
```

This is table fusion. It is not claimed to reproduce the historical
`IntaRNA-fuse.pl` biological construct workflow: its public help exposes RNA
folding and prediction parameters, but does not specify the construct sequence,
constraints, coordinate transformation, or energy recomputation sufficiently
for an independent scientific implementation.

## `intarnanew-svg`

This emits a self-contained SVG directly. `profile` plots any numeric CSV
columns and splits lines at `NA`; `heatmap` accepts IntaRNA `pMinE`-shaped data,
retains missing cells, and by default applies the documented convention of
clamping positive energies to zero. `regions` replaces the documented
`IntaRNA_plotRegions.R` workflow using `id1/start1/end1` or `id2/start2/end2`
columns. All labels and interactive cell titles are XML escaped.

```
intarnanew-svg profile target-minE.csv target.svg --x-column position --y-column minE
intarnanew-svg heatmap pair-minE.csv pair.svg
intarnanew-svg regions predictions.csv regions.svg --sequence 1
```

Output is deterministic for identical input and options. It intentionally
targets SVG only; rasterization is outside the standard C++ library.

## `intarnanew-mutate`

This is the prediction-independent CopomuS layer. Given base pairs, it
enumerates compensatory mutations using the documented `flip` or `any`
generators, supports GU/AU/CG wild-pair filters, validates explicit mutation
encodings, respects configured external coordinate origins, and materializes
the wild/mutant sequence combinations.

```
intarnanew-mutate -q GCAU -t CGUA --t-index-origin 10 --pairs '1&10' -g any
intarnanew-mutate -q GCAU -t CGUA --t-index-origin 10 \
  --mutation-encoding 'G1C&C10G'
```

CopomuS candidate selection (`mfe`, `mfeSO`) and ranking measures (`mfeCover`,
energy deltas, accessibility and energy profiles) depend on prediction results.
They are not assigned invented formulas here. The mutation API is designed to
feed a future in-process prediction runner for the `ww`, `wm`, `mw`, and `mm`
evaluations.

## Deliberately unspecified compatibility

Public documentation does not uniquely define:

- the optimizer, initialization, convergence policy, or failure policy used by
  historical SciPy distribution fitting;
- the exact random algorithm and seed stream used by the historical Python
  shuffle implementation;
- the biological construct and coordinate mapping used by
  `IntaRNA-fuse.pl`;
- ranking tie-breakers and every class boundary in historical CopomuS measures.

IntaRNAnew therefore publishes precise behavior above, tests it, and avoids
claiming byte-for-byte identity where the public contract cannot establish it.
