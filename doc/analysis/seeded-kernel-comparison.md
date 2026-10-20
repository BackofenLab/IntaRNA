# Seeded partition kernel comparison (PR #258, step 3)

Production uses the specialized stack-seed suffix automaton. Both candidates
passed independent exhaustive comparisons of raw Z, complete boundary partitions
and every pair numerator: all ten unrounded theory fixtures and 100 deterministic
heterogeneous-weight cases with missing edges, forbidden boundaries, mixed and
singleton occurrences, asymmetric spans and noLP. The tolerance is 2e-12 relative;
the oracle accumulates in long double. Mathematical fixture accessibility remains
unrounded, independently of native integer ED tests.

## Reproduce

GCC 14.4.0, Linux x86-64, release Autotools configuration, Z_type=double,
`-O3 -DNDEBUG`, Kokkos mdspan. The hidden Catch case contains every input and
uses deterministic weights and identical callback/checked-arithmetic contracts.

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /usr/bin/time -v \
  ./tests/runApiTests '[SeededBenchmark]'
```

Single-run milliseconds on the development host (these are kernel timings,
not native full-program speedups; rerun to account for shared-host noise):

| Case | Stack Z | Suffix Z | Stack Z+outside | Suffix Z+outside |
| --- | ---: | ---: | ---: | ---: |
| none | 0.000 | 0.000 | 0.009 | 0.000 |
| sparse | 59.462 | 43.934 | 156.972 | 140.072 |
| overlap | 109.494 | 77.083 | 304.623 | 249.862 |
| dense | 207.217 | 96.059 | 554.614 | 302.839 |
| long-stack | 1.003 | 0.879 | 3.057 | 1.466 |
| mixed | 196.146 | 99.959 | 522.997 | 310.235 |
| asymmetric | 65.849 | 37.285 | 182.743 | 121.412 |

The outside increment is the difference between each combined and partition-only
time. All corresponding Z values agree within 2e-12; pair masses were checked
before timing. The numerical kernels do not include native seed preprocessing.
Ordinary cases use 45 x 45 positions, spans 30 x 30, seed length four and loop
limits two. The long-stack case uses a 90-position diagonal (L up to 30); mixed
lengths are one through four; asymmetric spans are 12 x 24. No-seed input exits
before DP allocation. Sparse input has one seed. These are synthetic stress cases;
end-to-end native measurements belong to the final integration gate.

## Memory and decision

For a per-start dense rectangle of A cells and c=max(2, maximum seed length),
the suffix implementation stores (c+4)A partition values, or 2(c+4)A with the
outside pass, plus A bytes of validity and A size_t ending lengths. For c=4 and
A=900 this is 65,700 bytes partition-only and 123,300 with outside (eight-byte
Z_type/size_t). Maximal stacks use 4A or 8A partition values, A validity bytes
and up to 3L scratch partition values for the outside pass. The output matrix
adds n*m*sizeof(Z_type) for either kernel. Occurrences/eligible-start boxes are
separate O(number of retained seeds) metadata; neither candidate retains run
lists or boundary coefficients globally. These are physical rectangular-cell
counts, not just complementary vertices. Allocator overhead is not included.

The suffix candidate generally improves throughput, especially partition-only
and long-stack work. The extra bounded-span memory is acceptable for the intended
workload, so only suffix states ship in the library. Maximal stacks remain in a
test-only header for independent comparisons. This is not a claim of universal
speedup or lower memory use. Final integration also accounts for the global pair
matrix, private regional result, and legacy boundary-map storage when trackers
are enabled.
