# Profiling evidence audit

Status: **PASS**

56 case definitions; 224/224 raw records; 224 validated records; 0 command/schema failures.

4 raw collection-time metadata records corrected from preserved stdout. Original raw records remain unchanged.

## Findings

- No integrity, coverage, or equality issues found.

## Output comparisons

| Comparison | Complete | Equal | Required | Performance comparison valid |
|---|---|---|---|---|
| thread_scaling | True | True | True | True |
| prediction_windows | True | True | False | True |
| memory_windows_model_x | True | True | False | True |
| memory_windows_model_s | True | True | False | True |
| cached_accessibility_original_observational | True | True | False | False |
| cached_accessibility_matched_caps | True | True | True | True |
| blas_environment | True | True | True | True |
| duplicate_default_case | True | True | True | True |

Original cached/computed timings are observational: their effective maximum interaction lengths differ by one. Only the explicit target/query 149/99 matched pair supports the cache performance comparison. The JSON companion separately checks base-pair-list output equality.


## Largest relative timing ranges

| Case | n | Median s | Min–max s | Range / median |
|---|---:|---:|---:|---:|
| blas_twelve | 3 | 0.086083 | 0.082933–0.131682 | 56.6% |
| blas_unset | 3 | 0.147197 | 0.081080–0.147391 | 45.0% |
| accessibility_cap75_control | 3 | 0.576970 | 0.572174–0.597119 | 4.3% |
| ensemble_H | 3 | 0.023965 | 0.023815–0.024850 | 4.3% |
| composition_gc20 | 3 | 0.199414 | 0.197479–0.203681 | 3.1% |
| no_pair | 3 | 0.040492 | 0.039312–0.040493 | 2.9% |
| gc_repeat | 3 | 2.981876 | 2.962701–3.048332 | 2.9% |
| throughput_12t | 3 | 0.973767 | 0.969737–0.997370 | 2.8% |
| window_600 | 3 | 2.874218 | 2.858498–2.934834 | 2.7% |
| algorithm_gapped_seed | 3 | 0.071494 | 0.070764–0.072644 | 2.6% |
| throughput_6t | 3 | 1.103543 | 1.102933–1.130404 | 2.5% |
| ensemble_M | 3 | 0.125607 | 0.125202–0.128245 | 2.4% |
| algorithm_heuristic_s | 3 | 0.095471 | 0.094646–0.096931 | 2.4% |
| cache_computed_control | 3 | 0.821475 | 0.806360–0.825383 | 2.3% |
| algorithm_default | 3 | 0.043563 | 0.042843–0.043846 | 2.3% |

## Coverage

| Case | Warmup | Measured | Successful measured | Timeouts | Stable outputs |
|---|---:|---:|---:|---:|---|
| bio_fhla_oxys | 1 | 3 | 3 | 0 | True |
| bio_phob_gcvb | 1 | 3 | 3 | 0 | True |
| bio_ilve_gcvb | 1 | 3 | 3 | 0 | True |
| length_t100_q100 | 1 | 3 | 3 | 0 | True |
| length_t300_q100 | 1 | 3 | 3 | 0 | True |
| length_t1000_q100 | 1 | 3 | 3 | 0 | True |
| length_t3000_q100 | 1 | 3 | 3 | 0 | True |
| length_t10000_q100 | 1 | 3 | 3 | 0 | True |
| length_t300_q50 | 1 | 3 | 3 | 0 | True |
| length_t300_q200 | 1 | 3 | 3 | 0 | True |
| length_t300_q400 | 1 | 3 | 3 | 0 | True |
| length_t300_q800 | 1 | 3 | 3 | 0 | True |
| composition_gc20 | 1 | 3 | 3 | 0 | True |
| composition_gc50 | 1 | 3 | 3 | 0 | True |
| composition_gc80 | 1 | 3 | 3 | 0 | True |
| au_repeat | 1 | 3 | 3 | 0 | True |
| gc_repeat | 1 | 3 | 3 | 0 | True |
| no_pair | 1 | 3 | 3 | 0 | True |
| algorithm_default | 1 | 3 | 3 | 0 | True |
| algorithm_exact_x | 1 | 3 | 3 | 0 | True |
| algorithm_heuristic_s | 1 | 3 | 3 | 0 | True |
| algorithm_exact_s | 1 | 3 | 3 | 0 | True |
| algorithm_helix | 1 | 3 | 3 | 0 | True |
| algorithm_no_seed | 1 | 3 | 3 | 0 | True |
| algorithm_gapped_seed | 1 | 3 | 3 | 0 | True |
| algorithm_seed_only | 1 | 3 | 3 | 0 | True |
| algorithm_no_lonely_pairs | 1 | 3 | 3 | 0 | True |
| algorithm_loop0 | 1 | 3 | 3 | 0 | True |
| algorithm_loop30 | 1 | 3 | 3 | 0 | True |
| algorithm_out10 | 1 | 3 | 3 | 0 | True |
| algorithm_out100 | 1 | 3 | 3 | 0 | True |
| ensemble_H | 1 | 3 | 3 | 0 | True |
| ensemble_M | 1 | 3 | 3 | 0 | True |
| accessibility_computed | 1 | 3 | 3 | 0 | True |
| accessibility_none | 1 | 3 | 3 | 0 | True |
| accessibility_global | 1 | 3 | 3 | 0 | True |
| accessibility_window75 | 1 | 3 | 3 | 0 | True |
| accessibility_cap75_control | 1 | 3 | 3 | 0 | True |
| window_300 | 1 | 3 | 3 | 0 | True |
| window_600 | 1 | 3 | 3 | 0 | True |
| throughput_1t | 1 | 3 | 3 | 0 | True |
| throughput_2t | 1 | 3 | 3 | 0 | True |
| throughput_4t | 1 | 3 | 3 | 0 | True |
| throughput_6t | 1 | 3 | 3 | 0 | True |
| throughput_12t | 1 | 3 | 3 | 0 | True |
| memory_t3000_q800 | 1 | 3 | 3 | 0 | True |
| memory_window300 | 1 | 3 | 3 | 0 | True |
| cache_computed_control | 1 | 3 | 3 | 0 | True |
| cache_reused | 1 | 3 | 3 | 0 | True |
| blas_one_control | 1 | 3 | 3 | 0 | True |
| blas_unset | 1 | 3 | 3 | 0 | True |
| blas_twelve | 1 | 3 | 3 | 0 | True |
| memory_s | 1 | 3 | 3 | 0 | True |
| memory_s_window300 | 1 | 3 | 3 | 0 | True |
| cache_matched_computed | 1 | 3 | 3 | 0 | True |
| cache_matched_reused | 1 | 3 | 3 | 0 | True |

The JSON companion contains detailed provenance, corrections, scatter, and file checks.
