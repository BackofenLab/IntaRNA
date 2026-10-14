# ViennaRNA accessibility thread-safety audit

This audit accompanies [issue #265](https://github.com/BackofenLab/IntaRNA/issues/265).
The former critical section covered only `vrna_probs_window()`. That boundary
did not protect parameter construction: independent fold compounds can race
while ViennaRNA initializes them. Conversely, the sliding-window recurrence
uses compound-owned state and can run concurrently after synchronized setup.

## Reviewed upstream sources

The review covers the minimum supported ViennaRNA release, **2.4.14**, and the
**2.7.2** release used for the runtime experiments. Source links below are
pinned to immutable commits:

| Release | Commit | Scope |
| --- | --- | --- |
| 2.4.14 | [`0b8ca4c6604e7a8bea366af62d82f8f40cea49bf`](https://github.com/ViennaRNA/ViennaRNA/tree/0b8ca4c6604e7a8bea366af62d82f8f40cea49bf) | Minimum supported API and implementation |
| 2.7.2 | [`1ffec79f5e258896160f7362ced8263450f371dc`](https://github.com/ViennaRNA/ViennaRNA/tree/1ffec79f5e258896160f7362ced8263450f371dc) | Installed runtime and current parameter cache implementation |

This is a review of the paths IntaRNA calls, not a blanket claim that RNAlib,
its deprecated interfaces, or every future release is thread-safe. In
particular, compiling RNAlib with OpenMP does not make all its static state
thread-private.

## Call path and ownership

The relevant accessibility path is:

```text
AccessibilityVrna::fillByRNAplfold
  VrnaHandler::getModel -> vrna_md_copy
  vrna_fold_compound(PF | WINDOW)
    add_params -> vrna_params
    vrna_params_prepare -> vrna_exp_params
    sanitize_bp_span, sequence/constraint/matrix setup
  addConstraints
    vrna_constraints_add -> vrna_hc_add_from_db
    vrna_constraints_add_SHAPE, when requested
  vrna_fold_compound_prepare(PF | WINDOW)
  vrna_probs_window
    vrna_fold_compound_prepare(PF | WINDOW)
    allocate local auxiliary arrays
    update compound-owned hard/soft constraints and DP matrices
    callbackForStorage -> this accessibility object's ED matrix
  vrna_fold_compound_free
```

The explicit preparation before `vrna_probs_window()` is significant. The
window routine calls `vrna_fold_compound_prepare()` internally, and parameter
preparation can replace `exp_params` when its model differs from `params`.
Construction sanitizes span/window settings after initially constructing
parameters. Thus the synchronization rule must cover preparation as well as
the constructor. See the
[2.7.2 compound constructor and preparation](https://github.com/ViennaRNA/ViennaRNA/blob/1ffec79f5e258896160f7362ced8263450f371dc/src/ViennaRNA/fold_compound.c#L200-L276),
[preparation call chain](https://github.com/ViennaRNA/ViennaRNA/blob/1ffec79f5e258896160f7362ced8263450f371dc/src/ViennaRNA/fold_compound.c#L519-L626),
and [parameter replacement condition](https://github.com/ViennaRNA/ViennaRNA/blob/1ffec79f5e258896160f7362ced8263450f371dc/src/ViennaRNA/params/params.c#L427-L456).

After explicit preparation, the second preparation is idempotent with respect
to parameter construction for this path: the two parameter models match,
constraints do not change them, and no other worker accesses this compound.
Hard/soft constraint preparation and matrix scaling operate on its owned
storage. This reasoning assumes IntaRNA's single-sequence model, with no
external callbacks or custom grammar that mutate model settings. See
[hard constraints](https://github.com/ViennaRNA/ViennaRNA/blob/1ffec79f5e258896160f7362ced8263450f371dc/src/ViennaRNA/constraints/hard.c#L276-L310),
[soft constraints](https://github.com/ViennaRNA/ViennaRNA/blob/1ffec79f5e258896160f7362ced8263450f371dc/src/ViennaRNA/constraints/soft.c#L304-L356),
and [matrix preparation](https://github.com/ViennaRNA/ViennaRNA/blob/1ffec79f5e258896160f7362ced8263450f371dc/src/ViennaRNA/datastructures/dp_matrices.c#L299-L390).

## Shared state and required synchronization

| State or operation | Evidence | Consequence |
| --- | --- | --- |
| Parameter caches `p_pre`, `pf_pre`, and their initialization flags | [2.7.2 `params.c`](https://github.com/ViennaRNA/ViennaRNA/blob/1ffec79f5e258896160f7362ced8263450f371dc/src/ViennaRNA/params/params.c#L88-L217) enables `SPEEDUP_PARAMS`; both constructors rewrite cached model members and copy whole cached structures. Only `id` and `pf_id` are declared OpenMP thread-private. | Guard every IntaRNA entry to parameter construction with one shared mutex. Separate locks for accessibility and energy construction would not protect the same caches. |
| Parameter identifiers in older RNAlib | [2.4.14 declarations](https://github.com/ViennaRNA/ViennaRNA/blob/0b8ca4c6604e7a8bea366af62d82f8f40cea49bf/src/ViennaRNA/params/params.c#L88-L98), [identifier increment](https://github.com/ViennaRNA/ViennaRNA/blob/0b8ca4c6604e7a8bea366af62d82f8f40cea49bf/src/ViennaRNA/params/params.c#L552) | 2.4.14 lacks the newer caches, but its identifier is shared in builds without OpenMP. The same setup lock remains useful. |
| Global model defaults and legacy model variables | [2.7.2 `vrna_md_defaults_reset`](https://github.com/ViennaRNA/ViennaRNA/blob/1ffec79f5e258896160f7362ced8263450f371dc/src/ViennaRNA/model.c#L346-L432), [2.4.14 implementation](https://github.com/ViennaRNA/ViennaRNA/blob/0b8ca4c6604e7a8bea366af62d82f8f40cea49bf/src/ViennaRNA/model.c#L311-L384) | Passing a model to `vrna_md_defaults_reset()` changes process-wide defaults; it does not merely reset that local object. Remove the worker-side call when every fold already receives its explicit model. Guard parameter-file loading/default reads during handler setup. |
| `vrna_md_copy` | [2.7.2 implementation](https://github.com/ViennaRNA/ViennaRNA/blob/1ffec79f5e258896160f7362ced8263450f371dc/src/ViennaRNA/model.c#L196-L232) copies the supplied structure and its arrays. | Copying a stable handler model into worker-local storage needs no global lock. |
| Sliding-window DP, auxiliary arrays, and constraint updates | [2.7.2 window implementation](https://github.com/ViennaRNA/ViennaRNA/blob/1ffec79f5e258896160f7362ced8263450f371dc/src/ViennaRNA/partfunc/pf_window.c#L681-L958), [2.4.14 implementation](https://github.com/ViennaRNA/ViennaRNA/blob/0b8ca4c6604e7a8bea366af62d82f8f40cea49bf/src/ViennaRNA/LPfold.c#L669-L941) obtain arrays and models from the compound or allocate local helpers. | Independent prepared compounds can compute concurrently. Keep callback storage distinct and alive until each call returns. |
| Window compatibility state | [2.7.2 declarations](https://github.com/ViennaRNA/ViennaRNA/blob/1ffec79f5e258896160f7362ced8263450f371dc/src/ViennaRNA/partfunc/pf_window.c#L91-L110), [legacy wrapper](https://github.com/ViennaRNA/ViennaRNA/blob/1ffec79f5e258896160f7362ced8263450f371dc/src/ViennaRNA/partfunc/pf_window.c#L2170-L2243) | `backward_compat_compound` and `backward_compat` belong to the legacy wrapper path. The direct `vrna_probs_window()` call does not set them. Do not replace the direct API with a legacy wrapper under the same concurrency assumptions. |
| SHAPE D, Z, and W input | [2.7.2 implementation](https://github.com/ViennaRNA/ViennaRNA/blob/1ffec79f5e258896160f7362ced8263450f371dc/src/ViennaRNA/probing/SHAPE.c#L332-L404), [2.4.14 implementation](https://github.com/ViennaRNA/ViennaRNA/blob/0b8ca4c6604e7a8bea366af62d82f8f40cea49bf/src/ViennaRNA/constraints/SHAPE.c#L61-L120) | Setup loads input into the compound's soft constraints. W reads supplied pseudoenergies; this API does not perform stochastic optimization. Include constraint setup in the guarded preparation boundary. |
| Destruction | [2.7.2 `vrna_fold_compound_free`](https://github.com/ViennaRNA/ViennaRNA/blob/1ffec79f5e258896160f7362ced8263450f371dc/src/ViennaRNA/fold_compound.c#L117-L197) frees compound-owned arrays, constraints, and auxiliary data. | Ordinary accessibility destruction does not mutate the parameter caches. Each compound must have one owner. |
| Global PF compatibility pointer `pr` | [2.7.2 `vrna_pf`](https://github.com/ViennaRNA/ViennaRNA/blob/1ffec79f5e258896160f7362ced8263450f371dc/src/ViennaRNA/partfunc/partfunc.c#L208-L222), [2.4.14 implementation](https://github.com/ViennaRNA/ViennaRNA/blob/0b8ca4c6604e7a8bea366af62d82f8f40cea49bf/src/ViennaRNA/part_func.c#L129-L142) still write the process-global legacy pointer. | The accessibility conclusion does not justify unlocking global PF calls. Serialize the energy-model global-fold paths with the same mutex. |
| Lazy MFE SIMD dispatch | [2.7.2 dispatcher](https://github.com/ViennaRNA/ViennaRNA/blob/1ffec79f5e258896160f7362ced8263450f371dc/src/ViennaRNA/utils/higher_order_functions.c#L55-L128) replaces a shared function pointer on first use. | Keep global MFE under synchronization as well; this is a different path from the accessibility recurrence. |
| Diagnostics | [2.7.2 logger](https://github.com/ViennaRNA/ViennaRNA/blob/1ffec79f5e258896160f7362ced8263450f371dc/src/ViennaRNA/utils/log.c#L489-L552) lazily allocates callback storage; its default lock depends on the library's pthread build configuration. | OpenMP support alone is not proof for diagnostic paths. Install a dedicated logger lock through the public callback API before concurrent calls. |

## Boundary and limitations

The synchronization boundary must use a process-wide mutex shared by
`VrnaHandler`, `AccessibilityVrna`, and `InteractionEnergyVrna`; an OpenMP
critical section limited to one call site cannot cover the shared setup
state. RAII locking permits normal C++ exception unwinding and works when
IntaRNA itself is built without OpenMP.

The parallel boundary is **create, constrain, and prepare while
locked; compute the window recurrence with exclusive ownership of the
compound**. Keep global MFE/PF calls synchronized. RNAlib builds without
OpenMP use a conservative full-call fallback. Do not concurrently mutate
the handler model, reload global parameter tables through external RNAlib
calls, modify the same compound, or replace callback state while a call is
running. IntaRNA's mutex cannot coordinate direct RNAlib calls made by other
libraries in an embedding application.

This guarantee concerns the ViennaRNA boundary. Native-thread library users
must separately provide thread-safe IntaRNA/easylogging++ diagnostics, or leave
verbose/performance logging disabled. The CLI's OpenMP logging synchronization
is not a general native-thread logging contract.

For ViennaRNA 2.7 and later, handler initialization installs a process-lifetime
mutex callback once through `vrna_log_lock_set()`. This replaces the logger's
locking hook while preserving its output streams, levels, and message
callbacks. The logger mutex is distinct from the setup mutex because guarded
ViennaRNA operations can emit messages themselves. This protects diagnostics
even if RNAlib was built without pthread support. Embedding applications must
not reset or replace that hook while IntaRNA folds are running; callbacks must
not reenter folding. Earlier reviewed logging uses stdio directly rather than
the newer lazily initialized callback array. See the
[public locking callback](https://github.com/ViennaRNA/ViennaRNA/blob/1ffec79f5e258896160f7362ced8263450f371dc/src/ViennaRNA/utils/log.h#L61-L70)
and [hook setter](https://github.com/ViennaRNA/ViennaRNA/blob/1ffec79f5e258896160f7362ced8263450f371dc/src/ViennaRNA/utils/log.c#L284-L290).

The 2.4.14 source review establishes API availability and identifies its
state, but does not substitute for testing that release. Runtime comparisons
on 2.7.2 do not prove every intermediate release or build configuration safe.
Repeat the audit when changing the supported RNAlib range, adding model
mutations after preparation, changing SHAPE integration, or adding callbacks.

## Runtime evidence and reproduction

The diagnostic compares every unpaired probability returned by the callback
against a serial baseline, using independent compounds and callback storage.
It distinguishes fully concurrent lifecycles, serial creation followed by
parallel window calls, and synchronized preparation followed by parallel
window calls. Varying models matters: a shared cache can appear harmless when
all workers request identical parameters.

The committed [standalone diagnostic](vrna-window-threads.cpp) uses 24
deterministic sequences of varied lengths, windows, spans, and hard/soft
constraints. Its default profile also varies temperature and model flags.
`shared-model` fixes temperature to 37 degrees Celsius and uses the same
`noLP`/`noGUclosure` settings across cases, demonstrating that varying window
dimensions alone is sufficient to expose the cache race.

Build with a C++23 compiler, OpenMP, and RNAlib's pkg-config metadata available:

```bash
c++ -std=c++23 -O2 -fopenmp doc/analysis/vrna-window-threads.cpp \
  -o /tmp/vrna-window-threads $(pkg-config --cflags --libs RNAlib2)
/tmp/vrna-window-threads setup-locked 8 12
/tmp/vrna-window-threads setup-locked 8 30 shared-model
/tmp/vrna-window-threads window-only 8 12
```

The following modes intentionally exercise unsynchronized RNAlib state. Run
them as separate diagnostic processes, never as ordinary CI pass conditions:

```bash
/tmp/vrna-window-threads setup-only 8 30 shared-model
/tmp/vrna-window-threads window-locked 8 12
/tmp/vrna-window-threads unsafe 8 30 shared-model
```

Observed results used ViennaRNA 2.7.2's static library with OpenMP enabled,
GCC 14.4.0 (conda-forge build 14.4.0-7), and Linux x86-64 with CPU affinity
0-11. Runs requested eight workers; other build activity could share these
CPUs. No ViennaRNA 2.4.14 runtime was installed for the experiment.

| Mode/profile | Completed repetitions | Exact probability comparisons | Result |
| --- | ---: | ---: | --- |
| `setup-locked`, mixed model | 12 | 4,046,016 | No parameter, probability, or callback metadata mismatch |
| `setup-locked`, shared model | 30 | 10,115,040 | No parameter, probability, or callback metadata mismatch |
| `window-only`, mixed model | 12 | 4,046,016 | No parameter, probability, or callback metadata mismatch |
| `setup-only`, shared model | 1 | Not applicable | 17 compounds with mismatched parameter snapshots |
| `window-locked`, mixed model | 2 | 674,336 | 6,120 probability mismatches, one parameter mismatch, one metadata mismatch |
| `unsafe`, shared model | 1 | 337,168 | 6,120 probability mismatches, one parameter mismatch, one metadata mismatch |

A one-worker `unsafe` control completed three repetitions with 1,011,504
exact comparisons and no mismatches. Five fresh `setup-only` processes also
each reproduced parameter corruption. The failing modes stop at the first
iteration containing a mismatch, so repetition counts and mismatch counts
will vary with scheduling. These experiments locate a concrete failure
outside the original recurrence-only lock; finite passing stress runs are
not a race-detector proof or a performance measurement.

## Integrated regressions and existing SHAPE limitations

GCC 14 release and `--enable-debug` builds each passed `make tests -j4`:
all seven test programs/scripts passed, including 75,901 API assertions in
69 test cases. The source base was `c5f1d51` on `master` (the repository has no
`main` branch). `git diff --check` and shell syntax validation also passed.

`tests/AccessibilityVrna_test.cpp` compares all retained ED values (including
the extra band used for dangling ends), ES matrices, initiation energies and
both global ensemble energies against serial references. Twelve input/model
combinations vary sequence length, temperature, window, span, scaling and
hard constraints. Four native threads and four OpenMP workers each perform
48 jobs; assertions and exception reporting happen after workers join.

`tests/runAccessibilityThreads.sh` compares sorted predictions and every
per-sequence accessibility export at one and two workers, with SHAPE D and W
as well as unconstrained input. Unconstrained and D cases also request Eall
output, exercising overlapping target accessibility and global folding.

Two pre-existing serial failures in the installed ViennaRNA 2.7.2 path limit
SHAPE coverage. Z input crashes even without Eall output; W input crashes
when reported interactions request Eall. Both occur with the unmodified base
binary and the corrected binary at one thread. They are not introduced by
parallel window calculation, and are outside this synchronization change.
In particular, the existing `addConstraints()` passes window-specific options
also when called for global Eall, which is inappropriate for W's
`vrna_sc_set_up()` path. The regressions therefore exercise W without Eall and
do not claim working Z support.

A minimal serial Z reproduction, using either binary, is:

```bash
for i in $(seq 1 16); do printf '%s 0.5\n' "$i"; done > /tmp/reactivities.shape
IntaRNA --query=GGGGAAAACCCCUAGC --target=GGGGAAAACCCCUAGC --seedBP=3 \
  --qIntLenMax=8 --tIntLenMax=8 --qAccW=14 --tAccW=14 --qAccL=10 --tAccL=10 \
  --qShape=/tmp/reactivities.shape --tShape=/tmp/reactivities.shape \
  --qShapeMethod=Z --tShapeMethod=Z --outMode=C --outCsvCols=id1,id2,E,ED1,ED2 \
  --outMaxE=100 --threads=1 --default-log-file=/dev/null
```

An additional branch probe compiled the three changed implementation units
with `INTARNA_MULITHREADING=0` and without `-fopenmp`. Four native workers
matched serial ED, ES and Eall values for 48 jobs. A linker wrapper observed
four simultaneous `vrna_probs_window()` calls with the normal header and one
with an isolated header override removing `VRNA_WITH_OPENMP`. This verifies
both the parallel and serialized branches. It is **not** a runtime test of a
non-OpenMP RNAlib: the linked library remained the installed OpenMP-enabled
2.7.2, and unchanged IntaRNA objects came from the release archive.

## End-to-end measurement

The base and corrected binaries used the same GCC 14.4 release toolchain,
ViennaRNA 2.7.2 and options. The workload has 32 deterministic 400-nt queries,
a 6-nt target, and 120-nt accessibility windows. Kinetic interaction prediction
keeps this workload focused on accessibility. Every run exports all query ED
values as well as interaction structures, coordinates and energies.

On a Ryzen 5 7530U, affinity CPUs 8-11 exposes four logical CPUs on two physical
cores. One warmup preceded three alternating measured runs per configuration.
Other compilation activity ran on CPUs 0-7 during this first measurement, so
these numbers are workload observations rather than general speedup claims.

| Workers | Binary | Median wall time (range), s | Median peak RSS, KiB |
| ---: | --- | ---: | ---: |
| 1 | Base | 4.15 (4.11-5.07) | 18,376 |
| 1 | Corrected | 4.76 (4.51-5.92) | 18,392 |
| 4 | Base | 6.24 (4.71-6.48) | 30,396 |
| 4 | Corrected | 2.01 (1.98-2.03) | 30,564 |

The four-worker median improved by 3.10x with comparable peak memory. To
check the noisy single-worker result, four further alternating runs pinned
to CPU 8 after compilation finished measured 3.865 s (3.77-4.09) for the base
and 3.890 s (3.84-3.95) for the correction. That 0.025-s difference lies within
the observed run variation; it does not establish zero overhead.

All 24 runs matched the **base serial** reference exactly: sorted CSV output
(32 interaction rows) and all 32 ED exports (114,048 values per run). The
unmodified multiworker binary also happened to match on these inputs; its
passing benchmark does not negate the separate race reproducer.

Generate the inputs with Python 3:

```bash
mkdir -p /tmp/intarna-265-benchmark
python3 - <<'PY'
from pathlib import Path
from random import Random
out = Path('/tmp/intarna-265-benchmark')
rng = Random(265)
out.joinpath('queries.fa').write_text(''.join(
    '>query%02d\n%s\n' % (i + 1, ''.join(rng.choice('ACGU') for _ in range(400)))
    for i in range(32)))
out.joinpath('target.fa').write_text('>target\nGCGCGC\n')
PY
```

Use the following command with each binary, worker count, and a fresh output
directory for each repetition. Set `binary`, `threads`, and `out_dir` first;
adjust the CPU affinity to the local machine. The serial followup uses `8`
instead of `8-11`.

```bash
mkdir -p "$out_dir"
export OMP_DYNAMIC=FALSE OMP_PROC_BIND=TRUE OMP_PLACES=cores OPENBLAS_NUM_THREADS=1
/usr/bin/time -f '%e %U %S %M' -o "$out_dir/time.txt" \
  taskset -c 8-11 "$binary" \
  --query=/tmp/intarna-265-benchmark/queries.fa \
  --target=/tmp/intarna-265-benchmark/target.fa \
  --model=X --mode=K --seedBP=3 --qIntLenMax=8 --tIntLenMax=6 \
  --qAccW=120 --qAccL=120 --tAccW=6 --tAccL=6 --threads="$threads" \
  --outMode=C --outCsvCols=id1,start1,end1,id2,start2,end2,subseqDP,hybridDP,E,ED1,ED2 \
  --out="$out_dir/predictions.csv" --out="qAcc:$out_dir/query.ed"
```
