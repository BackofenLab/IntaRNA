# Mathematical audit of IntaRNA recurrences and partition functions

Audit date: 2026-08-13

## 1. Scope and evidence

This report keeps the three implementations in this workspace distinct:

- **Legacy** means the historical implementation in the [repository root](..).
- **Next** means a separately audited successor that retains much of the legacy predictor architecture; its source is not part of this repository.
- **IntaRNAnew** means this independent C++23 implementation. It is the only source tree modified by this audit.

A finding is labelled as follows:

- **reproduced**: a minimal executable case demonstrates the defect;
- **source-confirmed**: the recurrence or data flow proves the defect;
- **fixed with regression**: IntaRNAnew was changed and a focused oracle now protects the result;
- **proposed**: mathematically justified work that has not been implemented.

The audit covered monomer folding, interval-unpaired probabilities, interaction-path dynamic programming, site and global interaction partition functions, heuristic pruning, seed/noLP/noGU constraints, numerical stability, and result materialization.

## 2. Executive result

The IntaRNAnew recurrences have a sound decomposition, but the audit found and fixed four scientific defects:

1. both asymmetric Turner 1x2/2x1 internal-loop nucleotide orders were wrong in monomer and duplex scoring;
2. model-P interaction partition weights omitted terminal-pair and dangle energies;
3. valid Turner99 and Andronescu07 parameter files were rejected because their six-value Misc sections were parsed as exactly four values;
4. CSV P_E was reconstructed from a centikcal-truncated display energy and could exceed one instead of reporting the normalized site probability.

The main algorithmic breakthrough is an inside/outside theorem that generates every interval-unpaired partition for the scalar noncrossing model in \(O(n^3)\) time and \(O(n^2)\) memory. The previous all-interval path performed \(O(n^2)\) separate \(O(n^3)\) constrained folds, hence \(O(n^5)\) worst-case work. Queries are now \(O(1)\) after construction.

The performance investigation then found that the generalized interaction engine was doing path and ensemble work that the requested output could not observe. An output-aware scalar max-plus kernel now handles exact, unseeded, base-pair/model-S prediction when only ranked boundaries and energies are required; all broader configurations retain the full predictor.

On the three parity-gated biological workloads, the optimized executable is faster than Legacy in every case:

| Case | Legacy median | IntaRNAnew before | IntaRNAnew final | Final speedup | Final RSS |
|---|---:|---:|---:|---:|---:|
| fhlA/OxyS | 504.936 ms | 3409.157 ms | 36.066 ms | 14.000x | 5,120 KiB |
| ilvE/GcvB-ST | 4008.653 ms | 30093.911 ms | 208.473 ms | 19.229x | 5,120 KiB |
| phoB/GcvB | 3821.776 ms | 29566.255 ms | 205.155 ms | 18.629x | 5,120 KiB |

The geometric-mean speedup is **17.117x over Legacy** and **125.285x over the pre-redesign IntaRNAnew**. Section 10 gives the protocol, root cause, and memory evidence.

## 3. Notation and ensemble definitions

Internal indices are zero-based. The public inclusive interval \([l,r]\) corresponds to the mathematical half-open interval \([l,r+1)\).

For energy \(E\), temperature \(T\), and gas constant \(R\), define the Boltzmann map

\[
B(E)=\exp[-E/(RT)].
\]

The implementation stores partition values as logarithms and evaluates every addition with log-sum-exp:

\[
\operatorname{LSE}(x,y)
=\max(x,y)+\log\!\left(1+\exp[-|x-y|]\right).
\]

Thus zero weight is represented by \(-\infty\), multiplication becomes addition, and no explicit partition scaling is necessary.

A monomer partition includes the empty structure:

\[
Z_{\mathrm{mono}}=\sum_{\sigma\in\mathcal S}B(E(\sigma)),
\qquad
G_{\mathrm{mono}}=-RT\log Z_{\mathrm{mono}}.
\]

An interaction partition does **not** add a null/no-interaction state. It sums only sites and paths admitted by the configured domains, seed, energy, accessibility, noLP/noGU, length, and output filters. Reported site probabilities are conditional on this admitted interaction ensemble.

## 4. IntaRNAnew monomer recurrences

### 4.1 Scalar additive-pair model

Let \(Q[a,b)\) be the partition of noncrossing matchings on \([a,b)\), with

\[
Q[a,a)=1.
\]

Let \(u_j\in\{0,1\}\) state whether nucleotide \(j\) may remain unpaired, and let \(W_{ij}\ge 0\) be the Boltzmann weight of pair \((i,j)\). Forbidden pairs have \(W_{ij}=0\). Pair compatibility, minimum hairpin distance, maximum span, and x/b constraints are all encoded in \(W\); a forced-paired p position sets its unpaired indicator to zero.

Classifying a structure by the status of its rightmost nucleotide gives

\[
Q[a,b)
=
u_{b-1}Q[a,b-1)
+
\sum_{i=a}^{b-2}
Q[a,i)\,W_{i,b-1}\,Q[i+1,b-1).
\tag{1}
\]

The first term contains structures where \(b-1\) is unpaired. In every structure in the second term, \(i\) is the unique partner of \(b-1\); noncrossing then separates the prefix and enclosed interval. These classes are disjoint and exhaustive, so Equation (1) neither omits nor double-counts a structure.

IntaRNAnew uses this recurrence for NativeAccessibility and for base-pair folding without noLP/noGU-end state dependencies. Base-pair folding uses \(\log W_{ij}=1\), corresponding to an energy of \(-1\) at \(RT=1\); NativeAccessibility uses its sequence-dependent pair energy.

### 4.2 noLP and noGU-end states

Pair admissibility is not scalar under noLP or noGU-at-helix-ends. IntaRNAnew therefore retains two paired states:

- \(B^{A}_{ij}\): any structure enclosed by pair \((i,j)\);
- \(B^{S}_{ij}\): the subset whose immediate inner pair is \((i+1,j-1)\).

A pair needing right-stack support selects \(B^S\). An outer stack can instead certify its inner pair, which selects \(B^A\). A GU pair at a helix end is excluded unless both outer and inner stack context prove it is internal.

This stateful recurrence is source-reviewed as correct. It deliberately remains the fallback for per-interval constrained folds because substituting a scalar pair token would lose the stack context.

### 4.3 Turner nearest-neighbour folding

Write \(\ell(E)=-E/(RT)\). For a closing pair \((i,j)\), the Turner recurrence has the schematic form

\[
B^A_{ij}
=
\operatorname{LSE}\left(
  H_{ij},
  \{\,\ell(E_{\mathrm{int}}(i,j;k,l))+B^{\mathrm{sel}}_{kl}\,\}_{k,l},
  \ell(E_{\mathrm{multi-close}}(i,j))+M^{\ge2}[i+1,j)
\right),
\tag{2}
\]

where \(H_{ij}\) is the hairpin term and \(B^{\mathrm{sel}}\) chooses the state required by noLP/noGU context. The stack-only state is

\[
B^S_{ij}
=
\ell(E_{\mathrm{stack}}(i,j;i+1,j-1))
+
B^{\mathrm{sel}}_{i+1,j-1}.
\tag{3}
\]

The multiloop interval states \(M^0,M^1,M^{\ge2}\) record exactly zero, exactly one, or at least two branches. With \(U_a\) the multiloop-unpaired term and \(C_{a,r}\) a valid stem branch,

\[
M^0[a,b)=U_a+M^0[a+1,b),
\]

\[
M^1[a,b)
=
\operatorname{LSE}\left(
U_a+M^1[a+1,b),
\{\,C_{a,r}+M^0[r+1,b)\,\}_r
\right),
\]

\[
M^{\ge2}[a,b)
=
\operatorname{LSE}\left(
U_a+M^{\ge2}[a+1,b),
\{\,C_{a,r}+\operatorname{LSE}(M^1,M^{\ge2})[r+1,b)\,\}_r
\right).
\tag{4}
\]

Consequently, a multiloop closure in Equation (2) contains at least two branches. The exterior recurrence has a unique left-to-right tokenization into an unpaired nucleotide or a complete stem:

\[
X[r+1]
\mathrel{\oplus}=
X[r]+\ell(E_{\mathrm{unpaired}}(r)),
\]

\[
X[j+1]
\mathrel{\oplus}=
X[i]+B^{\mathrm{sel}}_{ij}+\ell(E_{\mathrm{exterior}}(i,j)).
\tag{5}
\]

Equations (2)--(5) are structurally sound. The implementation represents Vienna dangle-2 as direct-neighbour decoration and matches the ordinary unsmoothed Vienna partition convention, pf_smooth=0. ViennaRNA's default smoothed partition may therefore differ slightly even when MFE energies agree.

## 5. All-interval unpaired-partition theorem

### 5.1 Identity

Define the pair token

\[
R(i,j)=W_{ij}Q[i+1,j)
\]

and its outside context in ordinary arithmetic

\[
O_R(i,j)
=
\frac{\partial Q[0,n)}{\partial R(i,j)}.
\]

Let \(U=[a,b)\) be required to remain unpaired. With unit allowed-unpaired weights, its constrained partition is

\[
Z_U
=
Q[0,a)Q[b,n)
+
\sum_{\substack{i<a\\j\ge b}}
O_R(i,j)\,W_{ij}\,Q[i+1,a)\,Q[b,j).
\tag{6}
\]

The desired probability and opening energy are

\[
P(U\ \mathrm{unpaired})=\frac{Z_U}{Q[0,n)},
\qquad
ED(U)=-RT\log P(U\ \mathrm{unpaired}).
\tag{7}
\]

An interval containing a forced-paired p position has \(Z_U=0\) and is rejected before Equation (6).

### 5.2 Proof

Every noncrossing structure in which \(U\) is unpaired belongs to exactly one of two classes.

1. **No pair spans \(U\).** The structure separates at the interval boundaries, yielding
   \(Q[0,a)Q[b,n)\).

2. **At least one pair spans \(U\).** Among the spanning pairs there is a unique innermost pair \((i,j)\) with \(i<a\) and \(j\ge b\). Removing its pair token leaves the outside context \(O_R(i,j)\). Because this pair is innermost with respect to \(U\), its residual interior cannot contain another spanning pair and factorizes into \(Q[i+1,a)Q[b,j)\). Restoring \(W_{ij}\) gives the summand in Equation (6).

The classification is exhaustive, and uniqueness of the innermost spanning pair prevents double counting.

### 5.3 Cubic evaluation

The spanning sum is two triangular products over the log-sum-exp semiring. Define

\[
L[a,i]=Q[i+1,a),\quad
C[i,j]=O_R(i,j)W_{ij},\quad
R[j,b]=Q[b,j).
\]

Then the spanning contribution for every interval is

\[
S[a,b]=(L\otimes C\otimes R)[a,b].
\]

IntaRNAnew evaluates this as a right-end contraction followed by a left-end contraction. The inside pass, reverse outside pass, and each contraction are \(O(n^3)\); all tables occupy \(O(n^2)\). Every subsequent interval query is \(O(1)\).

Previously a single uncached interval query cost \(O(n^3)\). Generating all \(O(n^2)\) intervals therefore cost \(O(n^5)\) in the worst case. The new constructor costs \(O(n^3)\) total.

### 5.4 Scope

Equation (6) is exact for additive-pair noncrossing matchings with unit weights for allowed unpaired positions. It is used by:

- NativeAccessibility; and
- base-pair folding when noLP=false and noGUend=false.

It is **not** directly valid for Turner loop energies, noLP, or noGU-end, where a pair token has grammar state. Those cases retain the existing stateful/per-interval recurrence. A future extension must carry pair, stack, multiloop, and outside state indices.

## 6. IntaRNAnew interaction-path DP

An interaction path is an antiparallel sequence of base pairs

\[
p=((t_0,q_0),\ldots,(t_m,q_m)),
\qquad
t_0<\cdots<t_m,\quad q_0>\cdots>q_m.
\]

The core energy is

\[
E_{\mathrm{core}}(p)
=
E_{\mathrm{init}}
+
\sum_{k=1}^{m}
E_{\mathrm{loop}}\big((t_{k-1},q_{k-1}),(t_k,q_k)\big).
\tag{8}
\]

A DP cell fixes the terminal pair. A state additionally records the starting pair, seed automaton state, the noLP left-stack bit, and any path suffix required by explicit seeds or helix-block validation. Two paths with the same cell and state key have the same future transition set and the same site-dependent exterior energy. They may therefore be merged by

\[
\log W_{\mathrm{merged}}
=
\operatorname{LSE}(\log W_1,\log W_2),
\tag{9}
\]

while retaining the lower-energy path as the deterministic representative.

For site \(s\), let \(\mathcal P_s\) be its admitted paths and define

\[
E_{\mathrm{ext}}(s)
=
E_{\mathrm{end,L}}+E_{\mathrm{end,R}}
+E_{\mathrm{dangle,L}}+E_{\mathrm{dangle,R}}
+ED_T(s)+ED_Q(s)+E_{\mathrm{add}}.
\tag{10}
\]

Model P now computes

\[
Z_s
=
\exp[-E_{\mathrm{ext}}(s)/(RT)]
\sum_{p\in\mathcal P_s}
\exp[-E_{\mathrm{core}}(p)/(RT)],
\tag{11}
\]

\[
Z_I=\sum_s Z_s,\qquad
G_I=-RT\log Z_I,\qquad
P(s)=Z_s/Z_I.
\tag{12}
\]

The DP state's log weight contains Equation (8), so every term of Equation (10) must be applied exactly once at site completion. The pre-audit code omitted the four end/dangle terms.

Exact mode M performs no state-count pruning. Heuristic mode H retains protected and boundary-best states plus up to 96 additional candidates; 96 is not a strict total-state cap. Whenever pruning removes a positive path weight in model P,

\[
Z_H\le Z_M,\qquad G_H\ge G_M.
\tag{13}
\]

Thus model P is exact only in mode M for the represented constraints. IntaRNAens selects model P but otherwise inherits heuristic mode and default seed constraints; use mode M for an exact configured ensemble, and disable the seed too if an unconstrained all-interaction ensemble is intended.

## 7. Model and mode semantics

Model and prediction mode are independent; notably, model S and mode S do not mean the same thing.

| Model | Site contribution | Global partition |
|---|---|---|
| S, single-site | one MFE path per site | sum of site-MFE Boltzmann weights |
| X, seed extension | one MFE extension per site | valid like S without a required seed; deliberately unavailable for seeded X |
| B, helix blocks | retained paths satisfying block decomposition | heuristic retained-path ensemble; approximate |
| P, ensemble | sum of retained path weights per site | exact in mode M; underestimated after H pruning |

The accumulated log weight is not used to select or report representatives for models S and X. The selected-only scalar kernel avoids that work entirely for its supported model-S workload; the generic engine still retains it because the weight remains essential for model P and partitioned helix-block calculations.

CSV P_E now reports the normalized stored \(P(s)\). Before this audit it recomputed

\[
\exp[-\operatorname{trunc}_{0.01}(E_s)/(RT)-\log Z_I],
\]

which could exceed one and did not necessarily sum to one.

## 8. IntaRNAnew defects fixed

| Finding | Evidence before fix | Correction and regression |
|---|---|---|
| Turner int21 axes were permuted for both 1x2 and 2x1 orientations in monomer folding and duplex energy | RNAeval: AACGCAGCGU / pxxpxxxpxp should be 8.90, old result 7.80; mirrored AAUCGUAAGU case should be 10.10, old 9.00. Duplex UCGU&AAA and AAU&AAGU should have loop 3.70/full 8.80, old full 7.70. | Corrected nucleotide order in folding.cpp and energy.cpp; forced-structure and duplex tests added. |
| Model-P weights omitted terminal AU/end and dangle terms | One-pair A/U has \(E=4.10+0.50+0.50=5.10\), but old site \(G\) and log Z used 4.10. | Added both end and dangle sides at site completion; seed-only core normalization changed to init+loops to avoid double counting; single-AU and dangle-domain tests added. |
| Six-value Misc sections were rejected | Named Turner99 and Andronescu07 construction threw “Misc must contain 4 values.” | Parser accepts four or six values and reads LXC when present; exact GAAAC partition tests use hairpin 5.70 and 4.75 respectively. |
| CSV P_E used rounded display energy | A one-site AU ensemble could report 1.01636 although its normalized probability was exactly 1. | Text formatting and typed sorting now use Interaction::probability; a regression makes energy order disagree with probability order and verifies probability sorting. |

After the int21 correction, 700 randomized dangle-0 folds agreed with Vienna's partition oracle within approximately \(4.3\times10^{-7}\) kcal/mol. Dangle-2 agreed when Vienna was configured with pf_smooth=0.

## 9. Legacy and Next recurrence audit

### 9.1 Ensemble interaction recurrence

For fixed right boundary \(\mathbf j=(j_1,j_2)\), let \(H_{\mathbf j}(\mathbf i)\) be the hybrid partition beginning at \(\mathbf i=(i_1,i_2)\). Without noLP,

\[
H_{\mathbf j}(\mathbf j)=B(E_{\mathrm{init}}),
\]

\[
H_{\mathbf j}(\mathbf i)
=
\sum_{\mathbf k\in K(\mathbf i,\mathbf j)}
B(E_{\mathrm{loop}}(\mathbf i,\mathbf k))
H_{\mathbf j}(\mathbf k).
\tag{14}
\]

The admissible set \(K\) enforces complementarity, internal-loop limits, boundaries, and the fixed right endpoint. The completed site weight multiplies Equation (14) by all exterior site factors, including accessibility, terminal/end, dangle, and additive terms.

For noLP, a single pair is invalid, so \(H_{\mathbf j}(\mathbf j)\neq B(E_{\mathrm{init}})\). A clean two-state statement is:

- \(H_{\mathbf j}(\mathbf i)\): the leftmost pair still needs a right stack;
- \(A_{\mathbf j}(\mathbf p)\): pair \(\mathbf p\) already has left-stack support.

Then

\[
H_{\mathbf j}(\mathbf i)
=
B(E_{\mathrm{stack}}(\mathbf i,\mathbf i+\mathbf1))
A_{\mathbf j}(\mathbf i+\mathbf1),
\tag{15}
\]

\[
A_{\mathbf j}(\mathbf p)
=
\mathbf1_{\mathbf p=\mathbf j}B(E_{\mathrm{init}})
+
H_{\mathbf j}(\mathbf p)
+
\sum_{\mathbf k\in L(\mathbf p,\mathbf j)}
B(E_{\mathrm{loop}}(\mathbf p,\mathbf k))
H_{\mathbf j}(\mathbf k),
\tag{16}
\]

where \(L\) excludes a direct stack and requires at least one unpaired nucleotide.

### 9.2 Confirmed remaining defects

These findings affect both Legacy and Next unless marked otherwise. They were audited but not modified in this task.

| Status | Finding | Consequence |
|---|---|---|
| reproduced, source-confirmed | Heuristic noLP counts the direct stacked continuation and also admits the equivalent w1=w2=1 loop continuation, violating Equation (16). | Base-pair examples: length 4 exact Eall=-5.30 vs heuristic -5.46; length 5 exact -6.50 vs heuristic -6.73. The heuristic improperly has more weight than the exact recurrence. |
| reproduced, source-confirmed | Ensemble updateZ bypasses the base predictor's noGU-end and maxED filters and directly increments Zall. | One-base G/U reports Eall=-1 with noGU-end both false and true, although the true ensemble is empty; unchecked arithmetic can also overflow or propagate non-finite weights. |
| source-confirmed | Base-pair \(E_S\) uses full \(Q\), although its API defines substructures with at least one intramolecular pair. | Correct substructure energy is \(E_S=-RT\log(Q-1)\); whole-monomer energy remains \(-RT\log Q\). |
| reproduced, source-confirmed | PredictorMfeEns retains an unordered map keyed by all four site boundaries while the exact 2-D matrix is advertised as \(O(n_1n_2)\). | Retained space is \(O(n_1n_2L_1L_2)\), quartic at full length. Complementary toys used 9.6 MiB at n=20 and 106.2 MiB at n=50. |
| source-confirmed | Some seed predictors return before initZ on reuse. | Z/site state can leak from an earlier prediction. |
| source-confirmed | Out-of-range Nussinov getQb returns one rather than zero. | An invalid paired state contributes multiplicative identity. |
| source-confirmed, Legacy only | Seed extension adds independent Boltzmann factors instead of multiplying them; also duplicates adjacent noLP stacks and contains unsigned-overlap, one-past-end, and energy-vs-partition comparison faults. | Incorrect seed ensemble weights and boundary behavior. These seed defects were corrected in Next. |

For exact no-seed Legacy/Next prediction, each four-boundary site is completed once. A proposed optimization is to stream its energy directly into the top-k and Zall accumulators rather than retain the entire site map. This should restore the advertised \(O(n_1n_2)\) DP-space scale, but it has not been implemented here.

## 10. Performance evidence

### 10.1 All-interval accessibility

A 120-nt base-pair all-interval tPu workload emits 7,260 nonempty interval values.

| Version | Wall time | Peak RSS |
|---|---:|---:|
| repeated constrained folds | 4.35 s | 5,332 KiB |
| batched inside/outside recurrence | 0.01 s | 5,888 KiB |

This is over 400x faster, with a small table-related memory increase. Both outputs have SHA-256

~~~text
f81bd2eae8f08476e6cbd6a749227926914af630b04003b03a965c6bba0be864
~~~

The theoretical result is the more portable claim: all intervals changed from \(O(n^5)\) time to \(O(n^3)\) time and \(O(n^2)\) memory.

### 10.2 Interaction predictor

The durable biological matrix uses the public fhlA/OxyS, ilvE/GcvB-ST, and phoB/GcvB sequences with energy B, accessibility N, no seed, exact mode M, model S, interaction length at most 20, top three energy-sorted CSV rows, and one thread. The benchmark performs two warm-ups and nine alternating timed subprocess runs. It rejects a case before timing unless Legacy and IntaRNAnew stdout are byte-for-byte identical. Measurements below are medians on GCC 13.3.0 and an AMD Ryzen 5 7530U.

Before the redesign, this workload exposed an architectural mismatch:

- every DP state owned and repeatedly copied a full interaction path, seed vector, and log weight;
- every grid cell retained start-indexed states in vectors and rebuilt an unordered map while propagating them;
- model S paid for millions of log-sum-exp operations even though its representative ignores accumulated path weights;
- every terminal state was reevaluated and materialized into a map-held Interaction, followed by a global sort to return three rows;
- disabled accessibility eagerly constructed whole-monomer partition summaries even though the selected columns did not use them.

For fhlA/OxyS this caused 14.50 million state merges, 14.02 million log additions, 30.36 million heap allocations, 3.505 GB of cumulative heap traffic, and 359.8 MB peak live heap. The two larger cases each retained more than four million PathStates and materialized more than three million sites, producing about 2.4 GiB RSS.

For a fixed pairable start, write \(C[i,j]\) for the maximum number of pairs in a path ending at target offset \(i\) and reverse-query offset \(j\). In the base-pair energy model,

\[
C[0,0]=1,\qquad
C[i,j]=1+\max_{\substack{i-d_T\le a<i\\j-d_Q\le b<j}} C[a,b]
\]

at pairable endpoints, with zero denoting an unreachable state and \(d_T,d_Q\) equal to the configured loop reach plus one. The interaction energy is exactly \(-C[i,j]+E_{\mathrm{add}}\). Thus this selected-only workload needs neither a path nor a Boltzmann weight during propagation.

The implemented execution path:

- derives explicit requirements for complete sites, interaction partitions, monomer partitions, and traceback from the output plan;
- uses the scalar max-plus recurrence for exact, unseeded, unconstrained base-pair/model-S prediction;
- stores only \(\min(L_T,d_T)\) target rows plus query-side sliding maxima;
- streams candidates into the bounded top-k in the same energy/coordinate order as the generic predictor;
- constructs only the selected Interactions and skips unrequested monomer folding.

Direct library configurations remain conservative. The compact path also validates the actual accessibility provider objects and falls back for constraints, multiple regions/windows, structural output, partitions, seeds, noLP/noGU-end, and unsupported models.

Final parity-gated results:

| Case | Legacy median | IntaRNAnew median | Speedup | Legacy RSS | IntaRNAnew RSS |
|---|---:|---:|---:|---:|---:|
| fhlA/OxyS | 504.936 ms | 36.066 ms | 14.000x | 7,424 KiB | 5,120 KiB |
| ilvE/GcvB-ST | 4008.653 ms | 208.473 ms | 19.229x | 7,364 KiB | 5,120 KiB |
| phoB/GcvB | 3821.776 ms | 205.155 ms | 18.629x | 7,364 KiB | 5,120 KiB |
| geometric mean | | | **17.117x** | | |

Relative to the pre-redesign IntaRNAnew, final wall time improved by 94.526x, 144.354x, and 144.117x respectively, a 125.285x geometric mean. Peak RSS fell from 386,548/2,524,132/2,436,400 KiB to 5,120 KiB, reductions of 75.5x, 493.0x, and 475.9x.

Final allocation profiles are also independent of the millions of represented sites:

| Case | Allocation calls | Cumulative heap | Peak live heap |
|---|---:|---:|---:|
| fhlA/OxyS | 862 | 182,037 B | 119,753 B |
| ilvE/GcvB-ST | 847 | 183,552 B | 119,917 B |
| phoB/GcvB | 847 | 183,043 B | 119,707 B |

The generic path remains necessary for richer scientific outputs. Its PathState/hash representation is still the main performance target outside the compact kernel's deliberately narrow, output-observable contract.

## 11. Validation

The final scientific and ownership changes were checked with:

- the complete native make check suite;
- ASan and UBSan over the complete suite;
- all 17 executable-level Legacy fixtures byte-for-byte;
- independent exact predictor enumeration;
- 300 randomized compact-versus-full predictor comparisons, including unbounded/oversized length limits and additive energies;
- focused compact fallback tests for blocked providers, multiple regions, huge output counts, and floating-energy tie collapse;
- explicit Legacy parity checks for every non-structural compact-path CSV field, monomer/interaction partitions, and ensemble output;
- exhaustive native accessibility checks over all intervals and hard constraints;
- 2,400 randomized noncrossing pair graphs for n=1..8, totaling 21,155 interval-partition comparisons;
- RNAeval/Vienna forced-structure and duplex int21 oracles;
- Turner99 and Andronescu07 named-parameter partition oracles;
- exact model-P terminal and dangle tests;
- exact seeded and seed-only tests after destructive state moves;
- 50,016 edge/random comparisons proving the allocation-free comparator preserves the old decimal-string ordering;
- byte parity for the biological performance workloads.

Reproduction from IntaRNAnew:

~~~sh
INTARNANEW_PARAMETER_DIR=../.conda-env/share/ViennaRNA \
  make check

make sanitize

INTARNA_LEGACY_BIN=../src/bin/IntaRNA \
LD_LIBRARY_PATH=../.conda-env/lib \
  make legacy-check
~~~

## 12. Remaining work

Priorities after this audit are:

1. generalize output-aware compact prediction to nearest-neighbor energy, supported constraints, and on-demand traceback without weakening exact semantics;
2. replace the remaining generic PathState vectors and unordered maps with flat states and persistent backpointers;
3. skip unused log-weight accumulation for models S and X in generic configurations and make valid-path range access constant time;
4. derive a state-indexed outside recurrence for Turner/noLP/noGU interval probabilities instead of per-interval constrained folding;
5. make approximate model-P output explicit in heuristic mode, or make IntaRNAens default to exact mode;
6. implement and test the Legacy/Next fixes in Section 9, especially noLP double counting, updateZ filtering, \(Q-1\) for \(E_S\), and streamed exact-site accumulation.

The central results are now implemented and tested: scalar all-interval accessibility is an exact cubic computation, IntaRNAnew's model-P weights include every energy term exactly once, and selected-only exact base-pair prediction is faster than mature Legacy on every parity-gated biological workload.
