# ED cache round-trip review

The first cached-ED experiment has a real parameter confound: computed accessibility allows query/target interaction extents100/150, whereas reloading the saved files advertises99/149. Identical reported best-interaction coordinates and energy on this input do not make those search spaces equivalent.

## Source cause

- `Accessibility.cpp:68-82` writes length labels1 through `min(sequence length,getMaxLength())` and emits values through that maximum where the sequence prefix permits them. It does not write an extra length beyond the configured maximum.
- `AccessibilityFromStream.cpp:88-91` reads the final header length and unconditionally subtracts1, with the comment "max length == 1 smaller than max-accessibility-data-length due to dangling-end treatment".
- At `AccessibilityFromStream.cpp:92-100`, if this inferred value is below the requested cap, `availMaxLength` is reduced and the INFO message reports the change.
- `AccessibilityFromStream.h:158-165` returns `availMaxLength` as the effective maximum. Its `getED` at `:134-153` rejects regions longer than that effective maximum.
- The parser allocates and reads an extra length (`1+getMaxLength()` at `.cpp:105,139`). Thus an extra column can be parsed internally but still cannot be returned as a valid ED through `getED` when it exceeds the reduced advertised maximum.

This is a writer/reader limit asymmetry in the profiled source, not file truncation. The parser comment documents its intended dangling-end reserve; the actual downstream extent reduction must be included in profiling controls. No application fix was attempted as part of this profiling review.

## Inspection of the actual saved files

| File | Final header length | Data rows | Fields per row, including row index | Numeric entries in final length column |
|---|---:|---:|---:|---:|
| `inputs/cache/target.ed` |150|1000|151|851|
| `inputs/cache/query.ed` |100|100|101|1|

The target file's first valid length150 entry (row150) is35.76kcal/mol; its last row's length150 entry is34.58. The query row100 has a valid length100 entry28.47. Earlier impossible length columns contain `NA`, as expected. Therefore neither saved file is missing its longest advertised column.

## Correct matched experiment

Keep the original computed150/cached150-requested records and mark their different effective limits as a confound. Add a new paired experiment with these explicit flags on **both** variants:

```
--qIntLenMax=99 --tIntLenMax=149
```

For the computed variant, retain the original folding settings (`q/tAccW=150`, `q/tAccL=100`, energy parameters, sequences and constraints). For the cached variant, reuse the existing100/150-column files with `--qAcc=E --tAcc=E` and their exact paths. The reader's last-column-minus-one limits then equal the requested99/149 limits, so it need not change them.

Remove `--intLenMax=150` from both commands: the CLI rejects conflicting global and sequence-specific limits. Do not regenerate the cache at99/149 and then request99/149, since the same reader behavior would lower it again to98/148. Keep the existing files immutable and hashed.

The computed run may calculate fewer requested ED lengths while using the same folding ensemble; that is appropriate for comparing these explicitly chosen99/149 interaction limits. Preserve the original cache-generation time separately as the preparation cost for the reusable files. Require matching reported coordinates/energies and, if the later untimed guard includes it, matching base-pair lists; record any differences rather than assuming universal round-trip equivalence.

A matching150-target-cap experiment would require at least151 target columns. The100-nt query cannot provide a valid length101 column, so retaining its full100-nt advertised limit cannot be achieved through the current reader's header-minus-one logic with a normal self-generated cache. The99/149 control is the practical bounded comparison without changing application code or fabricating data.
