
# auxiliary tools

## [IntaRNA-fuse.py](IntaRNA-fuse.py)

Computes non-overlapping interaction sites and tests whether each suboptimal
site can coexist with the minimum-free-energy (MFE) interaction. It reruns
IntaRNA with the suboptimal site blocked on both strands. If the resulting MFE
overlaps the original MFE by at least two nucleotides on **each** strand, it
reports a fused interaction using the suboptimal site and that constrained MFE.
This preserves the compatibility criterion of the former Perl script.

Requires Python 3.6 or newer (standard library only), IntaRNA, and ViennaRNA's
`RNAfold`. `make install` installs `IntaRNA-fuse.py` alongside `IntaRNA`; no pip
installation or `--enable-intarnapvalue` is needed. The script looks for
executables beside itself and on `PATH`, and also checks `../src/bin/IntaRNA`
when run from the source tree. For an out-of-tree build, add its `src/bin`
directory to `PATH`. Make `RNAfold` available on `PATH` or beside the script.

```bash
# Installed script; -t is the target (seq1), -q is the query (seq2).
IntaRNA-fuse.py -t GGGGGGGAAAAAGGGGGGG -q CCCCCCCAAAAACCCCCCC -n 3

# From a built source checkout:
python3 python/IntaRNA-fuse.py -t GGGGGGGAAAAAGGGGGGG -q CCCCCCCAAAAACCCCCCC
```

The command replaces `perl/IntaRNA-fuse.pl`. Both `-t` and `-q` are required
literal sequences (A/C/G/U/T/N, case insensitive); FASTA files and multiple
sequences are not supported. Unlike the Perl prototype, supplied sequences are
used as given instead of being overwritten with an example pair.

| Option | Meaning | Default |
| --- | --- | --- |
| `-n` | Maximum number of original predictions, including the MFE; positive integer | `10` |
| `-l` | Maximum intramolecular base-pair span for both IntaRNA and RNAfold; positive integer | `100` |
| `-o` | Comma-separated CSV column names | `E,start1,end1,start2,end2,ED1,ED2,E_init,hybridDB` |
| `-h` | Show help without requiring either executable | |

`-o` must include `E,start1,end1,start2,end2,ED1,ED2,E_init`, in any order.
Additional IntaRNA CSV columns are passed through. Output is semicolon-separated
CSV on stdout, beginning with the original MFE and then each suboptimal row.
A compatible suboptimal row is immediately followed by its fused row. In that
row, `E`, `ED1`, and `ED2` describe the combined interaction; all other fields
contain `suboptimal:constrained-MFE` values. Coordinates are 1-based and
inclusive, in the original orientation of each input sequence. Fused rows are
additional to the `-n` predictions and are not re-sorted by energy. No original
predictions produces just the CSV header; no constrained prediction means no
fused row for that candidate.

For each strand, RNAfold computes its ensemble free energy with both sites
forced unpaired. The difference from the unconstrained ensemble gives the joint
accessibility cost `ED1` or `ED2`. The fused energy (all terms in kcal/mol) is:

```text
E_fused = (E_suboptimal - ED1_suboptimal - ED2_suboptimal)
        + (E_constrained_MFE - ED1_constrained_MFE - ED2_constrained_MFE)
        + ED1_joint + ED2_joint
```

As in the Perl implementation, this retains both interaction initiation terms;
it does not subtract an `E_init` term. Single-strand ensemble energies use
RNAfold's printed precision. This pairwise compatibility check does not search
all combinations of sites or optimize a full multi-site structure. Invalid
input, missing executables, failed commands, or malformed energy/CSV output
produce a nonzero exit status and a diagnostic on stderr, with no result CSV.

## [IntaRNApvalue](intarnapvalue)

A tool for calculating p-values of IntaRNA energy scores.

## [CopomuS](copomus)

The Compensatory Point Mutation Selector helps to identify the most promising
base pairs from an IntaRNA mfe prediction to verify the interaction via
compensatory mutation experiments.


