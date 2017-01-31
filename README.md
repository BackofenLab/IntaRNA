[![Build Status](https://travis-ci.org/BackofenLab/IntaRNA.svg?branch=master)](https://travis-ci.org/BackofenLab/IntaRNA)

# IntaRNA

**Efficient RNA-RNA interaction prediction incorporating accessibility and 
seeding of interaction sites**

During the last few years, several new small regulatory RNAs 
(sRNAs) have been discovered in bacteria. Most of them act as post-transcriptional 
regulators by base pairing to a target mRNA, causing translational repression 
or activation, or mRNA degradation. Numerous sRNAs have already been identified, 
but the number of experimentally verified targets is considerably lower. 
Consequently, computational target prediction is in great demand. Many existing 
target prediction programs neglect the accessibility of target sites and the 
existence of a seed, while other approaches are either specialized to certain 
types of RNAs or too slow for genome-wide searches.

IntaRNA, developed by
[Prof. Backofen's bioinformatics group at Freiburg University](http://www.bioinf.uni-freiburg.de),
is a general and fast approach to the 
prediction of RNA-RNA interactions incorporating both the accessibility of 
interacting sites 
as well as the existence of a user-definable seed interaction. We successfully applied 
IntaRNA to the prediction of bacterial sRNA targets and determined the exact 
locations of the interactions with a higher accuracy than competing programs. 

For testing or ad hoc use of IntaRNA, you can use its webinterface at the

==> **[Freiburg RNA tools IntaRNA webserver](http://rna.informatik.uni-freiburg.de/IntaRNA/)** <==


## Contribution

Feel free to contribute to this project by writing 
[Issues](https://github.com/BackofenLab/IntaRNA/issues) 
with feature requests or bug reports.

## Cite
If you use IntaRNA, please cite our 
[article](http://bioinformatics.oxfordjournals.org/content/24/24/2849):
```
doi: 10.1093/bioinformatics/btn544
```

<br /><br /><br /><br />
<a name="doc" />
# Documentation

## Overview

The following topics are covered by this documentation:

- [Installation](#install)
  - [Dependencies](#deps)
- [Usage and Parameters](#usage)
  - [Prediction modes, their features and emulated tools](#predModes)
  - [Accessibility and unpaired probabilities](#accessibility)
    - [Local versus global unpaired probabilities](#accLocalGlobal)
    - [Read/write accessibility from/to file or stream](#accFromFile)



<br /><br /><br /><br />
<a name="install" />
# Installation

<br /><br />
<a name="deps" />
## Dependencies

- compiler supporting C++11 standard and openmp
- GNU autotools (automake, autoconf, ..)
- [boost C++ library](http://www.boost.org/) version >= 1.50.0
- [Vienna RNA package](http://www.tbi.univie.ac.at/RNA/) version >= 2.3.0




<br /><br /><br /><br />
<a name="usage" />
# Usage and parameters


<br /><br />
<a name="predModes" />
## Prediction modes, their features and emulated tools

For the prediction of *minimum free energy interactions*, the following modes
and according features are supported and can be set via the `--mode` parameter.
The tiem and space complexities are given for the prediction of two sequences
of equal length *n*.

| Features | Heuristic `--mode=0` | Exact-SE `--mode=1` | Exact `--mode=2` |
| -------- | :------------------: | :-----------------: | :--------------: |
| Time complexity | O(*n*^2) | O(*n*^4) | O(*n*^4) |
| Space complexity | O(*n*^2) | O(*n*^2) | O(*n*^4) |
| Seed constraint | x | x | x |
| No seed constraint | x | x | x |
| Minimum free energy interaction | not guaranteed | x | x |
| Overlapping suboptimal interactions | x | x | x |
| Non-overlapping suboptimal interactions | x | - | x |

Given these features, we can emulate and extend a couple of RNA-RNA interaction
tools using IntaRNA.

**TargetScan** and **RNAhybrid** are approaches that predict the interaction hybrid with 
minimal interaction energy without consideratio whether or not the interacting 
subsequences are probably involved involved in intramolecular base pairings. Furthermore,
no seed constraint is taken into account.
This prediction result can be emulated (depending on the used prediction mode) 
by running IntaRNA when disabling both the seed constraint
as well as the accessibility integration using
```bash
# prediction results similar to TargetScan/RNAhybrid
IntaRNA [..] --noSeed --qAcc=N --tAcc=N
```

**RNAup** was one of the first RNA-RNA interaction prediction approaches that took the 
accessibility of the interacting subsequences into account while not considering the seed feature. 
IntaRNA's exact prediction mode is eventually an alternative implementation when disabling
seed constraint incorporation. Furthermore, the unpaired probabilities used by RNAup to score
the accessibility of subregions are covering the respective overall structural ensemble for each
interacting RNA, such that we have to disable accessibility computation based on local folding (RNAplfold)
using
```bash
# prediction results similar to RNAup
IntaRNA --mode=1 --noSeed --qAccW=0 --qAccL=0 --tAccW=0 --tAccL=0
```


<br /><br />
<a name="accessibility" />
## Accessibility and unpaired probabilities

Accessibility describes the availability of an RNA subsequence for intermolecular
base pairing. It can be expressed in terms of the probability of the subsequence
to be unpaired (its *unpaired probability* *Pu*).

A limited accessibility, i.e. a low unpaired probability, can be incorporated into
the RNA-RNA interaction prediction by adding according energy penalties. 
These so called *ED* values are transformed unpaired probabilities, i.e. the
penalty for a subsequence partaking in an interaction is given by *ED=-RT log(Pu)*, 
where *Pu* denotes the unpaired probability of the subsequence. Within the 
IntaRNA energy model, *ED* values for both interacting subsequences are considered.

Accessibility incorporation can be disabled for query or target sequences using
`--qAcc=N` or `--tAcc=N`, respectively.

A setup of `--qAcc=C` or `--tAcc=C` (default) enables accessibility computation 
using the Vienna RNA package routines for query or target sequences, respectively.


<a name="accLocalGlobal" />
### Local versus global unpaired probabilities

Exact computation of unpaired probabilities (*Pu* terms) is considers all possible
structures the sequence can adopt (the whole structure ensemble). This is referred
to as *global unpaired probabilities* as computed e.g. by **RNAup**.

Since global probability computation is (a) computationally demanding and (b) not
reasonable for long sequences, local RNA folding was suggested, which also enables
according *local unpaired probability* computation, as e.g. done by **RNAplfold**.
Here, a folding window of a defined length 'screens' along the RNA and computes
unpaired probabilities within the window (while only intramolecular base pairs 
within the window are considered).

IntaRNA enables both global as well as local unpaired probability computation.
To this end, the sliding window length has to be specified in order to enable/disable
local folding.

#### Use case examples global/local unpaired probability computation
The use of global or local accessibilities can be defined independently 
for query and target sequences using `--qAccW|L` and `--tAccW|L`, respectively.
Here, `--?AccW` defines the sliding window length (0 sets it to the whole sequence length)
and `--?AccL` defines the maximal length of considered intramolecular base pairs,
i.e. the maximal number of positions enclosed by a base pair
(0 sets it to the whole sequence length). Both can be defined
independently while respecting `AccL <= AccW`.
```bash
# using global accessibilities for query and target
IntaRNA [..] --qAccW=0 --qAccL=0 --tAccW=0 --qAccL=0
# using local accessibilities for target and global for query
IntaRNA [..] --qAccW=0 --qAccL=0 --tAccW=150 --qAccL=100
```


<a name="accFromFile" />
### Read/write accessibility from/to file or stream

It is possible to read precomputed accessibility values from file or stream to
avoid their runtime demanding computation. To this end, we support the following
formats

| Input format | produced by |
| ---- | --- |
| RNAplfold unpaired probabilities | `RNAplfold -u` or `IntaRNA --outPuFile*` |
| RNAplfold-styled ED values | `IntaRNA --outAccFile*` |

The **RNAplfold** format is a table encoding of a banded upper triangular matrix 
with band width l. First row contains a header comment on the data starting with
`#`. Second line encodes the column headers, i.e. the window width per column.
Every successive line starts with the index (starting from 1) of the window end
followed by a tabulator separated list for each windows value in increasing
window length order. That is, column 2 holds values for window length 1, column 
3 for length 2, ... . The following provides a short output/input 
example for a sequence of length 5 with a maximal window length of 3.

```
#unpaired probabilities
 #i$	l=1	2	3	
1	0.9949492	NA	NA	
2	0.9949079	0.9941056	NA	
3	0.9554214	0.9518663	0.9511048		
4	0.9165814	0.9122866	0.9090283		
5	0.998999	0.915609	0.9117766		
6	0.8549929	0.8541667	0.8448852		

```

#### Use case examples for read/write accessibilities and unpaired probabilities
If you have precomputed data, e.g. the file `plfold_lunp` with unpaired probabilities
computed by **RNAplfold**, you can run
```bash
# fill accessibilities from RNAplfold unpaired probabilities
IntaRNA [..] --qAcc=P --qAccFile=plfold_lunp

# fill accessibilities from RNAplfold unpaired probabilities via pipe
cat plfold_lunp | IntaRNA [..] --qAcc=P --qAccFile=STDIN
```
Another option is to store the accessibility data computed by IntaRNA for 
successive calls using 
```bash
# storing and reusing (target) accessibility data for successive IntaRNA calls
IntaRNA [..] --outPuFilet=intarna.target.pu
IntaRNA [..] --tAcc=P --tAccFile=intarna.target.pu

# piping (target) accessibilities between IntaRNA calls
IntaRNA [..] --outPuFilet=STDOUT | IntaRNA [..] --tAcc=P --tAccFile=STDIN
```
